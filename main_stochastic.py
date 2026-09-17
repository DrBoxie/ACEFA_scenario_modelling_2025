import numpy as np
from operator import itemgetter

def initial_conditions(params, states): 
    N, frac_S1 = itemgetter("N", "frac_S1")(params)
    
    # Everyone starts in the susceptible (unvaccinated) states
    S0_init = np.rint(N * (1-frac_S1)).astype(int)
    S1_init = np.rint(N * frac_S1).astype(int)
    
    # Any group whose size became negative due to the rounding above, is
    # floored at zero. This shouldn't happen mathematically, but rounding could cause this to happen (though still very unlikely).
    S0_init[S0_init < 0] = 0 
    S1_init[S1_init < 0] = 0 
    
    # Ensure the total population remains unchanged, by adding / removing
    # people from the largest group in S0 (as it would have the least 
    # proportional impact) to account for potential changes
    # that occured due to the roundings above
    diff = N - S0_init -S1_init
    
    if (diff != 0).any():
        max_age_group = np.argmax(N)
        S0_init[max_age_group] += np.sum(diff)
    
    return S0_init, S1_init

def current_pop(states, params, day): 
    # Extract today's population from the "states" dictionary
    
    compartments = params["compartments"]
    # Current population in each compartment
    return {c: states[c][:, day] for c in compartments}

def administer_vacc(states, incidences, params, day, nr_to_admin, vacc_post_inf, new_phis = False, laiv_ages = []):
   
    Unvacc_comps, Vacc_comps, Unvacc_to_Vacc = itemgetter("Unvacc_comps", "Vacc_comps", "Unvacc_to_Vacc")(params)
    if new_phis:
        laiv_idx = [params["age_groups"].index(age) for age in laiv_ages]
    
    curr = current_pop(states, params, day)                                     # Current population in each compartment
    sizes_unvacc_comps = np.vstack([curr[k] for k in Unvacc_comps])
    weights = sizes_unvacc_comps / sizes_unvacc_comps.sum(axis=0)
    raw = nr_to_admin * weights                                                  # Raw (non-rounded) vaccine allocation
    
    # int_vals contains the numbers of vaccines to administer today, in the format (compartment, age_group) where the compartments
    # appear in the same order as in the variable Unvacc_comps
    frac_part, int_vals = np.modf(raw)
    int_vals = int_vals.astype(int)
    remainder = (nr_to_admin - int_vals.sum(axis=0)[:]).astype(int)
    
    # Distributed the leftover vaccines to the compartments with the largest fractional remainders
    
    for j, r in enumerate(remainder):
        if r > 0:
            inds = np.argsort(-frac_part[:,j])[:r]
            int_vals[inds, j] += 1
    
    for i, c in enumerate(Unvacc_comps):
        states[c][:, day] -= int_vals[i, :]
        states[Unvacc_to_Vacc[c]][:, day] += int_vals[i, :]
        incidences[Unvacc_to_Vacc[c]][:, day] += int_vals[i, :]
        
        # If the compartment is not S0 or S1, than the vaccination was given to an individual who was already infected
        if c[0] != "S":
            # if the run is for calculating new phis for forward projections, only count those that are in the age groups that receive LAIV
            # otherwise count all of them
            if new_phis:
                vacc_post_inf[laiv_idx] += int_vals[i, laiv_idx].astype(int)
            else:
                vacc_post_inf += int_vals[i, :].astype(int)
        
    return vacc_post_inf
   
def tau_leap_step(states, incidences, params, day, nr_sub_int, p_ext, rng):
    """
    Perform one tau-leap sub-interval update for the current day.
    Remarks:
    - The p_ext in the arguments of the function is already the daily one, rather than the yearly one given in the initialisation of the run.
    - The code assumes that the gammas are constant in time and the same across all infectious compartments (though not necessarily across all age groups)
    """
    
    C, N, beta_full, phi_S1, phi_V0, phi_V1, gamma, half_life, compartments, S_comps, I_comps = itemgetter("C", "N", "beta", "phi_S1", "phi_V0", "phi_V1",\
                                                                                       "gamma", "half_life", "compartments", "S_comps", "I_comps" )(params)
    
    Susc0_comps = [c for c in S_comps if "0" in c]
    Susc1_comps = [c for c in S_comps if "1" in c]
    
    beta = beta_full[day]
    dt = 1.0 / nr_sub_int
    
    # Current population in each compartment
    curr = current_pop(states, params, day)

    omega = np.log(2) / half_life[0]
    
    # ---- Force of infection ----
    I_total = np.sum([curr[c] for c in I_comps], axis=0)
    # The 2 lines below are a calculation of lambda that catches cases where some compartments are empty correctly, avoiding div by 0 errors
    with np.errstate(divide='ignore', invalid='ignore'):
        I_div_N = np.divide(I_total, N, out=np.zeros_like(I_total), where=N != 0)
    lambda_i = beta * (C.T @ I_div_N)
    
    # ---- Binomial sampling for infections from S0 and V0 ----    
    n_from_no_prior_immun = np.array([curr[c] for c in Susc0_comps]).astype(int)
    p_inf_no_prior_immun = np.array([
        1 - np.exp(-(lambda_i + p_ext) * dt),
        1 - np.exp(-(lambda_i + p_ext) * phi_V0 * dt)
        ])
    infs_no_prior_immun = rng.binomial(n_from_no_prior_immun, p_inf_no_prior_immun)
    
    # ---- Binomial sampling for infections from S1 and V1 ----    
    n_from_prior_immun = np.array([curr[c] for c in Susc1_comps]).astype(int)
    
    inf_rate = np.array([
        (lambda_i + p_ext) * phi_S1, 
        (lambda_i + p_ext) * phi_V1
        ])
    
    total_rate = inf_rate + omega
    if np.any(total_rate == 0):
        raise ValueError("Error: rates for infections and / or waning cannot be zero!")
        
    p_leave = 1 - np.exp(-total_rate * dt)
    nrs_leave = rng.binomial(n_from_prior_immun, p_leave)
    
    infs_with_prior_immun = rng.binomial(nrs_leave, inf_rate / total_rate)
    nrs_wane = nrs_leave - infs_with_prior_immun
    
    # ---- Binomial sampling for recoveries ----    
    p_rec = 1 - np.exp(-gamma[0] * dt)
    rec = rng.binomial(np.array([curr[c] for c in I_comps]).astype(int), p_rec)

    # ---- Update compartments ----
    curr["S0"], curr["I_S0"] = curr["S0"] - infs_no_prior_immun[0], curr["I_S0"] + infs_no_prior_immun[0]
    curr["V0"], curr["I_V0"] = curr["V0"] - infs_no_prior_immun[1], curr["I_V0"] + infs_no_prior_immun[1]
    
    curr["I_S0"], curr["R_S"] = curr["I_S0"] - rec[0], curr["R_S"] + rec[0]
    curr["I_S1"], curr["R_S"] = curr["I_S1"] - rec[1], curr["R_S"] + rec[1]
    
    curr["I_V0"], curr["R_V"] = curr["I_V0"] - rec[2], curr["R_V"] + rec[2]
    curr["I_V1"], curr["R_V"] = curr["I_V1"] - rec[3], curr["R_V"] + rec[3]
    
    curr["S1"], curr["I_S1"], curr["S0"] = curr["S1"] - infs_with_prior_immun[0] - nrs_wane[0], curr["I_S1"] + infs_with_prior_immun[0], curr["S0"] + nrs_wane[0]
    curr["V1"], curr["I_V1"], curr["V0"] = curr["V1"] - infs_with_prior_immun[1] - nrs_wane[1], curr["I_V1"] + infs_with_prior_immun[1], curr["V0"] + nrs_wane[1]
    
    # Add today's infections and recoveries to the daily incidences overview
    incidences["S0"][:, day] += nrs_wane[0]
    incidences["V0"][:, day] += nrs_wane[1]
    
    incidences["I_S0"][:, day] += infs_no_prior_immun[0]
    incidences["I_V0"][:, day] += infs_no_prior_immun[1]
    
    incidences["I_S1"][:, day] += infs_with_prior_immun[0]
    incidences["I_V1"][:, day] += infs_with_prior_immun[1]
    
    incidences["R_S"][:, day] += rec[0] + rec[1]
    incidences["R_V"][:, day] += rec[2] + rec[3]
    
    # Write back updated values to states at current day
    for c in compartments:
        states[c][:, day] = curr[c]

def compute_R0(params): 
    # This calculated R0 with the initial value of beta during a run
    # and assuming all the gammas are the same (I'm using the value of gamma[0])
    
    C, gamma, beta_0, beta_1, N = itemgetter("C", "gamma", "beta_0", "beta_1", "N")(params)
    
    beta = beta_0
    gamma = gamma[0]
    
    N_a = N[:, np.newaxis]
    N_b = N[np.newaxis, :]
    
    K = (beta / gamma) * C * (N_a / N_b)
    
    # R0 = max(np.linalg.eigvals(K).real)
    R0 = np.max(np.abs(np.linalg.eigvals(K)))
    
    return R0
    
def simul(params, states, incidences, rng, new_phis = False, laiv_ages = []):
        
    tau, p_ext, nr_days, nr_age, compartments, daily_vacc_admin = itemgetter("tau", "p_ext", "nr_days", "nr_age", "compartments", "daily_vacc_admin")(params)
    
    # divide each day into subintervals depending on the size of tau, and for each time step calculate all the transitions that day
    nr_sub_int = int(1/tau)
    # the division of p_ext by nr_days is to get the daily external infection rate
    daily_ext_inf_rate = p_ext / nr_days
    
    sols = []
    R0_series = []
    vacc_post_inf_list = []
    
    # The initialisation of the run takes care of the initial conditions
    states["S0"][:, 0], states["S1"][:, 0] = initial_conditions(params, states)
    
    vacc_post_inf = np.zeros(nr_age, dtype=int)
    
    R0_run = []
    
    # # Calculate R0 for the day
    R0_run.append(compute_R0(params))  
    
    for day in range(nr_days):
        
        # First administer the day's vaccinations at the start of the day

        admin_today = daily_vacc_admin[day].astype(int)
        if np.sum(admin_today) > 0:
            vacc_post_inf = administer_vacc(states, incidences, params, day, admin_today, vacc_post_inf, new_phis = new_phis, laiv_ages = laiv_ages)
        
        # divide each day into 1/tau subintervals, and simulate transitions within it
        for sub_int in range(nr_sub_int):
            
            tau_leap_step(states, incidences, params, day, nr_sub_int, daily_ext_inf_rate, rng)
        
        # at the end of each day (except for the last one), transition the 
        # final populations of each compartments to be the starting 
        # population for the next day
        if day < nr_days-1:
            for c in compartments:
                states[c][:, day+1] = states[c][:, day]
    
    sols.append(states)
    R0_series.append(R0_run)
    vacc_post_inf_list.append(vacc_post_inf)
    
    return sols, R0_series, incidences, vacc_post_inf_list
     
