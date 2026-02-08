import numpy as np
import matplotlib.pyplot as plt
import os
from operator import itemgetter
import pandas as pd

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

def administer_iiv(states, incidences, params, day, nr_to_admin):
    
    # Administer LAIV
    
    Unvacc_comps, Vacc_comps, Unvacc_to_Vacc = itemgetter("Unvacc_comps", "Vacc_comps", "Unvacc_to_Vacc")(params)
    
    curr = current_pop(states, params, day)                                     # Current population in each compartment
    
    Unvacc = {k: curr[k] for k in Unvacc_comps}
    totals = np.sum([arr for arr in Unvacc.values()], axis=0)                   # Total unvaccinated population per age group

    allocation = {comp: arr * nr_to_admin / totals for comp, arr in Unvacc.items()}
    # allocation = {comp: arr * eta / totals for comp, arr in Unvacc.items()}
    
    # stack all arrays row-wise
    vac_not_rounded = np.array([allocation[k] for k in Unvacc_comps])
    vac_rounded = np.array([round_preserve_sum(vac_not_rounded[:, i]) for i in range(vac_not_rounded.shape[1])]).T
    
    for k, (u, v) in enumerate(zip(Unvacc_comps, Vacc_comps)):
        states[u][:, day] -= vac_rounded[k, :]
        states[v][:, day] += vac_rounded[k, :]
    
    for c in Unvacc_comps:
        incidences[Unvacc_to_Vacc[c]][:, day] += vac_rounded[Unvacc_comps.index(c) , :]
   
def tau_leap_step(states, incidences, params, day, nr_sub_int, p_ext, rng):    
    """
    Perform one tau-leap sub-interval update for the current day.
    Remarks:
    - The p_ext in the arguments of the function is already the daily one, rather than the yearly one given in the initialisation of the run.
    - The code assumes that the gammas are constant in time and the same across all infectious compartments (though not necessarily across all age groups)
    """
    
    C, N, beta_full, phi_S1, phi_V0, phi_V1, gamma, compartments, S_comps, I_comps = itemgetter("C", "N", "beta", "phi_S1", "phi_V0", "phi_V1",\
                                                                                       "gamma", "compartments", "S_comps", "I_comps" )(params)
    beta = beta_full[day]
    dt = 1.0 / nr_sub_int
    
    # Current population in each compartment
    curr = current_pop(states, params, day)
    
    # ---- Force of infection ----
    I_total = np.sum([curr[c] for c in I_comps], axis=0)
    # The 2 lines below are a calculation of lambda that catches cases where some compartments are empty correctly, avoidind div by 0 errors
    with np.errstate(divide='ignore', invalid='ignore'):
        I_div_N = np.divide(I_total, N, out=np.zeros_like(I_total), where=N != 0)
    lambda_i = beta * (C.T @ I_div_N)
    
    # ---- Binomial sampling for infections ----    
    infs = rng.binomial(np.array([curr[c] for c in S_comps]).astype(int), 
                        np.minimum(
                            np.vstack([
                                lambda_i + p_ext, 
                                lambda_i * phi_S1 + p_ext, 
                                lambda_i * phi_V0 + p_ext, 
                                lambda_i * phi_V1 + p_ext
                                ]) * dt, 
                            1
                            )
                        )
    
    # ---- Binomial sampling for recoveries ----    
    rec = rng.binomial(np.array([curr[c] for c in I_comps]).astype(int), np.vstack( [np.minimum(gamma * dt, 1)] * len(I_comps)))

    # ---- Update compartments ----
    curr["S0"], curr["I_S0"] = curr["S0"] - infs[0], curr["I_S0"] + infs[0]
    curr["S1"], curr["I_S1"] = curr["S1"] - infs[1], curr["I_S1"] + infs[1]
    curr["V0"], curr["I_V0"] = curr["V0"] - infs[2], curr["I_V0"] + infs[2]
    curr["V1"], curr["I_V1"] = curr["V1"] - infs[3], curr["I_V1"] + infs[3]
    
    curr["I_S0"], curr["R_S"] = curr["I_S0"] - rec[0], curr["R_S"] + rec[0]
    curr["I_S1"], curr["R_S"] = curr["I_S1"] - rec[1], curr["R_S"] + rec[1]
    
    curr["I_V0"], curr["R_V"] = curr["I_V0"] - rec[2], curr["R_V"] + rec[2]
    curr["I_V1"], curr["R_V"] = curr["I_V1"] - rec[3], curr["R_V"] + rec[3]
    
    # Add today's infections and recoveries to the daily incidences overview
    for c in I_comps:
        incidences[c][:, day] += infs[S_comps.index(c[2:]) , :]
    
        # The third character in the I_comps strings is either an S or a V.
        # Based on this, allocate the recovered incidences to the correct 
        # recovered state incidence overview
        if c[2] == "S":
            incidences["R_S"][:, day] += rec[I_comps.index(c) ,:]
        elif c[2] == "V":
            incidences["R_V"][:, day] += rec[I_comps.index(c) ,:]
        else:
            raise ValueError("Error: One of the I states seems to not have an S or a V as its third character in its name!")
    
    # Write back updated values to states at current day
    for c in compartments:
        states[c][:, day] = curr[c]

def round_preserve_sum(vec):
    # code snippet to round while keeping sum of vector fixed, assuming the the original sum of the vector was already an integer
    floored = np.floor(vec).astype(int)
    remainder = int(vec.sum() - floored.sum())  # number of 1s to distribute
    # Find indices with largest fractional parts
    frac_indices = np.argsort(vec - floored)[::-1]
    floored[frac_indices[:remainder]] += 1
    return floored

def compute_R0(params):
    # This calculated R0 with the initial value of beta during a run
    # and assuming all the gammas are the same (I'm using the value of gamma[0])
    # The commented out sections are kept as they could be relevant for the R_eff calculatio
    # should we want to implement it
    
    C, gamma, beta_full, N = itemgetter("C", "gamma", "beta", "N")(params)
    
    beta = beta_full[0]
    gamma = gamma[0]
    
    N_a = N[:, np.newaxis]
    N_b = N[np.newaxis, :]
    
    K = (beta / gamma) * C * (N_a / N_b)
    
    R0 = max(np.linalg.eigvals(K).real)
    
    return R0
    
def plotting_output(params, sols, R0_series, inf): 
    
    t_iiv_start, t_iiv_end, compartments, t_start, t_end, age_groups, nr_runs = itemgetter("t_iiv_start", "t_iiv_end", "compartments", "t_start", "t_end", "age_groups", "nr_runs")(params)

    t_eval = np.arange(t_start, t_end)                                          # full time interval

    labels = [f"Age group {ag}" for ag in age_groups]
    
    colors = plt.cm.Accent.colors                                               # Decent-looking colour scheme for figures

    for run in range(nr_runs):    
        
        fig_subfolder = os.path.join("figures stochastic", f"Run {run}")
        os.makedirs(fig_subfolder, exist_ok=True)
        csv_subfolder = os.path.join("csv stochastic", f"Run {run}")
        os.makedirs(csv_subfolder, exist_ok=True)
        
        states = sols[run]
        R0_run = R0_series[run]
        
        # Total per compartment over time
        for i, label in enumerate(compartments):
            data = states[label]    
            total = data.sum(axis=0)

            # # Get the maximum y-value
            # ymax = np.max(total)
    
            fig_filepath = os.path.join(fig_subfolder, f"{label}_stochastic.png")
            csv_filepath = os.path.join(csv_subfolder, f"{label}_stochastic.csv")
            
            plt.figure(figsize=(10, 6))
            plt.stackplot(t_eval, states[f"{label}"], labels=labels, colors=colors, alpha = 0.8)
            plt.xlim(left = 0)
            plt.axvspan(t_iiv_start, t_iiv_end, color='orange', alpha=0.2, label='Vaccination period')
            plt.xlabel("Day")
            plt.ylabel(f"Population in {label}")
            # plt.title(f"{label} by age group, stochastic run (R0 = {R0_run[0]:.2f})")
            plt.title(f"{label} by age group, stochastic run")                # This line is the same as the one above but without the R0 in case we don't want it in the title
            plt.legend(loc="upper right", frameon=False)
            plt.grid(False)
            plt.tight_layout()
            plt.savefig(fig_filepath)
            plt.close()
            
            # Save CSV
            df = pd.DataFrame(data.T, columns=age_groups)
            df["total"] = total
            df["day"] = t_eval
            df.to_csv(csv_filepath, index=False)

def simul(params, states, incidences, rng, new_phis = False, laiv_ages = []):
    # TODO: Check whether the variables new_phis and laiv_ages are still in use. If they are, it should be via either the laiv_phis
    # calculation, or in the forward projection
    
    # tau, nr_runs, p_ext, nr_days, nr_age, t_iiv_start, t_iiv_end, t_laiv_start, t_laiv_end, compartments, I_comps, N, daily_vacc_admin = itemgetter("tau", 
    #         "nr_runs","p_ext", "nr_days", "nr_age", "t_iiv_start", "t_iiv_end","t_laiv_start", "t_laiv_end", "compartments", "I_comps", 
    #         "N", "daily_vacc_admin")(params)
    
    tau, nr_runs, p_ext, nr_days, nr_age, compartments, daily_vacc_admin = itemgetter("tau", 
            "nr_runs","p_ext", "nr_days", "nr_age", "compartments", "daily_vacc_admin")(params)
    
    # divide each day into subintervals depending on the size of tau, and for each time step calculate all the transitions that day
    nr_sub_int = int(1/tau)
    # the division of p_ext by nr_days is to get the daily external infection rate
    daily_ext_inf_rate = p_ext / nr_days
    
    sols = []
    R0_series = []
    vacc_post_inf_list = []
    
    for run in range(nr_runs):
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
     