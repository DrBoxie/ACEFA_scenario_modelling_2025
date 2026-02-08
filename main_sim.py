import numpy as np
import os
from operator import itemgetter
import main_stochastic as stoc
# import main_determinstic as deter

def main(params, states, incidences, curr_dir, rng, outputfigs = False, det_run = 0, stoc_run = 1, new_phis = False, laiv_ages = []):
    # TODO: Check whether the variables new_phis and laiv_ages are still in use. If they are, it should be via either the laiv_phis
    # calculation, or in the forward projection
    if det_run:
        # states = create_states(params)
        run_det(params, curr_dir)
        
    if stoc_run:
        sols_stoc, R0_series, inc, vacc_post_inf = run_stoc(params, states, incidences, curr_dir, rng, outputfigs, new_phis = new_phis, laiv_ages = laiv_ages)
            
    return sols_stoc, R0_series, inc, vacc_post_inf

def init_params(nr_runs = 1):
    
    tau = 1                                                                     # tau-leaping step size. (Has no effect if stochastic isn't run)
    
    t_start, t_end = 0, 365
    t_iiv_start, t_iiv_end = 59, 226                                            # Start and end date of IIV program
    t_laiv_start, t_laiv_end = 59, 100                                          # Start and end date of LAIV program
    nr_days = t_end - t_start
    
    # age_groups = ["<1", "1-4", "5-11", "12-17", "18-64", "65-79", "80+"]        # 7 age groups
    age_groups = ["<1", "1" ,"2-4", "5-11", "12-17", "18-64", "65-79", "80+"]        # 8 age groups
    compartments = ["S0", "S1", "V0", "V1", "I_S0", "I_S1", "I_V0", "I_V1", "R_S", "R_V"]
    S_comps = ["S0", "S1", "V0", "V1"]                                          # Susceptible compartments
    I_comps = ["I_S0", "I_S1", "I_V0", "I_V1"]                                  # Infectious compartments
    R_comps = ["R_S", "R_V"]
    Unvacc_comps = ["S0", "S1", "I_S0", "I_S1", "R_S"]                          # Unvaccinated compartments
    Vacc_comps = ["V0", "V1", "I_V0", "I_V1", "R_V"]                            # Vaccinated compartments
    Unvacc_I_comps = ["I_S0", "I_S1"]                                           # Unvaccinated infectious compartments
    Vacc_I_comps = ["I_V0", "I_V1"]                                             # Vaccinated infectious compartments
    
    Unvacc_to_Vacc = dict(zip(Unvacc_comps, Vacc_comps))                        # Link each unvaccinated compartment to its vaccinated counterpart
    nr_age = len(age_groups)
    nr_comps = len(compartments)    
    
    tot_pop = 1000000
    
    # Age distribution with 7 age groups
    # frac_ages = np.array([0.010752583, 0.044639546, 0.08338837, 0.073682919, 0.614286160, 0.129143011, 0.044107411])
    
    # Age distribution with 8 age groups
    frac_ages = np.array([0.010752583, 0.010748152, 0.033891394, 0.08338837, 0.073682919, 0.614286159, 0.129143011, 0.044107411])

    # N holds the population for each of the age groups.    
    N = allocate_rounded(tot_pop, frac_ages)
    
    # Contact matrix with 7 age groups (aggregating the 1-2 and the 2-5 into one single group)
    # C = np.array([
    #     [0.375779488, 0.14313033,	0.052502325,	0.02248264,	0.093521173,	0.007517928,	0.001580057],
    #     [0.59290092,	2.726163835,	0.798894201,	0.181099599,	0.428477814,	0.133352435,	0.053152569],
    #     [0.405668348,	1.490154849,	8.467514786,	1.801643719,	0.677721934,	0.347731176,	0.189136569],
    #     [0.154319918,	0.300083216,	1.600482745,	9.463746013,	0.878385651,	0.189544108,	0.126039369],
    #     [5.355102935,	5.922918595,	5.022464015,	7.327712821,	10.38538854,	3.609060442,	2.608662633],
    #     [0.091321412,	0.391043379,	0.546670304,	0.3354363, 0.765615499,	1.753038513, 0.728460955],
    #     [0.006196919,	0.050324215,	0.096003292,	0.072016955,	0.178674775,	0.235198769,	0.129056892]
    # ])
    
    # Contact matrix with 8 age groups
    C = np.array([
        [0.375779488,	0.221425972,	0.117743201,	0.052502325,	0.02248264,	0.093521173,	0.007517928,	0.001580057],
        [0.224587717,	0.68020258,	0.275817706,	0.091685117,	0.032435929,	0.103449184,	0.015630453,	0.00444946],
        [0.368313203,	0.850640683,	2.837925261,	0.707209083,	0.14866367,	0.325028629,	0.117721982,	0.048703109],
        [0.405668348,	0.698447857,	1.746863481,	8.467514786,	1.801643719,	0.677721934,	0.347731176,	0.189136569],
        [0.154319918,	0.219504573,	0.326210601,	1.600482745,	9.463746013,	0.878385651,	0.189544108,	0.126039369],
        [5.355102935,	5.840197141,	5.949740779,	5.022464015,	7.327712821,	10.38538854,	3.609060442,	2.608662633],
        [0.091321412,	0.187192551,	0.457141402,	0.546670304,	0.3354363,	0.765615499,	1.753038513,	0.728460955],
        [0.006196919,	0.017204937,	0.061063042,	0.096003292,	0.072016955,	0.178674775,	0.235198769,	0.129056892]
    ])

     # Uncorrelated contact matrix for testing purposes
     # C = np.array([
     #     [5, 0, 0, 0, 0, 0, 0],
     #     [0, 5, 0, 0, 0, 0, 0],
     #     [0, 0, 5, 0, 0, 0, 0],
     #     [0, 0, 0, 5, 0, 0, 0],
     #     [0, 0, 0, 0, 5, 0, 0],
     #     [0, 0, 0, 0, 0, 5, 0],
     #     [0, 0, 0, 0, 0, 0, 5]
     # ])
    
    frac_S1 = np.full(nr_age, 0.1)
    
    # Everyone in S0 for testing purposes
    # frac_S1 = np.full(nr_age, 0)                                              # Fraction of population in each age group starting in S1

    # =============================================================================
    #   The equaiton for beta is: 
    #   beta(t) = beta_0 * [ 1 + beta_1 sin ( 2 pi (t - shift)/ 365)]
    #   where beta_1 controls the seasonality effect.
    # =============================================================================
    
    # The values of beta_0, beta_1, and shift below are dummy values which are
    # overwritten when running the ABC or the forward projection code
    beta_0 = 0.047
    beta_1= 0.3
    shift = 100
    
    t = np.linspace(t_start, t_end, num = nr_days)
    beta = beta_0 * (1 +  beta_1 * np.sin(( 2 * np.pi * (t - shift) ) / nr_days))
    
    # Sanity check on beta
    if np.any(beta < 0 ):
        raise ValueError("Error: Beta is not allowed to be negative!")
    
    base_phi = np.full(nr_age, 1)
    phi_S1 = np.full(nr_age, 0.7)
    phi_V0 = np.full(nr_age, 0.5)
    phi_V1 = np.full(nr_age, 0.3)
    
    gamma_i = np.full(nr_age, 1/3)

    target_cover = np.array([0.15, 0.15, 0.15, 0.10, 0.10, 0.20, 0.55, 0.55])         # target vaccination coverage in each age group
    
    iiv_rates = np.array([
        0.00266946, 0.00334557, 0.00737364, 0.03162889, 0.07259051,0.1482384,
        0.24145878, 0.36451686, 0.47168377, 0.56855491, 0.65137032, 0.75339126,
        0.81009063, 0.84931351, 0.88319599, 0.91812563, 0.93545698, 0.95122565,
        0.96457128, 0.97785923, 0.9846138 , 0.99014375, 0.99475859, 1.0
        ])

    laiv_rates = np.array([0.16666667, 0.33333333, 0.5, 0.66666667, 0.83333333, 1])    
    
    daily_vacc_admin = np.zeros((nr_days, nr_age))
    
    # Calculate the number of vaccinations administers on each day    
    daily_vacc_admin = allocate_vacc(daily_vacc_admin, target_cover, N, iiv_rates, t_iiv_start, t_iiv_end)

    p_ext = 0.0005                                                             # Percentage of the total population to be randomly infected from the outside throughout the whole run
    
    # Turn off external seeding of infection for testing purposes
    # p_ext = 0                                                  
    
    params = {"nr_runs": nr_runs,
              "tau": tau,
              "t_start": t_start,
              "t_end": t_end,
              "nr_days": nr_days,
              "t_iiv_start": t_iiv_start,
              "t_iiv_end": t_iiv_end,
              "t_laiv_start": t_laiv_start,
              "t_laiv_end": t_laiv_end,
              "age_groups": age_groups,
              "compartments": compartments,
              "S_comps": S_comps,
              "I_comps": I_comps,
              "R_comps": R_comps,
              "Unvacc_I_comps": Unvacc_I_comps,
              "Vacc_I_comps": Vacc_I_comps,
              "Unvacc_comps": Unvacc_comps,
              "Vacc_comps": Vacc_comps,
              "Unvacc_to_Vacc": Unvacc_to_Vacc,
              "nr_age": nr_age,
              "nr_comps": nr_comps,
              "tot_pop": tot_pop,
              "frac_ages": frac_ages,
              "frac_S1": frac_S1,
              "N": N,
              "C": C,
              "beta": beta,
              "base_phi": base_phi,
              "phi_S1": phi_S1,
              "phi_V0": phi_V0,
              "phi_V1": phi_V1,
              "gamma": gamma_i,
              "daily_vacc_admin": daily_vacc_admin,
              "p_ext": p_ext,
              "shift": shift,
              "laiv_rates": laiv_rates,
              "target_cover": target_cover
              }
    
    return params

def create_states(params):
    compartments, nr_age, nr_days = itemgetter("compartments", "nr_age", "nr_days")(params)
    states = {key: np.zeros((nr_age, nr_days)) for key in compartments}
    return states

def create_daily_incidence(params):
    compartments, nr_age, nr_days = itemgetter("compartments", "nr_age", "nr_days")(params)
    skip = {"S0", "S1"}
    daily_inc = {key: np.zeros((nr_age, nr_days)) for key in compartments if key not in skip}
    return daily_inc

def allocate_rounded(pop, frac):
    # Allocate population to group based on fractions, keeping their sum equal.
    # This is done by adding to / removing from the largest age group as needed.
    
    frac, N = np.modf(pop * frac)
    N = N.astype(int)
    
    diff = pop - np.sum(N)
    if diff != 0:
        indices = np.argsort(frac)[-diff:][::-1]
        N[indices] += 1
        
    return N

def allocate_vacc(daily_vacc_admin,target_cover, N, vacc_rates, t_vacc_start, t_vacc_end):
    total_to_admin = target_cover * N
    weekly_perc = np.diff(np.concatenate([[0], vacc_rates]))
    weekly_doses = weekly_perc[:, None] * total_to_admin[None, :]
    daily_doses = np.repeat(weekly_doses / 7, repeats = 7, axis = 0)
    
    duration_rollout = t_vacc_end - t_vacc_start
    
    for day in range(duration_rollout):
        
        # np.modf outputs a tuple where the second element is a vector with the integer part of the number of vaccination to administer ot that day to each age group
        # and the first element is a vector with the fractional part for each age group, which is then moved onward to the number of vaccines for the following day
        frac_admin, integ_admin = np.modf(daily_doses[day])
        
        daily_vacc_admin[t_vacc_start + day] = integ_admin.astype(int)
        daily_doses[day+1] += frac_admin
    
    daily_vacc_admin[t_vacc_start + duration_rollout] += np.round(daily_doses[duration_rollout]).astype(int)

    return daily_vacc_admin

def run_det(params, curr_dir):
    
    # This function is not in use anymore and may not work in the current setting without some tweaks
    
    # deter.main_run(params)
    os.makedirs("figures deterministic", exist_ok=True)
    os.makedirs("csv deterministic", exist_ok=True)

def run_stoc(params, states, incidences, curr_dir, rng, outputfigs, new_phis = False, laiv_ages = []):
    
    sols_stoc, R0_series, incidences, vacc_post_inf = stoc.simul(params, states, incidences, rng, new_phis = new_phis, laiv_ages = laiv_ages)
    
    if outputfigs:
        fig_folder_path = os.path.join(curr_dir, "figures stochastic")
        csv_folder_path = os.path.join(curr_dir, "csv stochastic")
        os.makedirs(fig_folder_path, exist_ok=True)
        os.makedirs(csv_folder_path, exist_ok=True)
        stoc.plotting_output(params, sols_stoc, R0_series, incidences)
    
    return sols_stoc, R0_series, incidences, vacc_post_inf

if __name__ == "__main__":
    
    current_dir = os.getcwd().lower()
    rng = np.random.default_rng()
    
    # Turn stochastic and / or deterministic runs on (1) or off (0)
    # Deterministic code has not been updated and need to be checked if it needs to be run
    det_run = 0
    stoc_run = 1
    outputfigs = False
    params = init_params()
    states = create_states(params)
    daily_inc = create_daily_incidence(params)
    main(params, states, daily_inc, current_dir, rng, outputfigs, det_run, stoc_run)
    