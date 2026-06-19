import numpy as np
import os
from operator import itemgetter
import main_stochastic as stoc

def main(params, states, incidences, rng, new_phis = False, laiv_ages = []): 
        
    prevs, R0_series, inc, vacc_post_inf = stoc.simul(params, states, incidences, rng, new_phis = new_phis, laiv_ages = laiv_ages)
                
    return prevs, R0_series, inc, vacc_post_inf

def init_params():
    
    tau = 1                                                                     # tau-leaping step size.
    
    term_admin_laiv = "term_2"                                                  # accepted values are "term_1" or "term_2"
    
    t_start, t_end = 0, 365
    t_iiv_start, t_iiv_end = 59, 226                                            # Start and end date of IIV program
    if term_admin_laiv == "term_1":
      t_laiv_start, t_laiv_end = 59, 100                                          # Start and end date of LAIV program
    elif term_admin_laiv == "term_2":  
      t_laiv_start, t_laiv_end = 110, 170                                          # Start and end date of LAIV program
    else:
      raise ValueError('Error: term_admin_laiv can only take on the values "term_1" or "term_2"!')
    nr_days = t_end - t_start
    
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
    
    tot_pop = 1_000_000
    
    # Age distribution with 8 age groups
    frac_ages = np.array([0.010752583, 0.010748152, 0.033891394, 0.08338837, 0.073682919, 0.614286159, 0.129143011, 0.044107411])

    # N holds the population for each of the age groups.    
    N = allocate_rounded(tot_pop, frac_ages)
    
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
    phi_S1 = np.full(nr_age, 0.7)                                               # dummy value to be replaced in ABC / LHS
    phi_V0 = np.full(nr_age, 0.5)                                               # dummy value to be replaced in ABC / LHS
    phi_V1 = np.full(nr_age, 0.3)                                               # dummy value to be replaced in ABC / LHS
    half_life = np.full(nr_age, 150)                                            # dummy value to be replaced in ABC / LHS
    
    gamma = np.full(nr_age, 1/3)

    target_cover = np.array([0.15, 0.15, 0.15, 0.10, 0.10, 0.20, 0.55, 0.55])         # target vaccination coverage in each age group
    
    iiv_rates = np.array([
        0.00266946, 0.00334557, 0.00737364, 0.03162889, 0.07259051,0.1482384,
        0.24145878, 0.36451686, 0.47168377, 0.56855491, 0.65137032, 0.75339126,
        0.81009063, 0.84931351, 0.88319599, 0.91812563, 0.93545698, 0.95122565,
        0.96457128, 0.97785923, 0.9846138 , 0.99014375, 0.99475859, 1.0
        ])
  
    if term_admin_laiv == "term_1":
      laiv_rates = np.array([0.16666667, 0.33333333, 0.5, 0.66666667, 0.83333333, 1])
    elif term_admin_laiv == "term_2":
      laiv_rates = np.array([0.1111111, 0.2222222, 0.3333333, 0.4444444, 0.5555556, 
      0.6666667,0.7777778, 0.8888889, 1])
    
    daily_vacc_admin = np.zeros((nr_days, nr_age))
    
    # Calculate the number of vaccinations administers on each day    
    daily_vacc_admin = allocate_vacc(daily_vacc_admin, target_cover, N, iiv_rates, t_iiv_start, t_iiv_end)

    p_ext = 0.001                                                                 # Percentage of the total population to be randomly infected from the outside throughout the whole run
    
    # Turn off external seeding of infection for testing purposes
    # p_ext = 0                                                  
    
    # format of scenarios dictionary is key = "name" and value = [scenario type, age groups receiving LAIV, coverage percentages for these groups, scenario name]
    scenarios = {
        "A": ["baseline", [], [], "baseline", "Status quo"],
        "B": ["pessimistic", ["5-11"], [0.4], "pessimistic 5-12", "Pessimistic_LAIV_5-12yo"], 
        "E": ["pessimistic", ["5-11", "12-17"], [0.4, 0.4], "pessimistic 5-18", "Pessimistic_LAIV_5-18yo"], 
        "C": ["mid", ["5-11"], [0.6], "mid 5-12", "Mid_LAIV_5-12yo"], 
        "F": ["mid", ["5-11", "12-17"], [0.6, 0.6], "mid 5-18" ,"Mid_LAIV_5-18yo"], 
        "D": ["optimistic", ["5-11"], [0.8], "optimistic 5-12", "Optimistic_LAIV_5-12yo"], 
        "G": ["optimistic", ["5-11", "12-17"], [0.8, 0.8], "optimistic 5-18", "Optimistic_LAIV_5-18yo"],
        "H": ["mid", ["2-4"], [0.4], "mid 2-5", "Mid_LAIV_2-5yo"], 
        "J": ["mid", ["2-4", "5-11"], [0.4, 0.6], "mid 2-12", "Mid_LAIV_2-12yo"], 
        "K": ["mid", ["2-4", "5-11", "12-17"], [0.4, 0.6, 0.6], "mid 2-18", "Mid_LAIV_2-18yo"],
        "L": ["mid", ["5-11"], [0.4], "mid low 5-12", "Mid_Low_LAIV_5-12yo"], 
        "M": ["mid", ["5-11", "12-17"], [0.4, 0.4], "mid Low 5-18" ,"Mid_Low_LAIV_5-18yo"],
        "N": ["mid", ["2-4", "5-11"], [0.4, 0.4], "mid Low 2-12", "Mid_Low_LAIV_2-12yo"],
        "P": ["mid", ["2-4", "5-11", "12-17"], [0.4, 0.4, 0.4], "mid Low 2-18", "Mid_Low_LAIV_2-18yo"],
        "Q": ["mid", ["5-11"], [0.2], "mid Extra Low 5-12", "Mid_Extra_Low_LAIV_5-12yo"],
        "R": ["mid", ["5-11", "12-17"], [0.2, 0.2], "mid Extra Low 5-18", "Mid_Extra_Low_LAIV_5-18yo"],
        }
    
    colors = {
        "A": "black",
        "B": "#1B9E77",
        "C": "#D95F02",
        "D": "#7570B3",
        "E": "#E7298A",
        "F": "#66A61E",
        "G": "#E6AB02",
        "H": "#A6761D",
        "J": "#666666",
        "K": "#1F78B4",
        "L": "#B15928",
        "M": "#A6CEE3",
        "N": "#B2DF8A",
        "P": "#FB9A99",
        "Q": "#FDBF6F",
        "R": "#CAB2D6"
        }
    
    labels = {
        "A": "baseline",
        "B": "pessimistic 5-12", 
        "E": "pessimistic 5-18", 
        "C": "central 5-12", 
        "F": "central 5-18", 
        "H": "central 2-5", 
        "J": "central 2-12",  
        "K": "central 2-18", 
        "D": "optimistic 5-12", 
        "G": "optimistic 5-18",
        "L": "low 5-12",
        "M": "low 5-18",
        "N": "low 2-12",
        "P": "low 2-18",
        "Q": "extra low 5-12",
        "R": "extra low 5-18",
        }
    
    params = {"tau": tau,
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
              "gamma": gamma,
              "half_life": half_life,
              "daily_vacc_admin": daily_vacc_admin,
              "p_ext": p_ext,
              "shift": shift,
              "laiv_rates": laiv_rates,
              "target_cover": target_cover,
              "scenarios": scenarios,
              "colors": colors,
              "labels": labels
              }
    
    return params

def create_states(params): 
    compartments, nr_age, nr_days = itemgetter("compartments", "nr_age", "nr_days")(params)
    states = {key: np.zeros((nr_age, nr_days)) for key in compartments}
    return states

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

if __name__ == "__main__":
    
    current_dir = os.getcwd().lower()
    rng = np.random.default_rng()
    
    params = init_params()
    states = create_states(params)
    compartments, nr_age, nr_days = itemgetter("compartments", "nr_age", "nr_days")(params)
    skip = {"S0", "S1"}
    daily_inc = {key: np.zeros((nr_age, nr_days)) for key in compartments if key not in skip}
    main(params, states, daily_inc, rng)
    
