# TODO: Replace linspace by arange in the discretisation of beta
# TODO: Add mean of forecast for each scenario for comparison
# TODO: Add tracking of vaccination while infectious
# TODO: Add delay to onset of protection after vaccination

# =============================================================================
# The full execution pipeline for the code is:
#    1. Run this file, and adapt the following variables as needed:
#       - ABC_output_results:       activate / suppress output of dataframes and plots of abc runs
#       - seed:                     only used to generate the seed list which goes into the code of each run
#       - seed2:
#       - parallel:                 True = parallel process, False = sequential process
#       - nr_cores:                 number of cores used for parallel processing. This variable is ignored if "parallel" is False
#       - nr_particles:             number of particles to run
#       - peak_time_range:          Acceptance range for peak time of epidemic for each run
#       - AR_range:                 Acceptance range for attack rate for each run
#       - VE_range:                 Acceptance range for vaccine efficacy for each run
#
#   In the list_to_sample_params of function, the following sampling ranges for the ABC samples are defined:

#       - beta_0_range
#       - beta_1_range
#       - shift
#       - half-life
#       - phi_S1_range
#       - phi_V0_range
#       - phi_V1_range
#       - frac_S1_range

#   Multiple model variables are defined in the "init_params" function in the main_sim.py file, and need to be changed there if the model is changed.
#   These are:
#       - tau:                  determines size of tau-leap
#       - t_start / t_end:      determine start and end day of the simulation
#       - age_groups:           determines size of tau-leap
#       - compartments:         determines size of tau-leap
#       - tot_pop:              total population
#       - frac_ages:            fraction of the population in each age group
#       - C:                    Contact matrix
#       - gamma:                recovery rate per age group
#       - target_cover:         target coverage percentage for each age group of vaccination
#       - iiv_rates:            weekly cumulative outroll percentage of IIV
#       - laiv_rates:           weekly cumulative outroll percentage of LAIV
#       - p_ext:                external force of infection
# =============================================================================

import os
import numpy as np
import pandas as pd
import main_sim as sim
import lhs_sampler as lhs
import create_plots as plt_crt
# import laiv_phis as laiv
# import forward_proj as forw
import copy
from operator import itemgetter
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time
import platform

ABC_output_results = True                                                       # Suppress outputting figures and dataframes for ABC (if False) for quicker runs

nr_runs_per_part = 2
accept_perc = 0.9
nr_allowed_rejected = nr_runs_per_part - int(np.ceil(nr_runs_per_part * accept_perc))

paral = True
cores = 12

waning = True

def classic_abc(parallel = True, nr_cores = 12):
    
    # Initialise the random number generator
    seed = 502                                                                  # seed is used to run the simulation
    seed2 = 463                                                                 # seed2 is only used to sample parameters in LHS
    rng = np.random.default_rng(seed)
    rng2 = np.random.default_rng(seed2)
    
    accepted_counter = 0
    
    # Number of samples to produce in the LHS
    nr_particles = 700
    seed_list = rng.integers(0, 1_000_000_000_000_000_000,                          # one seed value per unique combination of particle and run nr of that particle
                             size=nr_particles * nr_runs_per_part) 
    
    if platform.system() == "Windows":
        current_dir = os.getcwd().lower()
    else:
        current_dir = "/pvol"
    
    if waning:    
        fig_folder_path = os.path.join(current_dir, "ABC outputs with waning")
    else:
        fig_folder_path = os.path.join(current_dir, "ABC outputs no waning")
            
    os.makedirs(fig_folder_path, exist_ok=True)    
    
    # ranges for acceptance
    peak_time_range = [147, 289]
    AR_range = [0.1, 0.45]
    VE_range = [0.38, 0.5]
   
    # First function initialises the generic simulation parameters, pre LHS
    # The second one creates a dictionary with the parameters to sample in LHS
    # each with their corresponding upper and lower bounds
    # The third function returns an array where each row corresponds to 
    # to one set of sampled parameters (i.e. a particle) from the LHS
    params = sim.init_params()    
    nr_days, nr_age, pop, I_comps = itemgetter("nr_days", "nr_age", "tot_pop", "I_comps")(params)
    
    params_to_sample = list_to_sample_params(waning)
    sample_array, sample_keys = lhs.lhs(params_to_sample, nr_particles, rng2)
    
    # Initialise state for simulation
    empty_states = sim.create_states(params) 
    
    runs = [(idx, particle, seed_list[idx : nr_particles * nr_runs_per_part : nr_particles]) for idx, particle in enumerate(sample_array)]
    
    common_var_all_runs = {
        "nr_particles": nr_particles,
        "empty_states": empty_states,
        "params": params,
        "current_dir": current_dir,
        "sample_keys": sample_keys,
        "peak_time_range": peak_time_range,
        "AR_range": AR_range,
        "VE_range": VE_range,
        "accept_perc": accept_perc,
        "nr_allowed_rejected": nr_allowed_rejected
        }
    
    results_list = []
    
    if parallel:
        full_vars_runs = [run + (common_var_all_runs, True) for run in runs]
        with ProcessPoolExecutor(max_workers=nr_cores) as executor:
            results_list = list(tqdm(executor.map(run_particle_wrapper, full_vars_runs, chunksize=1), total = len(full_vars_runs), desc="Running particles"))
        
        accepted_counter = sum(res["accepted_part"] for res in results_list)
        
    else:
        for run in runs:
            results_list.append(run_particle(*run, common_var_all_runs, False, accepted_counter))            
            if results_list[-1]["accepted_part"]:
                accepted_counter += 1
        
    if accepted_counter == 0:
        raise ValueError("Error: No particles were accepted!")
                
    end_simul_time = time.time()
    simul_time = end_simul_time - start_time
    print(f"Simulation time was: {simul_time:.3f} seconds\n")
    print(f"There were {accepted_counter} accepted particles out of {nr_particles}.\n")
    
    # -----------------------------
    # OUTPUT RESULTS
    # -----------------------------
    
    if ABC_output_results:

        list_accepted_particles = []
        list_prev = []
        list_inc = []
        list_AR = []
        list_VE = []
        list_VE_alt = []
        list_vacc_post_inf = []
        list_R0_acc = []
        list_R0_rej = []
        
        # the [0] in the lines below is to select only the first run for each particle even for multi-run simulations
        for res in results_list:
            list_prev.append(res["prevalences"][0])
            list_inc.append(res["incidences"][0])
            list_AR.append(res["AR"][0])
            list_VE.append(res["VE"][0])
            list_VE_alt.append(res["VE_alt"][0])
            list_vacc_post_inf.append(np.sum(res["vacc_post_inf"][0]))
            if res["accepted_part"]:
                list_accepted_particles.append(res["idx"])
                list_R0_acc.append(res["R0"][0])
            else:
                list_R0_rej.append(res["R0"][0])

        print("Saving accepted particles.\n")
        acc_part_df = pd.DataFrame(sample_array[list_accepted_particles], columns=sample_keys, index = list_accepted_particles)
        acc_part_df.to_csv(os.path.join(fig_folder_path, "accepted_particles.csv"), index=True)
        acc_seed_list_df = pd.DataFrame(seed_list[list_accepted_particles], columns = ["Accepted_seeds"], index = list_accepted_particles)
        acc_seed_list_df.to_csv(os.path.join(fig_folder_path, "accepted_seeds.csv"), index=True)
        
        print("Saving all particles.\n")
        all_part_df = pd.DataFrame(sample_array, columns=sample_keys)
        all_part_df.to_csv(os.path.join(fig_folder_path, "all_particles.csv"), index=False)
        
        print("Saving infection incidence timeseries for accepted particles.\n")
        inf_inc_accepted = np.zeros((nr_days, len(list_accepted_particles)))
        for idx, part in enumerate(list_accepted_particles):
            # Total infectious incidence = sum over all infectious states and ages
            inf_inc_accepted[:, idx] = np.sum([np.sum(list_inc[part][c], axis=0) for c in I_comps], axis=0)
        acc_inf_inc_df = pd.DataFrame(inf_inc_accepted)
        acc_inf_inc_df.to_csv(os.path.join(fig_folder_path, "acc_inf_inc.csv"), index=False)

        fig_spag_plot = False
        fig_corr_accepted_particles = False
        fig_inf_per_comp_for_age_groups = False
        VE_plots = False
        VE_comparison_plot = False
        boxplot_vacc_post_inf = False
        R0_plot = False
        
        plt_crt.outputs(list_prev, list_inc, list_AR, list_VE, list_VE_alt, params, list_accepted_particles, list_vacc_post_inf, peak_time_range, 
                        AR_range, VE_range, sample_array, sample_keys, fig_spag_plot, fig_inf_per_comp_for_age_groups, VE_comparison_plot,
                        fig_corr_accepted_particles, VE_plots, boxplot_vacc_post_inf, params_to_sample, fig_folder_path, parallel, nr_cores, R0_plot, list_R0_acc, list_R0_rej)
        
def list_to_sample_params(waning):    
   
# =============================================================================
# List the parameters to sample via LHS, with their upper and lower bounds
# These replace the default ones given in the regular initialisation of the
# simulation parameters
# =============================================================================

    # beta_0_range = [0.023, 0.047]
    beta_0_range = [0.037, 0.055]
    beta_1_range = [0.05, 0.3]
    shift = [40, 190]
    phi_S1_range = [0.2, 0.8]
    phi_V0_range = [0.2, 0.8]
    phi_V1_range = [0.2, 0.8]
    frac_S1_range = [0.2, 0.8]
    if waning:
      half_life = [250, 1825]
    else:
      half_life = [365_000, 365_250]
    
    params_to_sample = {
        "beta_0": beta_0_range,
        "beta_1": beta_1_range,
        "shift": shift,
        "phi_S1": phi_S1_range, 
        "phi_V0": phi_V0_range, 
        "phi_V1": phi_V1_range,
        "frac_S1": frac_S1_range,
        "half_life": half_life 
        }    
    
    return params_to_sample
    
# -----------------------------
# Update params with particle
# -----------------------------
def update_params(params, particle, keys, it = None):
   # TODO: Check whether the it variable can ever be a non-default value, and if not, remove this
    nr_age, nr_days, t_start, t_end = itemgetter("nr_age", "nr_days", "t_start", "t_end")(params)
        
    beta_0_val = particle[keys.index("beta_0")]
    beta_1_val = particle[keys.index("beta_1")]
    shift = particle[keys.index("shift")]
    
    t = np.linspace(t_start, t_end, num = nr_days)
    beta = beta_0_val * (1 +  beta_1_val * np.sin(( 2 * np.pi * (t - shift) ) / nr_days))
    
    # Sanity check on beta
    if np.any(beta < 0 ):
        raise ValueError(f"Error in particle {it}: Beta is not allowed to be negative!")
    
    # In the following expressions, it is assumed that the fraction in S1, as well
    # as all the phi's, are constant across all age groups
    frac_S1 = np.full(nr_age, particle[keys.index("frac_S1")])
    phi_S1 = np.full(nr_age, particle[keys.index("phi_S1")])
    phi_V0 = np.full(nr_age, particle[keys.index("phi_V0")])
    phi_V1 = np.full(nr_age, particle[keys.index("phi_V1")])
    
    half_life = np.full(nr_age, particle[keys.index("half_life")])
    
    sample = {
        "beta": beta,
        "phi_S1": phi_S1,
        "phi_V0": phi_V0,
        "phi_V1": phi_V1,
        "frac_S1": frac_S1,
        "shift": shift,
        "half_life": half_life
        }
    
    # The generic values for the relevant parameters need to be overwritten by their values in the particle, 
    # which is what is done in the following line
    params.update(sample)
    
    return params

# -----------------------------
# 3) SIMULATOR
# -----------------------------
def simulate_model(params, states, incidences, rng):
    
    prevs, R0_series, incidences, vacc_post_inf = sim.main(params, states, incidences, rng)
    
    return prevs, R0_series, incidences, vacc_post_inf

# -----------------------------
# 4) SUMMARY STATISTICS
# -----------------------------
def compute_summary(prev, inc, params, peak_time_range, AR_range, VE_range, vacc_post_inf):
# =============================================================================
#     For each run, compute summary statistics to be used in ABC 
#     acceptance / rejection algorithm. If all summary statistics fall within
#     accepted range, add the particle to the accepted list.
# =============================================================================
    nr_days, t_end, unvacc_comps, vacc_comps, I_comps, frac_S1, phi_S1, phi_V0, phi_V1  = itemgetter("nr_days", "t_end", "Unvacc_comps", "Vacc_comps", 
                                                                                                     "I_comps", "frac_S1", "phi_S1", "phi_V0", "phi_V1")(params)
    
    daily_inf_from_unvacc = sum(inc[c].sum(axis=0) for c in ["I_S0", "I_S1"])
    daily_inf_from_vacc = sum(inc[c].sum(axis=0) for c in ["I_V0", "I_V1"])
    daily_inf = daily_inf_from_vacc + daily_inf_from_unvacc
    
    total_vacc_post_inf = np.sum(vacc_post_inf)
    
    # First we calculate the VE
    total_vaccinated_pop = 0
    for c in vacc_comps:
        total_vaccinated_pop += np.sum(prev[0][c][:, -1])
    total_vaccinated_pop -= total_vacc_post_inf
        
    total_unvaccinated_pop = 0
    for c in unvacc_comps:
        total_unvaccinated_pop += np.sum(prev[0][c][:, -1])
    total_unvaccinated_pop += total_vacc_post_inf
    
    AR_vacc = np.sum(daily_inf_from_vacc) / total_vaccinated_pop
    AR_unvacc = np.sum(daily_inf_from_unvacc) / total_unvaccinated_pop
    VE = 1 - (AR_vacc / AR_unvacc)
    
    VE_alt = 1 - (( (1-frac_S1[0]) * phi_V0[0] + frac_S1[0] * phi_V1[0] )/((1-frac_S1[0]) + frac_S1[0] * phi_S1[0]))
   
    peak_time = np.argmax(daily_inf)     
   
    total_infected = np.sum(daily_inf)
    AR = total_infected / params["tot_pop"] 
   
    # Basic sanity check
    if AR < 0 or AR > 1:
        raise ValueError("Error: AR is outside the valid range [0, 1]!") 
   
    VE_ok = min(VE_range) <= VE <= max(VE_range)
    peak_time_ok = min(peak_time_range) <= peak_time <= max(peak_time_range)
    AR_ok = min(AR_range) <= AR <= max(AR_range)
   
    if not (VE_ok and peak_time_ok and AR_ok):
        return False, AR, VE, VE_alt
    else:
        return True, AR, VE, VE_alt

def run_particle(idx, particle, seed_list, common_var_all_runs, parallel, accepted_counter = 0):
    
    nr_particles, empty_states, params, current_dir, sample_keys, peak_time_range, AR_range, VE_range,nr_allowed_rejected = itemgetter("nr_particles","empty_states", 
                                "params", "current_dir", "sample_keys", "peak_time_range", "AR_range", "VE_range", "nr_allowed_rejected")(common_var_all_runs)    
    
    acc_list = []
    AR_list = []
    VE_list = []
    VE_alt_list = []
    prev_list = []
    inc_list = []
    R0_list = []
    vacc_post_inf_list = []
    
    params_run = copy.deepcopy(params)

    # Update parameters for this particle
    params_run = update_params(params_run, particle, sample_keys, idx)
    
    accepted_part = True
    
    for it, seed in enumerate(seed_list):
        
        rng = np.random.default_rng(seed)
        
        # Copy states and params for this particle
        states = {k: v.copy() for k, v in empty_states.items()}
        incidences = {k: v.copy() for k, v in empty_states.items() if k not in {"S1"}}
    
        # Run simulation
        prevalences, R0_series, incidences, vacc_post_inf = simulate_model(params_run, states, incidences, rng)
    
        # Compute summary statistics & acceptance
        accepted, AR, VE, VE_alt = compute_summary(prevalences, incidences, params_run, peak_time_range, AR_range, VE_range, vacc_post_inf)
        
        acc_list.append(accepted)
        AR_list.append(AR)
        VE_list.append(VE)
        VE_alt_list.append(VE_alt)
        prev_list.append(prevalences)
        inc_list.append(incidences)
        R0_list.append(R0_series[0][0])
        vacc_post_inf_list.append(vacc_post_inf)
        
        nr_rejected = acc_list.count(False)
        if nr_rejected > nr_allowed_rejected:
            accepted_part = False
            break
    
    if (idx+1) % 50 == 0 and not parallel:
        print(f"Percentage of particles done: {(idx+1)*100 / nr_particles:.0f}%")
        print(f"Percentage of done particles accepted: {accepted_counter*100 / (idx+1):.0f}%\n")
            
    return {
        "idx": idx,
        "seed": seed_list,
        "particle": particle,
        "prevalences": prev_list,
        "incidences": inc_list,
        "accepted_list": acc_list,
        "vacc_post_inf": vacc_post_inf_list,
        "AR": AR_list,
        "VE": VE_list,
        "VE_alt": VE_alt_list,
        "R0": R0_list,
        "accepted_part": accepted_part
    }

def run_particle_wrapper(args):
    return run_particle(*args)

if __name__ == "__main__":
    
    start_time = time.time()
    
    classic_abc(parallel = paral, nr_cores = cores)
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Total initial ABC elapsed time: {elapsed:.3f} seconds\n")
    print("Simulation is done! Huzzah!!!\n")
    
    # start_time_laiv = time.time()

    # laiv.laiv_phis(parallel = paral, waning = waning, nr_cores = cores)

    # end_time_laiv = time.time()
    # elapsed = end_time_laiv - start_time_laiv
    # print(f"Simulation for new phi values elapsed time: {elapsed:.3f} seconds\n")
    print("Calculation of phi's is done! Booyah!!!\n")
    
    # start_time_forw = time.time()
    
    # forw.forw_proj(parallel = paral, nr_cores = cores)
    
    # end_time = time.time()
    # elapsed = end_time - start_time_forw
    # print(f"Forward project elapsed time: {elapsed:.3f} seconds\n")
    # print("Simulation is done! Forsooth, rejoice!!!\n")
    
    # print(f"Total time for full run: {end_time-start_time:.3f} seconds\n")
    
