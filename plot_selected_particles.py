
import os
import numpy as np
import pandas as pd
import main_sim as sim
import abc_code as abc
import create_plots as plt_crt
import copy
from operator import itemgetter
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time

local_parallel = True

def run_multi_accepted(parallel = True, nr_cores = 12):
    
    print("Running selected accepted particles for further analysis.\n")
    
    nr_parts = 6
    nr_runs_per_part = 50
    
    seed1 = 123568                                                              # seed used only to select which accepted particles to run
    seed2 = 88464                                                               # seed used only to produce seed list for the runs
    rng_choice = np.random.default_rng(seed1)
    rng_seeds_prod = np.random.default_rng(seed2)
    seed_list = rng_seeds_prod.integers(0, 1_000_000_000_000, size=nr_parts * nr_runs_per_part)    # this is the actual list of seeds used for the runs
    
    output_figs = True
    
    curr_dir = os.getcwd().lower()
    particle_addr = os.path.join(curr_dir, "ABC outputs RUN 17042026")
    PATH_IN = os.path.join(particle_addr, "accepted_particles.csv")
    output_folder_path = os.path.join(particle_addr, "Further analysis")

    os.makedirs(output_folder_path, exist_ok=True)
    
    particles = pd.read_csv(PATH_IN)
    # part_vars = list(particles.columns)
    
    selected_idx = rng_choice.choice(len(particles), size = nr_parts, replace = False)
    selected_parts = particles.iloc[selected_idx].to_numpy()
    param_names = sorted(particles.keys())

    params = sim.init_params()    
    # nr_days, nr_age, pop, I_comps = itemgetter("nr_days", "nr_age", "tot_pop", "I_comps")(params)

    # Initialise state for simulation
    empty_states = sim.create_states(params) 
    
    runs = []
    for idx_part, particle in enumerate(selected_parts):
        for run_of_part in range(nr_runs_per_part):
            idx = idx_part * nr_runs_per_part + run_of_part
            runs.append((idx, particle, seed_list[idx]))
    
    common_var_all_runs = {
        "nr_particles": nr_parts,
        "empty_states": empty_states,
        "params": params,
        "current_dir": curr_dir,
        "sample_keys": param_names
        }
    
    results_list = []
    
    if parallel:
        full_vars_runs = [run + (common_var_all_runs, True) for run in runs]
        with ProcessPoolExecutor(max_workers=nr_cores) as executor:
            results_list = list(tqdm(executor.map(run_particle_wrapper, full_vars_runs, chunksize=1), total = len(full_vars_runs), desc="Running particles"))
        
    else:
        for run in runs:
            results_list.append(run_particle(*run, common_var_all_runs, False))

    if output_figs:

        list_prev = []
        list_inc = []
        
        for res in results_list:
            list_prev.append(res["prevalences"][0])
            list_inc.append(res["incidences"][0])

        selected_spag_plot = True
        
        plt_crt.outputs(list_prev, list_inc, [], [], [], params, selected_idx, [], [], 
                        [], [], selected_parts, param_names, False, False, False,
                        False, False, False, {}, output_folder_path, parallel, nr_cores, False, [], [], selected_spag_plot, selected_idx)


def run_particle(idx, particle, seed, common_var_all_runs, parallel):
    
    empty_states, params, current_dir, sample_keys = itemgetter("empty_states", "params", "current_dir", "sample_keys")(common_var_all_runs)    
    
    prev_list = []
    inc_list = []
    R0_list = []
    vacc_post_inf_list = []
    
    params_run = copy.deepcopy(params)

    # Update parameters for this particle
    params_run = abc.update_params(params_run, particle, sample_keys, idx)
        
    rng = np.random.default_rng(seed)
    
    # Copy states and params for this particle
    states = {k: v.copy() for k, v in empty_states.items()}
    incidences = {k: v.copy() for k, v in empty_states.items() if k not in {"S0", "S1"}}

    # Run simulation
    prevalences, R0_series, incidences, vacc_post_inf = abc.simulate_model(params_run, states, incidences, current_dir, rng)
    
    prev_list.append(prevalences)
    inc_list.append(incidences)
    R0_list.append(R0_series[0][0])
    vacc_post_inf_list.append(vacc_post_inf)
            
    return {
        "idx": idx,
        "seed": seed,
        "particle": particle,
        "prevalences": prev_list,
        "incidences": inc_list,
        "vacc_post_inf": vacc_post_inf_list,
        "R0": R0_list,
    }

def run_particle_wrapper(args):
    return run_particle(*args)

if __name__ == "__main__":
    
    start_time = time.time()
    
    run_multi_accepted(parallel  = local_parallel)
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.3f} seconds\n")
    print("Simulation is done! Forsooth, rejoice!!!\n")
    
    
    
    