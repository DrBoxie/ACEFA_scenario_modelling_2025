# TODO: Correct figure code to compare different phi's

import numpy as np
import pandas as pd
import os
import main_sim as sim
import forward_proj as forw
from operator import itemgetter
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import time
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import copy

new_phis = True
local_parallel = True

def laiv_phis(new_phis = True, parallel = True, nr_cores = 12):
    
    print("\nStarting the simulation.\n")
    
    # This file uses 2 rng's
    # rng is used to run the simulation itself with the seeds coming from the initial ABC
    # rng2 is used at the very start of the simulation to generate the seeds with which the new phi's will be generated. And it is then used with those seeds to generate the new phi's for each particle
    seed = 147
    rng2 = np.random.default_rng(seed)
    
    make_figs = False

    VE_accept_range = {
        "pessimistic": [0.38, 0.50],                                                # same as baseline, so no need to run separately
        "mid" : [0.44, 0.69],
        "optimistic": [0.63, 0.88]
        }

    LAIV_ages = ['2-4', '5-11', '12-17']

    phi_V0_range = [0, 0.8]
    phi_V1_range = [0, 0.8]

    curr_dir = os.getcwd().lower()
    particle_addr = os.path.join(curr_dir, "ABC outputs")
    PATH_IN = os.path.join(particle_addr, "accepted_particles.csv")
    SEED_IN = os.path.join(particle_addr, "accepted_seeds.csv")
    VE_fig_out_addr = os.path.join(particle_addr, "acc_rej_phis_VE")

    results_list = []
    acc_pes_phi = []
    acc_mid_phi = []
    acc_opt_phi = []

    list_scenarios = []
    list_phis = []
    list_accepted = []

    particles_df = pd.read_csv(PATH_IN, index_col = 0)
    particles, part_vars = particles_df.to_numpy(), list(particles_df.columns)
    
    # part_vars = list(particles.columns)
    nr_parts = len(particles)
    seed_list = np.array(pd.read_csv(SEED_IN, index_col = 0))
    
    # Generate seeds we will use for the generation of the phi's (the first element of the tuple is for the pessimistic scenario, the second for the mid scenario, and the third for the optimistic one)
    phi_seed_list = list(zip(rng2.integers(0, 1_000_000_000_000, size = nr_parts), rng2.integers(0, 1_000_000_000_000, size = nr_parts), rng2.integers(0, 1_000_000_000_000, size = nr_parts)))
    
    params = sim.init_params()
    nr_days, nr_age, pop, t_end, unvacc_comps, vacc_comps, I_comps = itemgetter("nr_days", "nr_age", "tot_pop", "t_end", "Unvacc_comps", "Vacc_comps", "I_comps")(params)    
    
    # Initialise state for simulation
    empty_states = sim.create_states(params) 
    
    # Find index of age groups receiving LAIV
    laiv_idx = [params["age_groups"].index(age) for age in LAIV_ages]
    
    phi_V0_min, phi_V0_max = min(phi_V0_range), max(phi_V0_range)
    phi_V1_min, phi_V1_max = min(phi_V1_range), max(phi_V1_range)
    
    if parallel:
        
        full_vars_runs = [(idx, part, seed_list[idx][0], phi_seed_list[idx]) + (params, part_vars, phi_V1_min, phi_V1_max, phi_V0_max, LAIV_ages, laiv_idx, empty_states, curr_dir, VE_accept_range, True) for idx, part in enumerate(particles)]
        with ProcessPoolExecutor(max_workers=nr_cores) as executor:
            results_list = list(tqdm(executor.map(run_particle_wrapper, full_vars_runs, chunksize=1), total = len(full_vars_runs), desc="Running simulations for phi's"))
    
    else:
        for idx, part in enumerate(particles):
            
            results_list.append(individual_part_phi_calc(idx, part, seed_list[idx][0], phi_seed_list[idx], params, part_vars, phi_V1_min, phi_V1_max, phi_V0_max, LAIV_ages, laiv_idx, empty_states, curr_dir, VE_accept_range, False))
            
            if idx % 5 == 0 and idx > 0:
                print(f"Finished particle number {idx}.\n")
                
    for i in range(len(results_list)):
        acc_pes_phi.append(results_list[i][0])
        acc_mid_phi.append(results_list[i][1])
        acc_opt_phi.append(results_list[i][2])
        list_scenarios.append(results_list[i][3])
        list_phis.append(results_list[i][4])
        list_accepted.append(results_list[i][5])
        
    print("Simulation for phi's is finished.\n")
    # -----------------------------
    # OUTPUT RESULTS
    # -----------------------------
    
    if make_figs:
        print("Plotting figures.")
        VE_plots(list_scenarios, list_phis, list_accepted, VE_fig_out_addr)
    
    df_acc_pes_phis = pd.DataFrame(acc_pes_phi, columns=["phi_V0", "phi_V1"], index = particles_df.index)
    df_acc_mid_phis = pd.DataFrame(acc_mid_phi, columns=["phi_V0", "phi_V1"], index = particles_df.index)
    df_acc_opt_phis = pd.DataFrame(acc_opt_phi, columns=["phi_V0", "phi_V1"], index = particles_df.index)
    
    df_acc_pes_phis.to_csv(os.path.join(particle_addr, "acc_pes_phis.csv"), index=True)
    df_acc_mid_phis.to_csv(os.path.join(particle_addr, "acc_mid_phis.csv"), index=True)
    df_acc_opt_phis.to_csv(os.path.join(particle_addr, "acc_opt_phis.csv"), index=True)

def VE_plots(list_scenarios, list_phis, list_accepted, VE_fig_out_addr):

    # Mapping for x-axis positions
    scenario_to_x = {'pessimistic': 1, 'mid': 3, 'optimistic': 5}

    # Markers for the three values in list_phis
    markers = ['o', 's', '^']  # circle, square, triangle
    marker_labels = ['phi_S1', 'phi_V0', 'phi_V1']

    # Bigger horizontal offsets for separation
    offsets = [-0.5, 0, 0.5]

    fig = plt.figure(figsize=(8,5))

    for scenario, phis, accepted in zip(list_scenarios, list_phis, list_accepted):
        x_base = scenario_to_x[scenario]
        for j, phi in enumerate(phis):
            x = x_base + offsets[j] + np.random.uniform(-0.05, 0.05)  # small jitter
            colour = 'royalblue' if accepted else 'silver'
            plt.scatter(x, phi, marker=markers[j], color=colour, edgecolor=None, s=8, alpha=0.7)

    plt.xticks(list(scenario_to_x.values()), list(scenario_to_x.keys()))
    plt.ylim(-0.1, 1.1)
    plt.ylabel('Phi values')
    plt.title('Phi values by scenario')

    # Legend for markers
    marker_handles = [Line2D([0], [0], marker=m, color='w', label=l, markerfacecolor='white',
                             markeredgecolor='black', markersize=10) 
                      for m, l in zip(markers, marker_labels)]

    # Legend for accepted/rejected
    colour_handles = [
        Line2D([0], [0], marker='o', color='w', label='Accepted', markerfacecolor='royalblue', markeredgecolor='black', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Rejected', markerfacecolor='silver', markeredgecolor='black', markersize=10)
    ]

    # Combine legends
    plt.legend(handles=marker_handles + colour_handles, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    plt.tight_layout()

    plt.tight_layout()
    fig.savefig(VE_fig_out_addr, bbox_inches='tight')
    plt.close(fig)

def individual_part_phi_calc(idx, part, run_seed, phi_seeds, params, part_vars, phi_V1_min, phi_V1_max, phi_V0_max, laiv_ages, laiv_idx, empty_states, curr_dir, VE_accept_range, parallel):

    list_phis = []
    list_scenarios = []
    list_accepted = []
    
    for s in ["pessimistic", "mid", "optimistic"]:
        accepted = False
        
        params = forw.update_generic(params, part, part_vars)
        
        vacc_comps, unvacc_comps, phi_S1 = itemgetter("Vacc_comps", "Unvacc_comps", "phi_S1")(params)
        
        if s == "pessimistic":
            rng2 = np.random.default_rng(phi_seeds[0])
        elif s == "mid":
            rng2 = np.random.default_rng(phi_seeds[1])
        elif s == "optimistic":
            rng2 = np.random.default_rng(phi_seeds[2])
        
        iterator = 0
        
        while not accepted:
            
            u_V1 = rng2.uniform()
            phi1 = phi_S1[0] * u_V1

            u_V0 = rng2.uniform()
            phi0 = phi1 + (phi_V0_max - phi1) * u_V0                
            
            params["phi_V0"][laiv_idx] = phi0
            params["phi_V1"][laiv_idx] = phi1
            
            if not(0 <= phi0 <= 1 and 0 <= phi1 <= 1):
                raise ValueError("Error: Phi's have to be between 0 and 1!")
            
            states = {k: v.copy() for k, v in empty_states.items()}
            incidences = {k: v.copy() for k, v in empty_states.items() if k not in {"S0", "S1"}}
            params_run = copy.deepcopy(params)
            
            # Initialise the seed for each iteration
            rng = np.random.default_rng(run_seed)

            # Run simulation
            prev, R0_series, inc, vacc_post_inf = sim.main(params_run, states, incidences, curr_dir, rng, new_phis = new_phis, laiv_ages = laiv_ages)
            
            inf_from_unvacc = sum(np.sum(inc[c][laiv_idx, :]) for c in ["I_S0", "I_S1"])
            inf_from_vacc = sum(np.sum(inc[c][laiv_idx, :]) for c in ["I_V0", "I_V1"])
            
            total_vacc_post_inf = np.sum(vacc_post_inf)
            
            total_vaccinated_pop = 0
            for c in vacc_comps:
                total_vaccinated_pop += np.sum(prev[0][c][laiv_idx, -1])
            total_vaccinated_pop -= total_vacc_post_inf
                
            total_unvaccinated_pop = 0
            for c in unvacc_comps:
                total_unvaccinated_pop += np.sum(prev[0][c][laiv_idx, -1])
            total_unvaccinated_pop += total_vacc_post_inf
            
            AR_vacc = inf_from_vacc / total_vaccinated_pop
            AR_unvacc = inf_from_unvacc / total_unvaccinated_pop
            
            VE = 1 - (AR_vacc / AR_unvacc)
            
            if (min(VE_accept_range[s]) <= VE <= max(VE_accept_range[s])):
                accepted = True
                if s == "pessimistic":
                    acc_pes_phi = [phi0, phi1]
                elif s == "mid":
                    acc_mid_phi = [phi0, phi1]
                else:
                    acc_opt_phi = [phi0, phi1]
            
            list_scenarios.append(s)
            # It is assumed that the phi's have the same value for all age groups initially, and therefore
            # the 0th component is taken as the original value, before finding the new values for the forward projections
            list_phis.append([params["phi_S1"][0], phi0, phi1])
            list_accepted.append(accepted)
            
            iterator += 1
            
            if not parallel:
                if iterator % 20 == 0:
                    print(f"Particle {idx} in scenario {s} has been run {iterator} times\n")
            
    return acc_pes_phi, acc_mid_phi, acc_opt_phi, list_scenarios, list_phis, list_accepted

def run_particle_wrapper(args):
    return individual_part_phi_calc(*args)
    
if __name__ == "__main__":
    
    start_time = time.time()
    laiv_phis(parallel = local_parallel)
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.3f} seconds\n")
    print("Calculation of phi's is done! Booyah!!!\n")
    