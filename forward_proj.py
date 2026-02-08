# TODO: Check daily admin vacc helper function below to ensure new values are used of LAIV

import numpy as np
import pandas as pd
import os
import main_sim as sim
from operator import itemgetter
import matplotlib.pyplot as plt
import time
from matplotlib.ticker import StrMethodFormatter, FuncFormatter
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import copy

plt.close("all")                                                                # Close all figures that may be open in memory

local_parallel = True

def forw_proj(parallel = True, nr_cores = 4):
    
    print("Let's get this thing on the road...\n")
    
    output_figs = True
    output_dfs = True
    analyse_ARs = True
    
    curr_dir = os.getcwd().lower()
    particle_addr = os.path.join(curr_dir, "ABC outputs")
    PATH_IN = os.path.join(particle_addr, "accepted_particles.csv")
    SEED_IN = os.path.join(particle_addr, "accepted_seeds.csv")
    ACC_MID_IN = os.path.join(particle_addr, "acc_mid_phis.csv")
    ACC_OPT_IN = os.path.join(particle_addr, "acc_opt_phis.csv")
    ACC_INF_INC = os.path.join(particle_addr, "acc_inf_inc.csv")
    output_folder_path = os.path.join(curr_dir, "Forward projection")
    
    particles = pd.read_csv(PATH_IN)
    part_vars = list(particles.columns)
    nr_parts = len(particles)
    seed_list = np.array(pd.read_csv(SEED_IN))
    acc_mid_phis = np.array(pd.read_csv(ACC_MID_IN))
    acc_opt_phis = np.array(pd.read_csv(ACC_OPT_IN))
    baseline_inf_inc = np.array(pd.read_csv(ACC_INF_INC))
    
    os.makedirs(output_folder_path, exist_ok=True)
    
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
        "K": ["mid", ["2-4", "5-11", "12-17"], [0.4, 0.6, 0.6], "mid 2-18", "Mid_LAIV_2-18yo"]
        }
    
    params = sim.init_params()
    
    # Initialise generic state for simulation
    empty_states = sim.create_states(params) 
    
    records = []  # to store all rows for the final DataFrame
    
    results_list = []
    
    if parallel:
        full_vars_runs = [(part, seed_list[idx][0], acc_mid_phis[idx], acc_opt_phis[idx], scenarios, part_vars, params, empty_states, curr_dir) for idx, part in particles.iterrows()]
        with ProcessPoolExecutor(max_workers=nr_cores) as executor:
            results_list = list(tqdm(executor.map(run_particle_wrapper, full_vars_runs, chunksize=1), total = len(full_vars_runs), desc="Still doing some forecasting magic"))
    
    else:
        for idx, part in particles.iterrows():
            results_list.append(run_forward_proj_particle(part, seed_list[idx][0], acc_mid_phis[idx], acc_opt_phis[idx], scenarios, part_vars, params, empty_states, curr_dir))
            
            if idx % 5 == 0 and idx > 0:
                print(f"Forward simulation of particle {idx} out of {nr_parts} done.")
                print("")
    
    print("Simulation done.\n")
    
    if analyse_ARs:
        
        print("Doing some analyses of high vs low AR's.\n")
        
        AR_list = []
        AR_dict = {}
        
        for i in range(nr_parts):
            AR_dict = {}
            for k in scenarios:
                AR_dict[k] = results_list[i][1][k][:8] + results_list[i][1][k][8:]
            AR_list.append(AR_dict)
                
        total_AR_list = [
            {key: arr.sum() for key, arr in d.items()}
            for d in AR_list
            ]            
        
        sorted_indices = sorted(range(nr_parts), key=lambda i: total_AR_list[i]['A'], reverse=True)
        
        perc = int(np.ceil(0.1 * nr_parts))
        
        highest_AR = sorted_indices[:perc]
        
        lowest_AR = sorted_indices[-perc:]
        
        plot_AR_tot_inf(output_folder_path, scenarios, highest_AR, lowest_AR, results_list)
    
        plot_corr_for_AR(output_folder_path, particles, part_vars, highest_AR, lowest_AR)
    
        print("Finished the AR analyses. Hopefully that helps.\n")
    
    if output_figs:
        print("Time to go to the drawing board.\n")
        
        params = sim.init_params()
        nr_days, t_start, t_end = itemgetter("nr_days", "t_start", "t_end")(params)
        
        thin_space_formatter = FuncFormatter(lambda x, _: f"{int(x):,}".replace(",", "\u2009"))
        
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
            "K": "#1F78B4"
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
            "G": "optimistic 5-18"
            }
        
        days = np.arange(nr_days)
        
        max_baseline = np.max(baseline_inf_inc)
        
        plot_spaghetti_per_scenario(output_folder_path, max_baseline, thin_space_formatter, nr_parts, days, results_list, colors, labels)
        
        plot_AR_swarm(output_folder_path, scenarios, results_list, colors, labels, thin_space_formatter)
        
        print("Drawing is done!\n")
        
    if output_dfs:
        
        print("Now we're saving the particles and all the other goodies for further analysis down the line.\n")
        
        age_groups, N, nr_days, Unvacc_I_comps, Vacc_I_comps = itemgetter("age_groups", "N", "nr_days", "Unvacc_I_comps", "Vacc_I_comps")(params)
        
        # Output dataframes
        column_nms_AR = [age + "_unvacc" for age in age_groups] + [age + "_vacc" for age in age_groups]
        
        AR_per_scenario_age_strat = {s: np.zeros((nr_parts, 2*len(age_groups))) for s in scenarios}
        
        for idx, res in enumerate(results_list):
            for s in res[1]:
                AR_per_scenario_age_strat[s][idx, :] = res[1][s]
            
        particles.to_csv(os.path.join(output_folder_path, "particles.csv"), index=False)
        age_pops_df = pd.DataFrame([np.array(N)], columns = age_groups)
        age_pops_df.to_csv(os.path.join(output_folder_path, "age_pops.csv"), index=False)
        
        # AR_dict[scenario_nm].append(AR_outputs)
        AR_dfs = {s: pd.DataFrame(ARs, columns=column_nms_AR) for s, ARs in AR_per_scenario_age_strat.items()}
        for name, df in AR_dfs.items():
            df.to_csv(os.path.join(output_folder_path, f"{scenarios[name][3]}.csv"), index=False)       
       
        def vacc_string(vacc_stat):
            if vacc_stat == "Unvacc_I_comps":
                return "unvacc"
            elif vacc_stat == "Vacc_I_comps":
                return "vacc"
            else:
                raise ValueError("Error: Vacc_I_comps variable is not in the allowed list!")
       
        comp_dict_vacc_stat = {
            "Vacc_I_comps": Vacc_I_comps,
            "Unvacc_I_comps": Unvacc_I_comps
            } 
       
        for idx, part in particles.iterrows():
            # Run the scenario
            for s in scenarios:
               
                # --- Record infection incidence data per age and day ---
                for age_idx, age_label in enumerate(age_groups):
                    for day in range(nr_days):
                        for vacc_stat in ["Unvacc_I_comps", "Vacc_I_comps"]:
                            records.append({
                                "scenario": scenarios[s][4],
                                "setting": "temperate",
                                "simulation_index": idx + 1,
                                "age_group": age_label,
                                "horizon": day + 1,
                                "target": "infection_incidence_" + vacc_string(vacc_stat),
                                "value": sum(results_list[idx][2][s][c][age_idx, day] for c in comp_dict_vacc_stat[vacc_stat])
                            })

        # --- Convert to DataFrame and save ---
        df_incid = pd.DataFrame(records)
        df_incid.to_parquet(os.path.join(output_folder_path, "df_incid.parquet")) 

def plot_corr_for_AR(fig_folder_path, full_parts, sample_keys, highest_AR, lowest_AR):
    
    fig_filepath = os.path.join(fig_folder_path, "corr_parts_params_scenarios")
    
    beta_0_range = [0.037, 0.055]
    beta_1_range = [0.05, 0.3]
    shift = [40, 190]
    phi_S1_range = [0.2, 0.8]
    phi_V0_range = [0.2, 0.8]
    phi_V1_range = [0.2, 0.8]
    frac_S1_range = [0.2, 0.8]
    
    sample_ranges = {
        "beta_0": beta_0_range,
        "beta_1": beta_1_range,
        "shift": shift,
        "phi_S1": phi_S1_range, 
        "phi_V0": phi_V0_range, 
        "phi_V1": phi_V1_range,
        "frac_S1": frac_S1_range
        }    
    
    with plt.rc_context({
       "axes.titlesize": 40,
       "axes.labelsize": 36,
       "xtick.labelsize": 30,
       "ytick.labelsize": 30,
       "legend.fontsize": 36,
       "legend.title_fontsize": 38,
       "figure.titlesize": 48
       }):
        
        df = pd.DataFrame(full_parts, columns=sample_keys)
    
        # Create a boolean mask for accepted (largest AR) vs rejected (smallest AR) particles
        AR_size = np.full(len(df), "Medium")
        AR_size[highest_AR] = "High"
        AR_size[lowest_AR] = "Low"
        df["AR_category"] = AR_size
        
        plot_labs = [r"$\beta_0$", r"$\beta_1$", r"fraction $S^1$", r"$\phi^{S^1}$", r"$\phi^{V^0}$", r"$\phi^{V^1}$", "Shift"]

        n = len(sample_keys)
        
        col_high, col_mid, col_low = "royalblue", "silver", "firebrick"
        colour_map = {
            "High": col_high,
            "Medium": col_mid,
            "Low": col_low
            }      

        fig, axes = plt.subplots(n, n, figsize=(4*n,4*n))        
        
        fmt_1 = FuncFormatter(lambda x, _: f"{x:.1f}")
        fmt_2 = FuncFormatter(lambda x, _: f"{x:.2f}")
        fmt_int = FuncFormatter(lambda x, _: f"{int(x)}")
        
        for i, y in enumerate(sample_keys):
            for j, x in enumerate(sample_keys):
                ax = axes[i,j]
                if i < j:
                    ax.set_visible(False)
                    continue
                elif i == j:                    
                    bins = np.linspace(min(sample_ranges[x]), max(sample_ranges[x]), 21)

                    for cat in ["Low", "Medium", "High"]:
                        ax.hist(
                            df.loc[df["AR_category"] == cat, x],
                            bins=bins,
                            color=colour_map[cat],
                            alpha=0.5 if cat == "Medium" else 0.7,
                            edgecolor="none",
                            density=True  # optional but often very useful
                        )
                    
                    ax.set_xlim(min(sample_ranges[x]), max(sample_ranges[x]))
                else:
                    ax.scatter(df[x], df[y], color=df["AR_category"].map(colour_map), s=45, alpha=0.5, edgecolors='none')

                # Axis labels
                if i<n-1: 
                    ax.set_xticklabels([])
                else: 
                    ax.set_xlabel(plot_labs[j], labelpad=10)
                    if i == n-1 and j == 0:
                        # bottom-left: 2 decimals
                        ax.xaxis.set_major_formatter(fmt_2)
                
                    elif i == n-1 and j == n-1:
                        # bottom-right: integers only
                        ax.xaxis.set_major_formatter(fmt_int)
                
                    else:
                        # other bottom-row panels
                        ax.xaxis.set_major_formatter(fmt_1)
                    
                if j>0:
                    ax.set_yticklabels([])
                else: 
                    ax.set_ylabel(plot_labs[i], labelpad=10)
                
                # tick label size increase
                ax.tick_params(axis='x')
                ax.tick_params(axis='y')
                
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

        ax_bottom_left = axes[n-1, 0]
        ax_top_left    = axes[0, 0]
        
        # create a secondary y-axis
        ax_top_left_mirror = ax_top_left.twinx()
        
        # copy tick locations from bottom-left x-axis
        ticks = ax_bottom_left.get_xticks()
        
        ax_top_left_mirror.set_yticks(ticks)
        ax_top_left_mirror.set_ylim(ax_bottom_left.get_xlim())
        
        # Copy formatter
        ax_top_left_mirror.yaxis.set_major_formatter(ax_bottom_left.xaxis.get_major_formatter())
        
        # Make it look like a left axis
        ax_top_left_mirror.yaxis.set_label_position("left")
        ax_top_left_mirror.yaxis.tick_left()
        
        # Hide everything except labels
        ax_top_left_mirror.spines["right"].set_visible(False)
        ax_top_left_mirror.spines["top"].set_visible(False)
        ax_top_left_mirror.spines["bottom"].set_visible(False)
        
        # Hide the original histogram y-axis labels
        ax_top_left.tick_params(axis="y", labelleft=False)

        handles = [
            plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col_high, markersize=15, alpha=0.7, label="High"),
            plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col_mid, markersize=10, alpha=0.5, label="Medium"),
            plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col_low, markersize=15, alpha=0.7, label="Low")
        ]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.56,0.93), borderpad =1, frameon=True)

        plt.suptitle("Correlations of parameters — High AR vs Low AR", x=0.5, y=0.98)
        fig.align_xlabels()
        fig.align_ylabels()
        fig.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.95, hspace=0.35, wspace=0.35)        
        fig.savefig(fig_filepath)
        plt.close(fig)

def plot_AR_swarm(output_folder_path, scenarios, results_list, colors, labels, thin_space_formatter):

    # the two variables below are only used to add horizontal jitter to the plot
    width = 0.3
    rng2 = np.random.default_rng(seed=6978)
    
    # Block for AR for each scenario
    fig_filepath_AR_swarm_no_lines = os.path.join(output_folder_path, "AR_swarm_no_lines")
    fig_filepath_AR_swarm_lines = os.path.join(output_folder_path, "AR_swarm_lines")
    
    tick_label_size = 22
    axis_label_size= 24
    
    fig_lines, ax_lines = plt.subplots(figsize=(18,8))
    fig_no_lines, ax_no_lines = plt.subplots(figsize=(18,8))
    
    nr_parts = len(results_list)
    x_store = {idx: [] for idx in range(nr_parts)}
    y_store = {idx: [] for idx in range(nr_parts)}
    
    for x_pos, key in enumerate(scenarios, start=1):
        
        AR = np.zeros(nr_parts)
        
        if x_pos == 1:
            prev_key = key
        
        for idx, part in enumerate(results_list):
            part = part[0]
            AR[idx] += np.sum(part[key])
        
        x_jittered = x_pos + rng2.uniform(-width, width, nr_parts)
        
        ax_no_lines.scatter(x_jittered, AR, s = 20, alpha = 0.6, color=colors[key], edgecolors="none", zorder = 2)
        ax_lines.scatter(x_jittered, AR, s = 20, alpha = 0.6, color=colors[key], edgecolors="none", zorder = 2)
        med = np.median(AR)
        ax_no_lines.hlines(med, x_pos-0.35, x_pos + 0.35, linewidth=2, color="black", zorder = 3)
        ax_lines.hlines(med, x_pos-0.35, x_pos + 0.35, linewidth=2, color="black", zorder = 3)

        if x_pos > 1:
            for idx in range(nr_parts):
                ax_lines.plot([x_store[idx], x_jittered[idx]], [y_store[idx], AR[idx]], color = colors[prev_key], linewidth = 0.6, alpha = 0.3, zorder = 1)

        for idx in range(nr_parts):
            x_store[idx] = x_jittered[idx]
            y_store[idx] = AR[idx]
    
        prev_key = key
    
    for ax in (ax_no_lines, ax_lines):
        ax.tick_params(axis="both", which="major", labelsize=tick_label_size)
    
    ax_no_lines.set_xticks(range(1, len(scenarios)+1))
    ax_no_lines.set_xticklabels([labels[k] for k in scenarios], rotation=45, ha="right")
    ax_lines.set_xticks(range(1, len(scenarios)+1))
    ax_lines.set_xticklabels([labels[k] for k in scenarios], rotation=45, ha="right")
    
    ax_no_lines.set_ylabel("Attack rate", fontsize = axis_label_size)
    # ax_no_lines.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax_no_lines.yaxis.set_major_formatter(thin_space_formatter)
    ax_lines.set_ylabel("Attack rate", fontsize = axis_label_size)
    # ax_lines.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax_lines.yaxis.set_major_formatter(thin_space_formatter)
    
    # ax_no_lines.set_title("Attack rates per particle across scenarios")
    ax_no_lines.grid(False)
    # ax_lines.set_title("Attack rates per particle across scenarios")
    ax_lines.grid(False)
    
    fig_no_lines.tight_layout()
    fig_no_lines.savefig(fig_filepath_AR_swarm_no_lines)
    fig_lines.tight_layout()
    fig_lines.savefig(fig_filepath_AR_swarm_lines)
    plt.close(fig_no_lines)
    plt.close(fig_lines)

def plot_spaghetti_per_scenario(output_folder_path, max_baseline, thin_space_formatter, nr_parts, days, results_list, colors, labels):

    fig_filepath_compare_all_spag = os.path.join(output_folder_path, "compare_all__scenarios_spag")
    
    nr_rows, nr_cols = 4, 3
    fig_compare_all, axes = plt.subplots(nr_rows, nr_cols, figsize=(24, 14), sharex=True)
    axes = axes.flatten() 
    
    ax_map ={
        "A": axes[1], # 0 and 2 are skipped because the top row only contains the baseline figure
        "B": axes[3],       
        "C": axes[4],        
        "D": axes[5],        
        "E": axes[6],       
        "F": axes[7],        
        "G": axes[8],     
        "H": axes[9],
        "J": axes[10],
        "K": axes[11]
        }
    
    tick_label_size = 22
    legend_size = 22
    axis_label_size= 24
    
    axes[0].axis("off")  # left of baseline
    axes[2].axis("off")  # right of baseline
    
    for key, ax in ax_map.items():
        ax.set_ylim(0, max_baseline * 1.1)
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.yaxis.set_major_formatter(thin_space_formatter)
    
    for i in range(nr_parts):
        for key, ax in ax_map.items():
            ax.plot(days, results_list[i][0][key], color = colors[key], linewidth=0.8, alpha=0.6, label= labels[key])
    
    for i, ax in enumerate(axes):
        handles_fig, labels_fig = ax.get_legend_handles_labels()
        if handles_fig:
            handles_fig, labels_fig = handles_fig[0], labels_fig[0]
            ax.legend([handles_fig], [labels_fig], fontsize = legend_size)
            ax.tick_params(axis = "both", which = "major", labelsize = tick_label_size)
            
        row = i // 3
        
        # Determine if this ax is "active" (not turned off)
        if not ax.has_data():  # axes[0] and axes[2] in top row may be off
            continue
        
        # X-label: only on bottom row
        if row < nr_rows-1:
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Day", fontsize = axis_label_size)
        
        # Y-label: leftmost active subplot in row
        leftmost_ax = None
        for j in range(row*3, (row+1)*3):
            if axes[j].has_data():
                leftmost_ax = axes[j]
                break
        if ax != leftmost_ax:
            ax.set_ylabel("")
        else:
            ax.set_ylabel("Daily infection incidence", fontsize = axis_label_size)
    
    # fig_compare_all.suptitle("Forecast: daily infection incidence")
    fig_compare_all.tight_layout()
    fig_compare_all.savefig(fig_filepath_compare_all_spag)
    plt.close(fig_compare_all)

def plot_AR_tot_inf(output_folder_path, scenarios, highest_AR, lowest_AR, results_list):
    fig_filepath_AR_analysis_tot_inf = os.path.join(output_folder_path, "AR_analysis_total_infections")
    
    fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
    axes = axes.flatten()  # flatten to 1D array for easy indexing
    
    # Loop over keys
    for idx, key in enumerate(scenarios):
        ax = axes[idx]
    
        # Plot largest_AR in red
        for i in highest_AR:
            arr = results_list[i][0][key]
            ax.plot(arr, color='royalblue', alpha=0.4)
    
        # Plot smallest_AR in blue
        for i in lowest_AR:
            arr = results_list[i][0][key]
            ax.plot(arr, color='firebrick', alpha=0.4)
    
        ax.set_title(scenarios[key][3])
        ax.grid(False)
    
    axes[4].plot([], [], color='royalblue', label='Highest AR')
    axes[4].plot([], [], color='firebrick', label='Lowest AR')
    axes[4].legend()
    
    # Optional: common labels
    fig.text(0.5, 0.02, 'Time (days)', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'Value', va='center', rotation='vertical', fontsize=14)
    fig.suptitle('Infection incidence for highest vs lowest attack rates ', fontsize=16)
    
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])
    
    # fig.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    fig.savefig(fig_filepath_AR_analysis_tot_inf)
    plt.close(fig)        

def update_generic(params, particle, keys):
    
    nr_days, nr_age, t_start, t_end, shift = itemgetter("nr_days", "nr_age", "t_start", "t_end", "shift")(params)
    
    # Update values of parameters that are the same for all scenarios with this particle
    beta_0_val = particle[keys.index("beta_0")]
    beta_1_val = particle[keys.index("beta_1")]
    shift = particle[keys.index("shift")]
    
    t = np.linspace(t_start, t_end, num = nr_days)
    beta = beta_0_val * (1 +  beta_1_val * np.sin(( 2 * np.pi * (t - shift) ) / nr_days))
    
    # In the following expressions, it is assumed that the fraction in S1, as well
    # as phi S1, are constant across all age groups
    frac_S1 = np.full(nr_age, particle[keys.index("frac_S1")])
    phi_S1 = np.full(nr_age, particle[keys.index("phi_S1")])
    phi_V0 = np.full(nr_age, particle[keys.index("phi_V0")])
    phi_V1 = np.full(nr_age, particle[keys.index("phi_V1")])
    
    updated_vals = {
        "beta": beta,
        "phi_S1": phi_S1,
        "frac_S1": frac_S1,
        "phi_V0": phi_V0,
        "phi_V1": phi_V1,
        "shift": shift,
        }
    
    # The generic values for the relevant parameters need to be overwritten by their values in the particle, 
    # which is what is done in the following line
    params.update(updated_vals)
    
    return params

def update_daily_vacc(daily_vacc_admin, target_cover, t_laiv_start, t_laiv_end, laiv_rates, laiv_idx, N):
    
    new_daily_vacc_admin = sim.allocate_vacc(daily_vacc_admin, target_cover, N, laiv_rates, t_laiv_start, t_laiv_end)
    
    return daily_vacc_admin
    
def update_scenario_params(scenario, params, mid_phis, opt_phis):
    
    age_groups, laiv_rates, target_cover, t_laiv_start, t_laiv_end, daily_vacc_admin, phi_V0, phi_V1, N  = itemgetter("age_groups", "laiv_rates", "target_cover", "t_laiv_start", "t_laiv_end", 
                                                                                                                   "daily_vacc_admin", "phi_V0", "phi_V1", "N")(params)
    
    laiv_idx = [age_groups.index(age) for age in scenario[1]]
    
    # Change the phi_V0 and pi_V1 values for those groups that receive the LAIV
    if scenario[0] == "mid":
        phi_V0[laiv_idx] = mid_phis[0]
        phi_V1[laiv_idx] = mid_phis[1]
    elif scenario[0] == "optimistic":
        phi_V0[laiv_idx] = opt_phis[0]
        phi_V1[laiv_idx] = opt_phis[1]
    
    # Adjust target coverage of vaccination for current scenario
    target_cover[laiv_idx] = scenario[2]    
    
    daily_vacc_admin[:,laiv_idx] = 0
    # For the ages groups receiving the LAIV, adjust the daily_vacc_admin matrix to reflect LAIV timing and administration rate
    daily_vacc_admin[:, laiv_idx] = update_daily_vacc(daily_vacc_admin[:, laiv_idx], target_cover[laiv_idx], t_laiv_start, t_laiv_end, laiv_rates, laiv_idx, N[laiv_idx])
    
    sample = {
        "phi_V0": phi_V0,
        "phi_V1": phi_V1,
        "target_cover": target_cover,
        "daily_vacc_admin": daily_vacc_admin
        }
    
    # The generic values for the relevant parameters need to be overwritten by their values in the particle, 
    # which is what is done in the following line
    params.update(sample)
    
    return params

def run_forward_proj_particle(part, seed, mid_phis, opt_phis, scenarios, part_vars, params, empty_states, curr_dir):
    
    part = part.to_numpy()
    
    I_comps, nr_days = itemgetter("I_comps", "nr_days")(params)
    
    # Update the particle parameters that do not depend on the scenario, but only on the particle
    params = update_generic(params, part, part_vars)
    params_orig = copy.deepcopy(params)
    
    scenario_inf = {}
    scenario_AR_per_age = {}
    full_output = {}
    
    # Run the scenario
    for idx, s in enumerate(scenarios):
        
        # Initialise the rng for each iteration
        rng = np.random.default_rng(seed)
        
        params = copy.deepcopy(params_orig)
        
        # Update with scenario-specific values
        params = update_scenario_params(scenarios[s], params, mid_phis, opt_phis)        
        
        states = {k: v.copy() for k, v in empty_states.items()}
        incidences = {k: v.copy() for k, v in empty_states.items() if k not in {"S0", "S1"}}
        
        # Run simulation
        prev, R0_series, inc, vacc_post_inf = sim.main(params, states, incidences, curr_dir, rng)
        
        # TODO: Adjust code to export vacc post inf in an age stratified way. I already did this in the abc_code.py and in the laiv_phis.py files
        
        # Total infectious incidence = sum over all infectious states and ages
        inf_inc = np.sum([np.sum(inc[c], axis=0) for c in I_comps], axis=0)
        # Total AR per age group separated between vaccinated and unvaccinated
        total_inf_unvacc = sum(np.sum(inc[key], axis=1) for key in ["I_S0", "I_S1"])
        total_inf_vacc = sum(np.sum(inc[key], axis=1) for key in ["I_V0", "I_V1"])
        
        scenario_inf[s] = inf_inc
        scenario_AR_per_age[s] = np.concatenate([total_inf_unvacc, total_inf_vacc])
        full_output[s] = inc
        
    return scenario_inf, scenario_AR_per_age, full_output

def run_particle_wrapper(args):
    return run_forward_proj_particle(*args)
        
if __name__ == "__main__":
    
    start_time = time.time()
    
    forw_proj(parallel  = local_parallel)
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.3f} seconds\n")
    print("Simulation is done! Forsooth, rejoice!!!\n")