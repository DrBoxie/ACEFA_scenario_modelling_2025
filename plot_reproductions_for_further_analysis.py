
import numpy as np
import os
import pandas as pd
import pyarrow.dataset as ds
import matplotlib.pyplot as plt
import main_sim as sim
from operator import itemgetter
from matplotlib.ticker import FuncFormatter
import time
from matplotlib.lines import Line2D
# from matplotlib.ticker import StrMethodFormatter
# import gzip
# import pickle

def plot_reproduction_calcs():
    
    curr_dir = os.getcwd().lower()
    files_folder_path = os.path.join(curr_dir, "Forward projection")
    filename_full_df = os.path.join(files_folder_path, "df_incid_summed.parquet")
    dataset = ds.dataset(filename_full_df, format = "parquet")
    
    selected_parts = [55, 57, 67, 40]
    filtered_ds = ds.field("simulation_index").isin(selected_parts)
    
    dfs = []
    
    scanner = dataset.scanner(filter = filtered_ds, batch_size = 50_000)
    
    i = 1
    for batch in scanner.to_batches():
        if batch.num_rows > 0:
            dfs.append(batch.to_pandas())
        if i % 10 == 0:
            print(f"Processed {i} batches.")
        i += 1
        
    df = pd.concat(dfs, ignore_index = True)
    print("")
    
    params = sim.init_params()
    nr_days, t_start, t_end = itemgetter("nr_days", "t_start", "t_end")(params)
    
    col_map = {
        55: "most_neg",
        57: "no_change",
        67: "med_red",
        40: "max_red"
        }
    
    colors = {
        "most_neg": "#D95F02",
        "no_change": "#666666",
        "med_red": "#1F78B4",
        "max_red": "#66A61E"
        }
    
    particle_labels = {
        "most_neg": "Most negative reduction",
        "no_change": "No change",
        "med_red": "Median reduction",
        "max_red": "Maximal reduction"}
    
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
    
    age_list = ['<1', '1', '2-4', '5-11', '12-17', '18-64', '65-79', '80+']
    
    scenario_nms = {k: v[4] for k, v in scenarios.items()}
    scenario_labels = {k: v[3].replace("mid", "central") for k, v in scenarios.items()}
    
    legend_handles = [
        Line2D([0], [0], color=colors["most_neg"], lw=2, label=particle_labels["most_neg"]),
        Line2D([0], [0], color=colors["no_change"], lw=2, label=particle_labels["no_change"]),
        Line2D([0], [0], color=colors["med_red"], lw=2, label=particle_labels["med_red"]),
        Line2D([0], [0], color=colors["max_red"], lw=2, label=particle_labels["max_red"])
    ]
    
    thin_space_formatter = FuncFormatter(lambda x, _: f"{int(x):,}".replace(",", "\u2009"))
    
    days = np.arange(nr_days)
    
    # plot_spaghetti_per_scenario(files_folder_path, thin_space_formatter, days, df, colors, col_map, particle_labels, legend_handles, scenario_nms, scenario_labels)
    
    # plot_spaghetti_age_strat(files_folder_path, thin_space_formatter, days, df, scenario_nms, scenario_labels, col_map, colors, particle_labels, age_list)
    
    # plot_AR_swarm(files_folder_path, scenarios, df, colors, col_map, particle_labels, thin_space_formatter, legend_handles)
    
    # plot_AR_age_strat_swarm(files_folder_path, scenarios, df, colors, col_map, thin_space_formatter, legend_handles, age_list)
    
    plot_particles_correlation_enriched(files_folder_path)

def plot_particles_correlation_enriched(fig_folder_path):
    
    fig_filepath = os.path.join(fig_folder_path, "corr_parts_params")
    
    
    
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
    
        # Create a boolean mask for accepted vs rejected particles
        accepted_mask = np.zeros(len(df), dtype=bool)
        accepted_mask[list_accepted_particles] = True
        
        # Add a column for particle status
        df["status"] = np.where(accepted_mask, "Accepted", "Rejected")
        
        plot_labs = [r"$\beta_0$", r"$\beta_1$", r"fraction $S^1$", r"$\phi^{S^1}$", r"$\phi^{V^0}$", r"$\phi^{V^1}$", "Shift"]
        
        df_acc = df[df["status"]=="Accepted"]
        df_rej = df[df["status"]=="Rejected"]
    
        n = len(sample_keys)
        col_acc, col_rej = "royalblue", "silver"
    
        fig, axes = plt.subplots(n, n, figsize=(4*n,4*n))        
        
        for i, y in enumerate(sample_keys):
            for j, x in enumerate(sample_keys):
                ax = axes[i,j]
                if i < j:
                    ax.set_visible(False)
                    continue
                elif i == j:
                    ax.hist(df_acc[x], bins=20, color=col_acc, alpha=0.6, edgecolor="none")
                    ax.set_xlim(min(sample_ranges[x]), max(sample_ranges[x]))
                    # ax.set_visible(False)
                else:
                    ax.scatter(df_rej[x], df_rej[y], color=col_rej, s=15, alpha=0.6, edgecolors='none')
                    ax.scatter(df_acc[x], df_acc[y], color=col_acc, s=15, alpha=0.8, edgecolors='none')
    
                # Axis labels
                if i<n-1: 
                    ax.set_xticklabels([])
                else: 
                    ax.set_xlabel(plot_labs[j], labelpad=10)
                    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))
                if j>0:
                    ax.set_yticklabels([])
                else: 
                    ax.set_ylabel(plot_labs[i], labelpad=10)
                
                # tick label size increase
                ax.tick_params(axis='x')
                ax.tick_params(axis='y')
                
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
    
        handles = [
            plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col_acc, markersize=10, alpha=0.6, label="Accepted"),
            plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=col_rej, markersize=10, alpha=0.6, label="Rejected")
        ]
        fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.56,0.93), borderpad =1, frameon=True)
    
        plt.suptitle("ABC Particle Sampling — Accepted vs Rejected", x=0.5, y=0.98)
        # fig.text(0.5, 0.97, "ABC Particle Sampling — Accepted vs Rejected", ha='center', va='top', fontsize=48)
        fig.align_xlabels()
        fig.align_ylabels()
        fig.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.95, hspace=0.35, wspace=0.35)        
    
        fig.savefig(fig_filepath)
        plt.close(fig)

    print("Particle correlation plot done.\n")


def plot_AR_age_strat_swarm(output_folder_path, scenarios, df, colors, col_map, thin_space_formatter, legend_handles, age_list):
    
    print("Plotting age stratified AR swarm plots.\n")
    
    # the two variables below are only used to add horizontal jitter to the plot
    width = 0.3
    rng2 = np.random.default_rng(seed=56978)
    
    fig_filepath_AR_swarm_age_strat = os.path.join(output_folder_path, "AR_swarm_age_stratified")
    
    tick_label_size = 18
    axis_label_size = 22
    
    # calculate AR for each unique run per age group
    df_age = df.groupby(
        ["scenario", "simulation_index", "age_group", "run_nr"], as_index=False)["value"].sum()
    
    nr_age = len(age_list)
    # nr_runs = df_age["run_nr"].nunique()
    
    fig_lines, axes = plt.subplots(
        nrows=nr_age,
        ncols=1,
        figsize=(18, 3.2 * nr_age),
        sharex=True
    )
    
    for age_idx, age in enumerate(age_list):
    
        ax = axes[age_idx]
    
        df_age_sub = df_age[df_age["age_group"] == age]
    
        x_store_prev, y_store_prev = {}, {}
        x_store_new, y_store_new = {}, {}
    
        for x_pos, key in enumerate(scenarios, start=1):
    
            AR_scen = df_age_sub[
                df_age_sub["scenario"] == scenarios[key][4]
            ]
    
            for part in AR_scen["simulation_index"].unique():
    
                AR_part = AR_scen[
                    AR_scen["simulation_index"] == part
                ].sort_values("run_nr")
    
                x_jittered = x_pos + rng2.uniform(-width, width, len(AR_part))
                vals = AR_part["value"].to_numpy()
    
                ax.scatter(
                    x_jittered,
                    vals,
                    s=20,
                    alpha=0.6,
                    color=colors[col_map[part]],
                    edgecolors="none",
                    zorder=2
                )
    
                x_store_new[part] = x_jittered
                y_store_new[part] = vals
    
            if x_pos > 1:
                for part in AR_scen["simulation_index"].unique():
                    for run_idx in range(len(x_store_new[part])):
                        ax.plot(
                            [x_store_prev[part][run_idx], x_store_new[part][run_idx]],
                            [y_store_prev[part][run_idx], y_store_new[part][run_idx]],
                            color=colors[col_map[part]],
                            linewidth=0.6,
                            alpha=0.3,
                            zorder=1
                        )
    
            for part in AR_scen["simulation_index"].unique():
                x_store_prev[part] = x_store_new[part]
                y_store_prev[part] = y_store_new[part]
    
        ax.tick_params(axis="both", which="major", labelsize=tick_label_size)
        ax.set_ylabel(str(age), fontsize=axis_label_size)
        ax.yaxis.set_major_formatter(thin_space_formatter)
        ax.grid(False)
    
    axes[-1].set_xticks(range(1, len(scenarios) + 1))
    axes[-1].set_xticklabels(
        [scenarios[k][3] for k in scenarios],
        rotation=45,
        ha="right"
    )
    
    fig_lines.supylabel("Attack rate", fontsize=axis_label_size)
    
    fig_lines.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.99, 0.98),
        fontsize=18,
        frameon=True
    )
    
    fig_lines.tight_layout()
    fig_lines.savefig(fig_filepath_AR_swarm_age_strat)
    plt.close(fig_lines)
    
def plot_AR_swarm(output_folder_path, scenarios, df, colors, col_map, labels, thin_space_formatter, legend_handles):
    
    print("Plotting AR swarm plots.\n")
    
    # the two variables below are only used to add horizontal jitter to the plot
    width = 0.3
    rng2 = np.random.default_rng(seed=56978)
    
    fig_filepath_AR_swarm_lines = os.path.join(output_folder_path, "AR_swarm_lines_special_particles")
    
    tick_label_size = 22
    axis_label_size= 24
    
    fig_lines, ax_lines = plt.subplots(figsize=(18,8))
    
    figures = [(fig_lines, fig_filepath_AR_swarm_lines)]
    
    # nr_parts = df["simulation_index"].nunique()
    nr_runs = df["run_nr"].nunique()
    # nr_points = nr_parts * nr_runs
    
    x_store_prev, y_store_prev = {}, {}
    x_store_new, y_store_new = {}, {}    
    
    # After the following line, df contains the AR per age group for each run of each particle
    df = df.groupby(["scenario", "simulation_index", "age_group", "run_nr"], as_index = False)["value"].sum()
    
    # AR contains the total AR rate across all age groups for each run of each particle
    AR = df.groupby(["scenario", "simulation_index", "run_nr"], as_index = False)["value"].sum()
    
    for x_pos, key in enumerate(scenarios, start=1):
        
        AR_scen = AR[AR["scenario"] == scenarios[key][4]]
        
        for part in AR_scen["simulation_index"].unique():
        
            x_jittered = x_pos + rng2.uniform(-width, width, nr_runs)
            vals = AR_scen.loc[AR_scen["simulation_index"] == part, "value"].to_numpy()
            ax_lines.scatter(x_jittered, vals, s = 20, alpha = 0.6, color=colors[col_map[part]], edgecolors="none", zorder = 2)
            
            x_store_new[part] = x_jittered
            y_store_new[part] = vals
        
        if x_pos > 1:
            for part in AR_scen["simulation_index"].unique():
                for run in range(nr_runs):
                    ax_lines.plot([x_store_prev[part][run], x_store_new[part][run]], [y_store_prev[part][run], y_store_new[part][run]], color = colors[col_map[part]], linewidth = 0.6, alpha = 0.3, zorder = 1)
        
        for part in AR_scen["simulation_index"].unique():
            x_store_prev[part] = x_store_new[part]
            y_store_prev[part] = y_store_new[part]
        
    ax_lines.tick_params(axis="both", which="major", labelsize=tick_label_size)
    ax_lines.set_xticks(range(1, len(scenarios)+1))
    ax_lines.set_xticklabels([])
    # ax.set_xticklabels([labels[k] for k in scenarios], rotation=45, ha="right")
    ax_lines.set_ylabel("Attack rate", fontsize = axis_label_size)
    ax_lines.yaxis.set_major_formatter(thin_space_formatter)
    ax_lines.grid(False)
    
    ax_lines.set_xticklabels([scenarios[k][3] for k in scenarios], rotation=45, ha="right")
    
    fig_lines.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor = (0.99, 0.8),
        fontsize = 18,
        frameon=True
        )
    
    for fig, filepath in figures:
        fig.tight_layout()
        fig.savefig(filepath)
        plt.close(fig)

def plot_spaghetti_per_scenario(output_folder_path, thin_space_formatter, days, df, colors, col_map, labels, legend_handles, scenario_nms, scenario_labels):
    
    print("Plotting spaghetti plots.\n")
    
    filename_spag = "spag_special_particles"
    
    fig_filepath_spag = os.path.join(output_folder_path, filename_spag)
    
    nr_rows, nr_cols = 4, 3
    fig_compare_all, axes = plt.subplots(nr_rows, nr_cols, figsize=(23, 16), sharex=True)
    axes = axes.flatten() 
    
    max_baseline = 13_000
    
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
        ax.yaxis.set_major_formatter(thin_space_formatter)
        ax.text(
            0.03, 0.95, 
            scenario_labels[key],
            transform = ax.transAxes,
            fontsize = 18,
            va = "top",
            ha = "left"
            )
    
    # sums infection incidences across all age groups
    df = df.groupby(["scenario", "simulation_index", "horizon", "run_nr"], as_index = False)["value"].sum()
    
    for key, ax in ax_map.items():
        sc_inf = df[ df["scenario"] == scenario_nms[key]]
        
        for (part, run), inf_inc in sc_inf.groupby(["simulation_index", "run_nr"]):
            inf_inc = inf_inc.sort_values("horizon")
            ax.plot(inf_inc["horizon"], inf_inc["value"], color = colors[col_map[part]], linewidth = 0.7, alpha = 0.4)
    
    for i, ax in enumerate(axes):
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
        elif row in [1, 3]:
            ax.set_ylabel("")
        else:
            ax.set_ylabel("Daily infection incidence", fontsize = axis_label_size)
    
    fig_compare_all.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor = (0.01, 0.96),
        fontsize=legend_size,
        frameon=True
        )
    
    fig_compare_all.suptitle("Representative particle behaviour: all runs", fontsize = 28, y = 0.99)
    fig_compare_all.tight_layout(rect=[0,0,0.97,0.97])
    fig_compare_all.savefig(fig_filepath_spag)
    plt.close(fig_compare_all)


def plot_spaghetti_age_strat(output_folder_path, thin_space_formatter, days, df, scenario_nms, scenario_labels, col_map, colors, labels, age_list):
    
    print("Plotting spaghetti plots per age group.\n")
    
    scenarios_to_plot = ["A", "F", "G"]
    selected_sc_nms = [scenario_nms[s] for s in scenarios_to_plot]
    
    parts_to_plot = ["most_neg", "med_red"]
    selected_parts = [k for k,v in col_map.items() if v in parts_to_plot]
    
    df = df[ (df["scenario"].isin(selected_sc_nms)) & (df["simulation_index"].isin(selected_parts)) ]
    
    fig_filepath_spag = os.path.join(output_folder_path, "spag_age_strat_special_particles")
    
    nr_age = len(age_list)
    nr_scens = len(scenarios_to_plot)
    nr_runs = df["run_nr"].nunique()
    
    fig_age_scen_spag, axes = plt.subplots(nr_age, nr_scens, figsize=(24, 32), sharex=True)
    
    tick_label_size = 22
    legend_size = 22
    axis_label_size= 24
    
    legend_handles = [
        Line2D(
            [0], [0],
            color = colors[p],
            lw = 2,
            label = labels[p]
            )
        for p in parts_to_plot
        ]
    
    for sc_idx, scen in enumerate(scenarios_to_plot):
        for age_idx, age in enumerate(age_list):
            ax = axes[age_idx, sc_idx]
            ax.yaxis.set_major_formatter(thin_space_formatter)
            
            scen_age_dat = df[ (df["scenario"] == scenario_nms[scen]) & (df["age_group"] == age)]
            
            for part in selected_parts:
                for run in range(nr_runs):
                    fig_dat = scen_age_dat[ (scen_age_dat["simulation_index"] == part) & (scen_age_dat["run_nr"] == run)].sort_values("horizon")
                    ax.plot(days, fig_dat["value"], color = colors[col_map[part]], linewidth = 0.6, alpha = 0.2)
    
    for i, ax in enumerate(axes.flat):
        ax.tick_params(axis = "both", which = "major", labelsize = tick_label_size)
                    
        row = i // 3
        
        # x-label only on bottom row
        if row < nr_age-1:
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Day", fontsize = axis_label_size)
        
        # y-label only on leftmost column
        if i % 3 == 0:
            ax.set_ylabel(f"{age_list[i // 3]}", fontsize = axis_label_size, fontweight = "bold")
        else:
            ax.set_ylabel("")
    
    # add titles to the columns of figure
    for sc_idx, scen in enumerate(scenarios_to_plot):
        axes[0, sc_idx].set_title(scenario_labels[scen], fontsize = 26, pad = 18, fontweight = "bold")
    
    fig_age_scen_spag.text(
        0.03, 0.5,
        "Daily infection incidence",
        ha = "center",
        va = "center",
        rotation = 90,
        fontsize = axis_label_size
        )
    
    fig_age_scen_spag.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor = (0.5, 0.01),
        ncol = len(selected_parts),
        # fontsize=legend_size,
        frameon=True,
        prop={"size": legend_size, "weight": "bold"}
        )
    
    fig_age_scen_spag.suptitle("Age stratified infection incidence for particles of interest", fontsize = 30, y = 0.985)
    fig_age_scen_spag.tight_layout(rect=[0.04, 0.03, 0.98,0.97])
    fig_age_scen_spag.savefig(fig_filepath_spag)
    plt.close(fig_age_scen_spag)

if __name__ == "__main__":
    
    start_time = time.time()
    
    plot_reproduction_calcs()
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.3f} seconds\n")
    print("Plotting of special particles done. Yippie-ki-yay!!!\n")