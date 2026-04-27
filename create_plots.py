import numpy as np
import os
# import psutil
from operator import itemgetter
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter, FuncFormatter, PercentFormatter
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.signal import savgol_filter                                          # Only used for the smoothing of figures (when that option is activated)
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from tqdm import tqdm
import time

plt.close("all")

# Helper function for spagetthi plots layout
def style_spaghetti_plot(params, ax, cumul = False, fig = None, accept_rng = None, t_start = 0, t_end = 365, thin_space_formatter = None, title = None, ylab = None, vert_or_horiz = None, selected_parts = False):
    ax.set_xlabel("Day")
    ax.set_ylabel(ylab, labelpad=10)
    if vert_or_horiz == "vert":
        ax.axvspan(min(accept_rng), max(accept_rng), facecolor="orange", alpha=0.4, label="Acceptance range", edgecolor="none")
    elif vert_or_horiz == "horiz":
        ax.axhspan(min(accept_rng) * params["tot_pop"], max(accept_rng) * params["tot_pop"], facecolor="orange", alpha=0.4, label="Acceptance range", edgecolor="none")
    ax.set_xlim(t_start, t_end)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    ax.yaxis.set_major_formatter(thin_space_formatter)
    leg = ax.legend(loc="upper left" if cumul else "upper right")
    if selected_parts:
        for line in leg.get_lines():
            line.set_linewidth(2.5)
            line.set_alpha(1)
    ax.grid(False)
    fig.suptitle(title)
    fig.tight_layout()

def VE_plots(VE_list, VE_range, list_accepted_particles, fig_folder_path):

    fig_filepath_VE_bp = os.path.join(fig_folder_path, "VE_boxplot")
    
    VEs = np.array(VE_list)
    accepted_mask = np.zeros(len(VEs), dtype=bool)
    accepted_mask[list_accepted_particles] = True
    
    accepted_VEs = VEs[accepted_mask]
    rejected_VEs = VEs[~accepted_mask]

    accept_lower, accept_upper = min(VE_range), max(VE_range)
    
    # Separate data
    data_to_plot = [accepted_VEs, rejected_VEs]
    labels = ['Accepted', 'Rejected']
    colors = ['royalblue', 'silver']
    
    fig, ax = plt.subplots(figsize=(5,5))
    bp = ax.boxplot(data_to_plot, patch_artist=True, labels=labels)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        # patch.set_zorder(2)
    
    ax.axhspan(accept_lower, accept_upper, color='orange', alpha=0.4, zorder=0, label='Acceptance range')
    
    rng = np.random.default_rng(246)                                            # rng only used to create horizontal jitter on the dots
    for i, y in enumerate(data_to_plot):
        jitter_accepted = rng.normal(1, 0.06, size=sum(accepted_mask))
        jitter_rejected = rng.normal(2, 0.06, size=sum(~accepted_mask))
        ax.scatter(jitter_accepted, VEs[accepted_mask], color='royalblue', alpha=0.4, s=5, label='Accepted points')
        ax.scatter(jitter_rejected, VEs[~accepted_mask], color='silver', alpha=0.4, s=5, label='Rejected points')
    
    ax.set_ylabel("Vaccine effectiveness (VE)")
    ax.set_title("Accepted vs rejected VE's")
    
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys())
    
    plt.tight_layout()
    fig.savefig(fig_filepath_VE_bp)
    plt.close(fig)

def vacc_post_inf_plot(list_vacc_post_inf, list_accepted_particles, fig_folder_path):
    
    fig_filepath_vacc_post_inf_bp = os.path.join(fig_folder_path, "boxplot_vacc_post_inf")
    
    vacc_post_inf = np.array(list_vacc_post_inf)
    accepted_mask = np.zeros(len(vacc_post_inf), dtype=bool)
    accepted_mask[list_accepted_particles] = True
    
    accepted_runs = vacc_post_inf[accepted_mask]
    rejected_runs = vacc_post_inf[~accepted_mask]
    
    # Separate data
    data_to_plot = [accepted_runs, rejected_runs]
    labels = ['Accepted', 'Rejected']
    colors = ['royalblue', 'silver']
    
    fig, ax = plt.subplots(figsize=(8,5))
    bp = ax.boxplot(data_to_plot, patch_artist=True, labels=labels)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        # patch.set_zorder(2)
    
    rng = np.random.default_rng(246)                                            # rng only used to create horizontal jitter on the dots
    for i, y in enumerate(data_to_plot):
        jitter_accepted = rng.normal(1, 0.06, size=sum(accepted_mask))
        jitter_rejected = rng.normal(2, 0.06, size=sum(~accepted_mask))
        ax.scatter(jitter_accepted, vacc_post_inf[accepted_mask], color='royalblue', alpha=0.4, s=5, label='Accepted points')
        ax.scatter(jitter_rejected, vacc_post_inf[~accepted_mask], color='silver', alpha=0.4, s=5, label='Rejected points')
    
    ax.set_ylabel("Vaccinated after infection")
    ax.set_title("Number of vaccinated after infection")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))

    
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys())
    
    plt.tight_layout()
    fig.savefig(fig_filepath_vacc_post_inf_bp)
    plt.close(fig)

def plot_R0(list_R0_acc, list_R0_rej, fig_folder_path):
    
    data = [list_R0_acc, list_R0_rej]
    labels = ["Accepted", "Rejected"]
    dot_colors = ["royalblue", "silver"]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Boxplots without facecolor
    bp = ax.boxplot(data, patch_artist=False, widths=0.4, positions=[1, 2])
    
    # Add points with small horizontal jitter
    jitter_strength = 0.05
    for i, (lst, color) in enumerate(zip(data, dot_colors), start=1):
        x = np.random.normal(i, jitter_strength, size=len(lst))
        ax.scatter(x, lst, color=color, alpha=0.7, s=20)
    
    for i, label in enumerate(labels, start=1):
       ax.text(i, ax.get_ylim()[1] + 0.05*(ax.get_ylim()[1]-ax.get_ylim()[0]), label, ha='center', va='bottom', fontsize=12)
    
    # Remove vertical tick marks (x-axis ticks)
    ax.set_xticks([])
    
    # Remove the box around the plot if desired
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # ax.set_ylim(0, 2)
    
    # Add suptitle
    fig.suptitle("R0 values for accepted and rejected particles", fontsize=16)
    
    # Adjust layout
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig_filepath_R0_plot = os.path.join(fig_folder_path, "R0_plot")
    fig.savefig(fig_filepath_R0_plot)
    plt.close(fig)

def VE_comparison_plot(VE_list, VE_alt_list, list_accepted_particles, AR_list, fig_folder_path):
    fig_filepath_comparison_VE = os.path.join(fig_folder_path, "VE_comparison")
    
    VE_list, VE_alt_list = np.array(VE_list), np.array(VE_alt_list)
    
    acc_list_for_plot = set(list_accepted_particles)
    AR_list = np.array(AR_list)
    
    VE_accepted = [VE_list[i] for i in range(len(VE_list)) if i in acc_list_for_plot]
    VE_alt_accepted = [VE_alt_list[i] for i in range(len(VE_list)) if i in acc_list_for_plot]
    AR_acc = [AR_list[i] for i in range(len(VE_list)) if i in acc_list_for_plot]
    
    VE_rejected = [VE_list[i] for i in range(len(VE_list)) if i not in acc_list_for_plot]
    VE_alt_rejected = [VE_alt_list[i] for i in range(len(VE_list)) if i not in acc_list_for_plot]
    AR_rej = [AR_list[i] for i in range(len(VE_list)) if i not in acc_list_for_plot]
    
    fig, ax = plt.subplots()
    
    sc1 = ax.scatter(VE_rejected, VE_alt_rejected, c = AR_rej, marker='^', s = 8, cmap= "plasma", label='Rejected', alpha = 0.5)
    sc2 = ax.scatter(VE_accepted, VE_alt_accepted, c = AR_acc , marker='o', cmap = "plasma" ,label='Accepted', alpha = 0.8)

    ax.set_xlabel('"True" VE')
    ax.set_ylabel("Approximate VE")
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.legend()
    
    cbar = fig.colorbar(sc2, ax=ax)
    cbar.set_label("AR value")
    
    fig_filepath_comparison_VE = os.path.join(fig_folder_path, "VE_comparison")
    fig.savefig(fig_filepath_comparison_VE)
    plt.close(fig)
    
def compute_spag_line(idx, incidences, I_comps):
    total_inf = np.sum([np.sum(incidences[c], axis=0) for c in I_comps], axis=0)
    cumul_inf = np.cumsum(total_inf)
    return idx, total_inf, cumul_inf

def plot_particle(part_id, incidences, *, nr_parts, nr_accepted, fig_folder_path, params, days, AR_range, log_scale = False, smooth_fig = False, perc_fig = False):
    
    nr_days, I_comps, t_start, t_end, ages = itemgetter("nr_days", "I_comps", "t_start", "t_end", "age_groups")(params)
    
    I_comps_greek = ["I^{S^0}", "I^{S^1}", "I^{V^0}", "I^{V^1}"]    # Used for the titles of the subplots
    palette_age = dict(zip(ages, sns.color_palette("deep")))
    
    thin_space_formatter = FuncFormatter(lambda x, _: f"{int(x):,}".replace(",", "\u2009"))
    
    fig_subfolder_path = os.path.join(fig_folder_path, f"Particle {part_id}")
    os.makedirs(fig_subfolder_path, exist_ok=True)
    fig_filepath_parent = os.path.join(fig_subfolder_path, f"Particle_{part_id}")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes=axes.flatten()
    
    for comp_id, comp in enumerate(I_comps):
        ax = axes[comp_id]
        ax.set_title(f"${I_comps_greek[comp_id]}$")
        
        for age_id, age in enumerate(ages):
            y = incidences[comp]
            if perc_fig:
                y = 100 * y / params["N"][:, None]    
            
            if smooth_fig:
                # The two variables below are used for the smoothing function
                window_length = 11  # must be odd
                polyorder = 2
                
                smoothed = savgol_filter(y[age_id, :], window_length, polyorder)
                ax.plot(days, smoothed, color=palette_age[age], linewidth=1.8)
            else:
                ax.plot(days, y[age_id,:], color=palette_age[age], label=age, linewidth=1.8)    
        
        if log_scale:
            ax.set_yscale("log")
            if perc_fig:
                ax.set_ylim(bottom=0.0001)
            else:
                ax.set_ylim(bottom=0.9)
        else:
            ax.set_ylim(bottom=0)
        
        if perc_fig:
            ax.yaxis.set_major_formatter(PercentFormatter(xmax = 1.0))
        else:
            ax.yaxis.set_major_formatter(thin_space_formatter)
        
        ax.minorticks_off()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axvspan(params["t_iiv_start"], params["t_iiv_end"], facecolor='yellow', alpha=0.2, edgecolor="none")
        ax.set_xlim(t_start, t_end)            
        ax.grid(False)
    
        if comp_id in [2, 3]:
            ax.set_xlabel("Day")

        if comp_id in [0, 2]:
            if log_scale and perc_fig:
                ax.set_ylabel("Daily Infection Incidence\n (% of age group)\n (log scale)", labelpad=10)
            elif log_scale and not perc_fig:
                ax.set_ylabel("Daily Infection Incidence\n (log scale)", labelpad=10)
            elif not log_scale and perc_fig:
                ax.set_ylabel("Daily Infection Incidence\n (% of age group)", labelpad=10)
            else:
                ax.set_ylabel("Daily Infection Incidence", labelpad=10)
    
    fig.align_ylabels()
    handles = [plt.Line2D([0], [0], color=palette_age[age], lw=1.8) for age in ages]
    labels = ages.copy()
    handles.append(mpatches.Patch(color='yellow', alpha=0.2))
    labels.append("Vaccine\nrollout")
    fig.suptitle(f"Infection incidence ({nr_accepted} out of {nr_parts} particles accepted)", y=0.97)

    fig_filepath = fig_filepath_parent + "_inf_per_comp.png"
    plt.tight_layout(rect=[0, 0, 0.87, 1])
    fig.legend(handles, labels,  loc='center left', bbox_to_anchor=(0.86, 0.5), frameon=False)
    
    fig.savefig(fig_filepath, dpi=300)
    plt.close(fig)
    
    # Figure of daily infections across all I states
    
    fig, axes = plt.subplots(figsize=(12, 8))
    
    tot_inf_per_age = np.sum([incidences[c] for c in I_comps], axis = 0)
    
    # Output df of incidences for use as baseline in the forward projections
    inc_df_filepath = fig_filepath_parent + "_inc_df.csv"
    part_df = pd.DataFrame(tot_inf_per_age.T, columns = ages)
    part_df.to_csv(inc_df_filepath, index = False)
    
    for age_id, age in enumerate(ages):
        axes.plot(days, tot_inf_per_age[age_id ,:], color=palette_age[age], label=age, linewidth=1.8)
        
    axes.minorticks_off()
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.axvspan(params["t_iiv_start"], params["t_iiv_end"], facecolor='yellow', alpha=0.2, edgecolor="none")
    axes.set_xlim(t_start, t_end)            
    axes.grid(False)
    axes.set_ylim(bottom=0)
    
    axes.set_xlabel('Day')
    axes.set_ylabel('Daily infections', labelpad=10)
    fig.suptitle("Total daily infection incidence", y=0.97)

    fig_filepath = fig_filepath_parent + "_total_daily_inf.png"
    plt.tight_layout(rect=[0, 0, 0.87, 1])
    fig.legend(handles, labels,  loc='center left', bbox_to_anchor=(0.86, 0.5), frameon=False)

    fig.savefig(fig_filepath, dpi=300)
    plt.close(fig)
    
    # Figure of cumulative I state time series
    I_cumul_per_age = np.cumsum(tot_inf_per_age, axis=1)
    I_cumul = np.sum(I_cumul_per_age, axis=0)
    
    fig, axes = plt.subplots(figsize=(12, 8))
    for age_id, age in enumerate(ages):
        axes.plot(days, I_cumul_per_age[age_id ,:], color=palette_age[age], label=age, linewidth=1.8)
    
    axes.plot(days, I_cumul, label="Total cumulative infection", linewidth=2.5, color="black")
    
    handles = [plt.Line2D([0], [0], color=palette_age[age], lw=1.8) for age in ages]
    labels = ages.copy()
    handles.append(plt.Line2D([0], [0], color="black", lw=2.5))
    labels.append("Total cumul.\ninfection")
    handles.append(mpatches.Patch(color='orange', alpha=0.2))
    labels.append("Acceptance\nrange")
        
    axes.minorticks_off()
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.axhspan(min(AR_range) * params["tot_pop"], max(AR_range) * params["tot_pop"], facecolor='orange', alpha=0.2, edgecolor="none")

    axes.set_xlim(t_start, t_end)            
    axes.grid(False)
    axes.set_ylim(bottom=0)
    axes.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    axes.yaxis.set_major_formatter(thin_space_formatter)
    
    axes.set_xlabel('Day')
    axes.set_ylabel('Cumulative infections', labelpad=10)
    fig.suptitle("Cumulative infections", y=0.94)

    fig_filepath = fig_filepath_parent + "_cumul_inf.png"
    plt.tight_layout(rect=[0, 0, 0.82, 1])
    
    fig.legend(handles, labels,  loc='center left', bbox_to_anchor=(0.82, 0.5), frameon=False)

    fig.savefig(fig_filepath, dpi=300)
    plt.close(fig)
    
    return part_id

def plot_selected_parts_spag(params, incidences, nr_parts, nr_runs_per_part, parallel, nr_cores, fig_folder_path, thin_space_formatter, list_selected_parts):
    
    nr_days, I_comps, t_start, t_end, ages = itemgetter("nr_days", "I_comps", "t_start", "t_end", "age_groups")(params)
    days = np.arange(nr_days)
    
    # Infection incidence spagetthi plot
    fig_filepath_inf_inc_selected_parts = os.path.join(fig_folder_path, "Spagetthi_particle_inf_inc_selected_parts")
    fig_filepath_cumul_inf_selected_parts = os.path.join(fig_folder_path, "Spagetthi_particle_cumul_inf_selected_parts")
    
    palette = list(plt.get_cmap("Dark2").colors)
    rng_col = np.random.default_rng(4328921)                                    # seed only used to select colours from palette
    selected_colours = rng_col.choice(palette, size = nr_parts, replace = False)
    
    if parallel:
        with ProcessPoolExecutor(max_workers=nr_cores) as executor:
            futures = [executor.submit(compute_spag_line, i, incidences[i], I_comps) for i in range(len(incidences))]
            results = [f.result() for f in as_completed(futures)]
    else:
        results = []
        for i, inc_dict in enumerate(incidences):
            # Total infectious incidence = sum over all infectious states and ages
            total_inf_inc = np.sum([np.sum(inc_dict[c], axis=0) for c in I_comps], axis=0)
            cumul_inf = np.cumsum(total_inf_inc)
            results.append((i, total_inf_inc, cumul_inf))
    
    results.sort(key=lambda x: x[0])
    
    fig_spag_inf_inc_selected_parts, ax_spag_inf_inc_selected_parts = plt.subplots(figsize=(12, 8))
    fig_spag_cumul_inf_selected_parts, ax_spag_cumul_inf_selected_parts = plt.subplots(figsize=(12, 8))
    
    for idx, total_inf, cumul_inf in results:
        
        particle_idx = idx // nr_runs_per_part
        colour = selected_colours[particle_idx]
        label = (f"Particle {list_selected_parts[idx // nr_runs_per_part]}" if idx % nr_runs_per_part == 0 else None)
        
        ax_spag_inf_inc_selected_parts.plot(
            days, total_inf, 
            color=colour,
            linewidth=0.8,
            alpha=0.2, 
            label=label
            )
        
        ax_spag_cumul_inf_selected_parts.plot(
            days, cumul_inf, 
            color=colour, 
            linewidth= 0.8, 
            alpha=0.2, 
            label=label
            )
    
    style_spaghetti_plot(params, ax_spag_inf_inc_selected_parts, cumul = False, fig = fig_spag_inf_inc_selected_parts, t_start = t_start, t_end = t_end, thin_space_formatter = thin_space_formatter, title="Infection incidence", ylab = "Daily Infection Incidence", selected_parts = True)
    style_spaghetti_plot(params, ax_spag_cumul_inf_selected_parts, cumul = True, fig = fig_spag_cumul_inf_selected_parts, t_start = t_start, t_end = t_end, thin_space_formatter = thin_space_formatter, title="Cumulative infection incidence", ylab = "Cumulative Infection Incidence", selected_parts = True)
    
    fig_spag_inf_inc_selected_parts.savefig(fig_filepath_inf_inc_selected_parts)
    fig_spag_cumul_inf_selected_parts.savefig(fig_filepath_cumul_inf_selected_parts)
    plt.close(fig_spag_inf_inc_selected_parts)
    plt.close(fig_spag_cumul_inf_selected_parts)

def outputs(prevalences, incidences, AR_list, VE_list, VE_alt_list, params, list_accepted_particles, list_vacc_post_inf, peak_time_range, AR_range, VE_range,
            full_parts, sample_keys, fig_spag_plot, fig_inf_per_comp_for_age_groups, fig_VE_comparison_plot, fig_corr_accepted_particles, fig_VE_plots, boxplot_vacc_post_inf, 
            sample_ranges, fig_folder_path, parallel, nr_cores, R0_plot, list_R0_acc, list_R0_rej, selected_spag_plot = False, list_selected_parts = []): 
    
    print("Plotting figures:\n")
    
    plt.close("all")                                                            # clear plt memory
    
    # nr_age, ages, nr_days, t_start, t_end, I_comps = itemgetter("nr_age", "age_groups", "nr_days", "t_start", "t_end" ,"I_comps")(params)
    nr_days, I_comps, t_start, t_end, ages = itemgetter("nr_days", "I_comps", "t_start", "t_end", "age_groups")(params)
    days = np.arange(nr_days)
    
    plt.rcParams.update({
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "legend.title_fontsize": 16,
    "figure.titlesize": 20
    })
    
    thin_space_formatter = FuncFormatter(lambda x, _: f"{int(x):,}".replace(",", "\u2009"))
    
    # =============================================================================
    #          Spagetthi plot particles
    # =============================================================================
    
    if fig_spag_plot:
        start_spag_fig = time.time()
        
        # Infection incidence spagetthi plot
        fig_filepath_inf_inc = os.path.join(fig_folder_path, "Spagetthi_particle_inf_inc")
        fig_filepath_cumul_inf = os.path.join(fig_folder_path, "Spagetthi_particle_cumul_inf")
        
        if parallel:
            with ProcessPoolExecutor(max_workers=nr_cores) as executor:
                futures = [executor.submit(compute_spag_line, i, incidences[i], I_comps) for i in range(len(incidences))]
                results = [f.result() for f in as_completed(futures)]
        else:
            results = []
            for i, inc_dict in enumerate(incidences):
                # Total infectious incidence = sum over all infectious states and ages
                total_inf_inc = np.sum([np.sum(inc_dict[c], axis=0) for c in I_comps], axis=0)
                cumul_inf = np.cumsum(total_inf_inc)
                results.append((i, total_inf_inc, cumul_inf))
                
        results.sort(key=lambda x: (x[0] in list_accepted_particles, x[0]))
        
        fig_spag_inf_inc, ax_spag_inf_inc = plt.subplots(figsize=(12, 8))
        fig_spag_cumul_inf, ax_spag_cumul_inf = plt.subplots(figsize=(12, 8))
        
        for idx, total_inf, cumul_inf in results:
            color = "royalblue" if idx in list_accepted_particles else "silver"
            linewidth = 1.3 if idx in list_accepted_particles else 0.5
            alpha = 0.6 if idx in list_accepted_particles else 0.3
            label = "Accepted particle" if idx == list_accepted_particles[0] else ""
            
            ax_spag_inf_inc.plot(days, total_inf, color=color, linewidth=linewidth, alpha=alpha, label=label)
            ax_spag_cumul_inf.plot(days, cumul_inf, color=color, linewidth=linewidth, alpha=alpha, label=label)
        
        style_spaghetti_plot(params, ax_spag_inf_inc, False, fig_spag_inf_inc, peak_time_range, t_start, t_end, thin_space_formatter, title=f"Infection incidence\n({len(list_accepted_particles)} out of {len(incidences)} particles accepted)", ylab = "Daily Infection Incidence", vert_or_horiz = "vert")
        style_spaghetti_plot(params, ax_spag_cumul_inf, True, fig_spag_cumul_inf, AR_range, t_start, t_end, thin_space_formatter, title=f"Cumulative infection incidence\n({len(list_accepted_particles)} out of {len(incidences)} particles accepted)", ylab = "Cumulative Infection Incidence", vert_or_horiz = "horiz")
        
        fig_spag_inf_inc.savefig(fig_filepath_inf_inc)
        fig_spag_cumul_inf.savefig(fig_filepath_cumul_inf)
        plt.close(fig_spag_inf_inc)
        plt.close(fig_spag_cumul_inf)
        
        print("Spagetthi plots done.\n")
        
        end_spag_fig = time.time()
        elapsed = end_spag_fig - start_spag_fig
        print(f"Spaghetti plot took {elapsed:.3f} seconds to produce.\n")
    
    # =============================================================================
    #          Correlation plot of accepted parameters
    # =============================================================================
    
    if fig_corr_accepted_particles:
        start_correl_fig = time.time()
        
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
    
        end_correl_fig = time.time()
        elapsed = end_correl_fig - start_correl_fig
        print(f"Correlation plot took {elapsed:.3f} seconds to produce.\n")
    
    # =============================================================================
    #          Infection incidence plots for accepted particles
    #          One plot for each I compartment for all age groups
    # =============================================================================
    
    if fig_inf_per_comp_for_age_groups:
        start_individual_fig = time.time()
        
        if parallel:
            
            plot_particle_prefilled = partial(plot_particle, nr_parts = len(incidences), nr_accepted = len(list_accepted_particles), fig_folder_path=fig_folder_path, days = days, AR_range = AR_range, params= params)
            
            with ProcessPoolExecutor(max_workers=nr_cores) as executor:
                list(tqdm(executor.map(plot_particle_prefilled, list_accepted_particles, [incidences[i] for i in list_accepted_particles], chunksize=1), total = len(list_accepted_particles), desc="Plotting particles"))
        else:
            
            iterator = 0 
            for part_id in list_accepted_particles:
                iterator += 1
                
                plot_particle(part_id, incidences[part_id], nr_parts=len(incidences), nr_accepted=len(list_accepted_particles), fig_folder_path = fig_folder_path, params=params, days=days, AR_range=AR_range)
                
                if iterator % 5 == 0:
                    print(f"Infection incidence for {iterator} out of {len(list_accepted_particles)} accepted particles plotted.")
                    print("")

            end_individual_fig = time.time()
            elapsed = end_individual_fig - start_individual_fig
            print(f"Individual plots took {elapsed:.3f} seconds to produce.\n")

    if fig_VE_plots:
        start_VE_fig = time.time()
        VE_plots(VE_list, VE_range, list_accepted_particles, fig_folder_path)
        print("VE boxplot done.\n")
        
        end_VE_fig = time.time()
        elapsed = end_VE_fig - start_VE_fig
        print(f"VE plot took {elapsed:.3f} seconds to produce.\n")
        
    # =============================================================================
    #          Boxplot of vaccination post infection
    # =============================================================================
    
    if boxplot_vacc_post_inf:
        start_boxplot_vacc_post_inf_fig = time.time()
        
        vacc_post_inf_plot(list_vacc_post_inf, list_accepted_particles, fig_folder_path)
        print("Boxplot vaccination post infection done.\n")
        
        end_boxplot_vacc_post_inf_fig = time.time()
        elapsed = end_boxplot_vacc_post_inf_fig - start_boxplot_vacc_post_inf_fig
        print(f"Boxplot of vaccination post infection took {elapsed:.3f} seconds to produce.\n")
    
    if fig_VE_comparison_plot:
        start_VE_comparison_fig = time.time()
        
        VE_comparison_plot(VE_list, VE_alt_list, list_accepted_particles, AR_list, fig_folder_path)
        print("Comparison plot between VE and alternative VE done.\n")
        
        end_VE_comparison_fig = time.time()
        elapsed = end_VE_comparison_fig - start_VE_comparison_fig
        print(f"Comparison plot between VE and alternative VE took {elapsed:.3f} seconds to produce.\n")
    
    if R0_plot:
        start_R0_plot = time.time()
        
        plot_R0(list_R0_acc, list_R0_rej, fig_folder_path)
        print("R0 boxplot done.\n")
        
        end_R0_plot = time.time()
        elapsed = end_R0_plot - start_R0_plot
        print(f"Boxplot of R0 took {elapsed:.3f} seconds to produce.\n")
        
    if selected_spag_plot:
        start_selected_spag_plot = time.time()
        
        nr_parts = len(list_accepted_particles)
        nr_runs_per_part = len(incidences) // nr_parts
        plot_selected_parts_spag(params, incidences, nr_parts, nr_runs_per_part, parallel, nr_cores, fig_folder_path, thin_space_formatter, list_selected_parts)
        print("Selected spaghetti plot done.\n")
        
        end_selected_spag_plot = time.time()
        elapsed = end_selected_spag_plot - start_selected_spag_plot
        print(f"Selected spaghetti plot took {elapsed:.3f} seconds to produce.\n")
    
    
if __name__ == "__main__":
    
    # Dummy values
    prevalences = 1
    incidences = 1
    params = 1
    list_accepted_particles = 1
    peak_time_range = 1
    AR_range = 1 
    full_parts = 1
    sample_keys = 1
    fig_spag_plot = 1
    fig_inf_per_comp_for_age_groups = 1
    fig_inf_per_age_group_for_comps = 1
    fig_corr_accepted_particles = 1
    
    outputs(prevalences, incidences, params, list_accepted_particles, peak_time_range, AR_range, full_parts, sample_keys, fig_spag_plot, fig_inf_per_comp_for_age_groups, fig_inf_per_age_group_for_comps, fig_corr_accepted_particles)
    
    plt.close("all")
        