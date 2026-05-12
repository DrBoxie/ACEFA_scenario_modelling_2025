
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import main_sim as sim
from operator import itemgetter
from matplotlib.ticker import FuncFormatter
import time
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import duckdb

plt.close("all")

# if "new" than plots are for the outputs after updates of phi sampling, if "old" than it's for those before the updates
old_vs_new_plots = "new"

def plot_reproduction_calcs(old_vs_new_plots):
    
    curr_dir = os.getcwd().lower()
    files_folder_path = os.path.join(curr_dir, "Forward projection")
    filename_full_df = os.path.join(files_folder_path, "df_incid_summed.parquet")
    filename_max_timing = os.path.join(files_folder_path, "df_max_timing.parquet")
    df_max_timing = pd.read_parquet(filename_max_timing)
    filename_df_reduction_pct = os.path.join(files_folder_path, "mean_reduction_df_scenario_F.csv")
    df_red_only_mid_5_18 = pd.read_csv(filename_df_reduction_pct, index_col = 0)
    
    if old_vs_new_plots == "old":
        selected_parts = [55, 57, 67, 40]
        col_map = {
            55: "most_neg",
            57: "no_change",
            67: "med_red",
            40: "max_red"
            }
    elif old_vs_new_plots == "new":
        selected_parts = [20, 27, 55, 141]
        col_map = {
            20: "most_neg",
            27: "no_change",
            55: "med_red",
            141: "max_red"
            }
    else:
        raise ValueError("Error: old_vs_new_plots value can only be 'old' or 'new' !.")
    
    selected_parts_str = ", ".join(map(str, selected_parts))
    df_selected = duckdb.sql(f"""
        SELECT *
        FROM read_parquet('{filename_full_df}')
        WHERE simulation_index IN ({selected_parts_str})
    """).df()
    
    params = sim.init_params()
    nr_days, t_start, t_end, scenarios, age_groups = itemgetter("nr_days", "t_start", "t_end", "scenarios", "age_groups")(params)
    
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
    
    param_ranges = {
        "beta_0": (0.037, 0.055),
        "beta_1": (0.05, 0.3),
        "shift": (40, 190),
        "frac_S1": (0.2, 0.8),
        "phi_S1": (0.2, 0.8),
        "phi_V0": (0.2, 0.8),
        "phi_V1": (0.2, 0.8),
        "phi_V0_pes": (0, 0.8),
        "phi_V1_pes": (0, 0.8),
        "phi_V0_mid": (0, 0.8),
        "phi_V1_mid": (0, 0.8),
        "phi_V0_opt": (0, 0.8),
        "phi_V1_opt": (0, 0.8),
        }
    
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
    
    # plot_spaghetti_per_scenario(files_folder_path, thin_space_formatter, days, df_selected, colors, col_map, particle_labels, legend_handles, scenario_nms, scenario_labels)
    
    # plot_spaghetti_age_strat(files_folder_path, thin_space_formatter, days, df_selected, scenario_nms, scenario_labels, col_map, colors, particle_labels, age_groups)
    
    # plot_AR_swarm(files_folder_path, scenarios, df_selected, colors, col_map, particle_labels, thin_space_formatter, legend_handles)
    
    # plot_AR_age_strat_swarm(files_folder_path, scenarios, df_selected, colors, col_map, thin_space_formatter, legend_handles, age_groups)
    
    # plot_particles_correlation_enriched(files_folder_path, df_red_only_mid_5_18)        
    
    plot_shift_vs_peak_time(files_folder_path, df_max_timing, df_red_only_mid_5_18, thin_space_formatter, scenario_labels, scenarios, param_ranges)
    
    # plot_beta_vs_peak_time(files_folder_path, filename_full_df, df_red_only_mid_5_18, thin_space_formatter, scenario_labels, scenarios, param_ranges)
    
    # plot_scatter_plot_params(files_folder_path, df_red_only_mid_5_18, param_ranges)
    
    # plot_selected_betas(files_folder_path, df_red_only_mid_5_18, selected_parts, col_map, colors)

def plot_shift_vs_peak_time(output_folder_path, df_max_timing, df_parts, thin_space_formatter, scenario_labels, scenarios, param_ranges):
    
    print("Plotting shift vs infection incidence peak time.\n")
    
    nr_parts = len(df_parts)
    
    filename_shift_vs_peak = "shift_vs_peak_time"
    filename_beta_peak_vs_peak_time = "beta_peak_vs_peak_time"
    
    fig_filepath_shift_peak = os.path.join(output_folder_path, filename_shift_vs_peak)
    fig_filepath_beta_peak_vs_peak_time = os.path.join(output_folder_path, filename_beta_peak_vs_peak_time)
    
    nr_rows, nr_cols = 4, 3
    fig_shift_vs_peak, axes = plt.subplots(nr_rows, nr_cols, figsize=(20, 22), sharex=True)
    fig_beta_peak_vs_peak_time, axes_beta = plt.subplots(nr_rows, nr_cols, figsize=(20, 22), sharex=True)
    axes = axes.flatten() 
    axes_beta = axes_beta.flatten()
    axes[0].axis("off")  # left of baseline
    axes[2].axis("off")  # right of baseline
    axes_beta[0].axis("off")  # left of baseline
    axes_beta[2].axis("off")  # right of baseline
    
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
    
    ax_beta_map ={
        "A": axes_beta[1], # 0 and 2 are skipped because the top row only contains the baseline figure
        "B": axes_beta[3],       
        "C": axes_beta[4],        
        "D": axes_beta[5],        
        "E": axes_beta[6],       
        "F": axes_beta[7],        
        "G": axes_beta[8],     
        "H": axes_beta[9],
        "J": axes_beta[10],
        "K": axes_beta[11]
        }
    
    df_max_by_part = {
        part: group
        for part, group in df_max_timing.groupby("simulation_index")
    }
    
    fmt_shift = FuncFormatter(lambda x, _: f"{x:.0f}")
    fmt_beta_peak = FuncFormatter(lambda x, _: f"{x:.0f}")
    fmt_peak = FuncFormatter(lambda y, _: f"{y:.0f}")
    
    strip_height = 0.16
    
    shift_range = [0, 365] # use this range instead of the actual range of shift to make the figure into squares
    beta_peak_range = [0, 365] # use this range instead of the actual range of beta to make the figure into squares
    peak_day_range = [0, 365]

    # pct_norm = Normalize(
    #     vmin=df_parts["pct_reduction"].min(),
    #     vmax=df_parts["pct_reduction"].max()
    # )
    # cmap_pct = "coolwarm_r"

    pct_min = df_parts["pct_reduction"].min()
    pct_max = df_parts["pct_reduction"].max()
    
    # pct_norm = TwoSlopeNorm(
    #     vmin=pct_min,
    #     vcenter=0,
    #     vmax=pct_max
    # )
    
    pct_norm = Normalize(
        vmin = pct_min,
        vmax = pct_max
        )
    
    zero_pos = (0 - pct_min) / (pct_max - pct_min)
    
    cmap_pct = LinearSegmentedColormap.from_list(
        "red_gray_blue_green",
        [
            (0.00, "#7F0000"),   # strongly negative: dark red
            (zero_pos, "#B0B0B0"),   # zero: gray
            (min(zero_pos + 0.01, 1), "#2166AC"), # slightly positive: blue
            (1.00, "#1A9850"),   # strongly positive: green
        ]
    )
    
    tick_label_size = 22
    axis_label_size= 24
    
    for key in ax_map:
        ax = ax_map[key]
        ax_beta = ax_beta_map[key]
        
        ax.set_xlim(peak_day_range)
        ax.set_ylim(shift_range)
        
        ax_beta.set_xlim(peak_day_range)
        ax_beta.set_ylim(beta_peak_range)
        
        ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)
        ax_beta.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)
        
        ax.xaxis.set_major_formatter(fmt_peak)
        ax.yaxis.set_major_formatter(fmt_shift)
        
        ax_beta.xaxis.set_major_formatter(fmt_peak)
        ax_beta.yaxis.set_major_formatter(fmt_beta_peak)
        
        ax.add_patch(
            Rectangle(
                (0, 1.0),
                1,
                strip_height,
                transform=ax.transAxes,
                facecolor="white",
                edgecolor="black",
                linewidth=1.2,
                clip_on=False
            )
        )
        
        ax_beta.add_patch(
            Rectangle(
                (0, 1.0),
                1,
                strip_height,
                transform=ax_beta.transAxes,
                facecolor="white",
                edgecolor="black",
                linewidth=1.2,
                clip_on=False
            )
        )
        
        ax.text(
            0.5,
            1.0 + strip_height / 2,
            scenario_labels[key],
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
            clip_on=False
        )
        
        ax_beta.text(
            0.5,
            1.0 + strip_height / 2,
            scenario_labels[key],
            transform=ax_beta.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
            clip_on=False
        )
        
        ax.set_box_aspect(1)
        # ax.set_aspect("equal", adjustable="box")
        ax.plot(
        [0, 365],
        [0, 365],
        color="black",
        linestyle="--",
        linewidth=1
    )
        
        ax_beta.set_box_aspect(1)
        # ax.set_aspect("equal", adjustable="box")
        ax_beta.plot(
        [0, 365],
        [0, 365],
        color="black",
        linestyle="--",
        linewidth=1
    )
    
    for part in range(nr_parts):
        
        particle_params = df_parts.loc[df_parts["particle"] == part].iloc[0]

        shift = particle_params["shift"]
        
        beta_0 = particle_params["beta_0"]
        beta_1 = particle_params["beta_1"]
        t_start = 0
        t_end = 365
        nr_days = t_end - t_start
        t = np.linspace(t_start, t_end, num = nr_days)
        beta = beta_0 * (1 +  beta_1 * np.sin(( 2 * np.pi * (t - shift) ) / nr_days))
        beta_max_day = np.argmax(beta) + 1
        
        pct_reduction = particle_params["pct_reduction"]
    
        df_max_part = df_max_by_part[part]
        
        for key in ax_map:
        
            ax = ax_map[key]
            ax_beta = ax_beta_map[key]
    
            peak_days = (
                df_max_part[df_max_part["scenario"] == scenarios[key][4]]
                .sort_values("run_nr")
            )
    
            ax.scatter(
                peak_days["peak_day"].to_numpy(),
                np.full(len(peak_days), shift),
                c=np.full(len(peak_days), pct_reduction),
                cmap=cmap_pct,
                norm=pct_norm,
                s=7,
                alpha = 1,
            )
            
            ax_beta.scatter(
                peak_days["peak_day"].to_numpy(),
                np.full(len(peak_days), beta_max_day),
                c=np.full(len(peak_days), pct_reduction),
                cmap=cmap_pct,
                norm=pct_norm,
                s=7,
                alpha = 1,
            )
    
        if part % 10 == 0:
            print(f"Particle {part} plotted peak day plots.\n")
    
    for i in range(len(axes)):
    
        ax = axes[i]
        ax_beta = axes_beta[i]
    
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=tick_label_size
        )
        
        ax_beta.tick_params(
            axis="both",
            which="major",
            labelsize=tick_label_size
        )
    
        row = i // nr_cols
    
        if not ax.has_data():
            continue
    
        if row < nr_rows - 1:
            ax.set_xlabel("")
            ax_beta.set_xlabel("")
        else:
            ax.set_xlabel("peak day", fontsize=axis_label_size)
            ax_beta.set_xlabel("peak day", fontsize=axis_label_size)
    
        leftmost_ax = None
    
        for j in range(row * nr_cols, (row + 1) * nr_cols):
            if axes[j].has_data():
                leftmost_ax = axes[j]
                break
    
        if ax != leftmost_ax:
            ax.set_ylabel("")
            ax_beta.set_ylabel("")
        elif row in [1, 3]:
            ax.set_ylabel("")
            ax_beta.set_ylabel("")
        else:
            ax.set_ylabel("shift", fontsize=axis_label_size)
            ax_beta.set_ylabel("Beta peak day", fontsize=axis_label_size)
    
    sm = plt.cm.ScalarMappable(
        norm=pct_norm,
        cmap=cmap_pct
    )
    sm.set_array([])
    
    cbar_ax = fig_shift_vs_peak.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig_shift_vs_peak.colorbar(sm, cax=cbar_ax)
    
    cbar_beta_ax = fig_beta_peak_vs_peak_time.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar_beta = fig_beta_peak_vs_peak_time.colorbar(sm, cax=cbar_beta_ax)
    
    cbar.set_label("Percent reduction", fontsize=axis_label_size)
    cbar.ax.tick_params(labelsize=tick_label_size)
    
    cbar_beta.set_label("Percent reduction", fontsize=axis_label_size)
    cbar_beta.ax.tick_params(labelsize=tick_label_size)
    
    fig_shift_vs_peak.suptitle(
        "Peak time vs shift",
        fontsize=28,
        y=0.99
    )
    
    fig_beta_peak_vs_peak_time.suptitle(
        "Peak time vs beta peak day",
        fontsize=28,
        y=0.99
    )
    
    fig_shift_vs_peak.subplots_adjust(
        left=0.07,
        right=0.88,
        bottom=0.08,
        top=0.93,
        wspace=0.08,
        hspace=0.20
    )
    
    fig_beta_peak_vs_peak_time.subplots_adjust(
        left=0.07,
        right=0.88,
        bottom=0.08,
        top=0.93,
        wspace=0.08,
        hspace=0.20
    )
    
    fig_shift_vs_peak.savefig(fig_filepath_shift_peak, dpi=300)
    plt.close(fig_shift_vs_peak)
    fig_beta_peak_vs_peak_time.savefig(fig_filepath_beta_peak_vs_peak_time, dpi=300)
    plt.close(fig_beta_peak_vs_peak_time)
    
    print("Finished plotting beta vs infection incidence peak time.\n")

def plot_scatter_plot_params(output_folder_path, full_df, param_ranges):
    
    print("Plotting parameter scatter plots")
    
    fig_filepath = os.path.join(output_folder_path, "params_scatter_shared_axis")
    
    pct_red = full_df[["particle", "pct_reduction"]].copy()
    df = full_df.drop(columns = ["scenario", "old_particle_nr", "pct_reduction"])
    df.columns = df.columns.str.replace("LAIV_", "", regex = False)
    df = df.set_index("particle")
    
    fig, ax_base = plt.subplots(figsize = (18, 8))
    
    x = np.arange(len(param_ranges))
    
    ax_base.set_xticks(x)
    ax_base.set_xticklabels(df.columns, rotation = 45, ha = "right")
    
    ax_base.set_xlim(-0.5, len(param_ranges) - 0.5)
    ax_base.set_yticks([])
    
    axes = []
    
    for i, p in enumerate(df.columns):

        ax = ax_base.twinx()
    
        # Move axis to the correct x-position
        ax.spines["right"].set_position(("axes", i / len(param_ranges)))
    
        # Hide everything except the spine
        ax.spines["left"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
    
        low, high = param_ranges[p]
        ax.set_ylim(low, high)
    
        ax.set_ylabel(p, rotation=90)
    
        axes.append(ax)
    
    # -------------------------------------------------
    # Plot each row
    # -------------------------------------------------
    
    for _, row in df.iterrows():
    
        for i, p in enumerate(df.columns):
    
            axes[i].plot(
                [x[i]],
                [row[p]],
                marker="o",
                markersize=3,
                alpha=0.3,
                color="black"
            )
    
    # -------------------------------------------------
    # Final touches
    # -------------------------------------------------
    
    ax_base.set_title("Parameters on their native scales (multi-axis plot)")
    
    fig.tight_layout()
    
    fig.savefig(fig_filepath, dpi=300)
    plt.close(fig)
    
    print("Finished plotting parameter scatter plots")

def plot_selected_betas(output_folder_path, df, selected_parts, col_map, colors):
    
    print("Plotting selected betas.\n")
    
    fig_filepath = os.path.join(output_folder_path, "betas_plot")
    
    t_start = 0
    t_end = 365
    nr_days = t_end - t_start
    
    t = np.linspace(t_start, t_end, num = nr_days)
    
    # title_map = {
    #     55: "Most negative reduction",
    #     57: "No reduction",
    #     67: "Median reduction",
    #     40: "Maximal reduction"
    # }
    
    title_map = {
        20: "Most negative reduction",
        27: "No reduction",
        55: "Median reduction",
        141: "Maximal reduction"
    }
    
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    
    for ax, part in zip(axes, selected_parts):
        
        part_params = df.loc[df["particle"] == part].iloc[0]
        beta_0 = part_params["beta_0"]
        beta_1 = part_params["beta_1"]
        shift = part_params["shift"]
    
        beta = beta_0 * (1 +  beta_1 * np.sin(( 2 * np.pi * (t - shift) ) / nr_days))
        
        color_key = col_map[part]
        color = colors[color_key]
        
        ax.plot(t, beta, linewidth=2, color = color)

        ax.set_ylabel(r"$\beta(t)$")
    
        ax.set_title(
            rf"{title_map[part]}: "
            rf"$\beta_0={beta_0:.3f}$, "
            rf"$\beta_1={beta_1:.3f}$, "
            rf"shift$={shift:.1f}$"
        )
    
    axes[-1].set_xlabel("Time")
    
    fig.suptitle(r"Seasonal transmission rate $\beta(t)$", fontsize=16)
    
    fig.tight_layout()
    
    fig.savefig(fig_filepath)
    plt.close(fig)
    
    print("Finished plotting selected betas.\n")

def plot_beta_vs_peak_time(output_folder_path, filename_full_df, df_parts, thin_space_formatter, scenario_labels, scenarios, param_ranges):
    
    print("Plotting beta vs infection incidence peak time.\n")
    
    nr_parts = duckdb.sql(f"""
        SELECT MAX(simulation_index)
        FROM read_parquet('{filename_full_df}')
    """).fetchone()[0] + 1
    
    nr_runs = duckdb.sql(f"""
        SELECT MAX(run_nr)
        FROM read_parquet('{filename_full_df}')
    """).fetchone()[0] + 1
    
    filename_beta_peak = "beta_vs_peak_time"
    
    fig_filepath_beta_peak = os.path.join(output_folder_path, filename_beta_peak)
    
    nr_rows, nr_cols = 4, 3
    fig_beta_vs_peak, axes = plt.subplots(nr_rows, nr_cols, figsize=(23, 16), sharex=True)
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
    
    scenario_ids_sql = ", ".join(
        f"'{scenarios[key][4]}'"
        for key in ax_map.keys()
    )
    
    df_max = duckdb.sql(f"""
        WITH daily AS (
            SELECT
                simulation_index,
                scenario,
                run_nr,
                horizon,
                SUM(value) AS daily_value
            FROM read_parquet('{filename_full_df}')
            WHERE scenario IN ({scenario_ids_sql})
            GROUP BY
                simulation_index,
                scenario,
                run_nr,
                horizon
        )
        SELECT
            simulation_index,
            scenario,
            run_nr,
            arg_max(horizon, daily_value) AS peak_day
        FROM daily
        GROUP BY
            simulation_index,
            scenario,
            run_nr
    """).df()
    
    df_max_by_part = {
        part: group
        for part, group in df_max.groupby("simulation_index")
    }
    
    fmt_beta_0 = FuncFormatter(lambda x, _: f"{x:.3f}")
    fmt_beta_1 = FuncFormatter(lambda y, _: f"{y:.3f}")
    
    strip_height = 0.16
    
    peak_norm = Normalize(vmin = 0, vmax = 365)
    cmap_peak = "OrRd"
    
    beta_0_range = param_ranges["beta_0"]
    beta_1_range = param_ranges["beta_1"]
    
    width_beta_0 = 0.0005
    width_beta_1 = 0.0005
    rng = np.random.default_rng(seed=36558)
    
    jitter_beta_0 = rng.uniform(-width_beta_0, width_beta_0, nr_runs)
    jitter_beta_1 = rng.uniform(-width_beta_1, width_beta_1, nr_runs)
    
    tick_label_size = 22
    axis_label_size= 24
    
    axes[0].axis("off")  # left of baseline
    axes[2].axis("off")  # right of baseline
    
    for key, ax in ax_map.items():
        ax.set_xlim(beta_0_range)
        ax.set_ylim(beta_1_range)
        
        ax.xaxis.set_major_formatter(fmt_beta_0)
        ax.yaxis.set_major_formatter(fmt_beta_1)
        
        ax.add_patch(
            Rectangle(
                (0, 1.0),
                1,
                strip_height,
                transform=ax.transAxes,
                facecolor="white",
                edgecolor="black",
                linewidth=1.2,
                clip_on=False
            )
        )
        
        ax.text(
            0.5,
            1.0 + strip_height / 2,
            scenario_labels[key],
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            fontweight="bold",
            clip_on=False
        )
    
    for part in range(nr_parts):
        
        particle_params = df_parts.loc[df_parts["particle"] == part].drop(columns = ["scenario", "old_particle_nr"]).iloc[0]
        
        jittered_beta_0 = particle_params["beta_0"] + jitter_beta_0
        jittered_beta_1 = particle_params["beta_1"] + jitter_beta_1
        
        df_max_part = df_max_by_part.get(part)
        
        if df_max_part is None:
            continue
        
        for key, ax in ax_map.items():
            
            peak_days = df_max_part[df_max_part["scenario"] == scenarios[key][4]].sort_values("run_nr")
            
            ax.scatter(
                jittered_beta_0,
                jittered_beta_1,
                c = peak_days["horizon"].to_numpy(),
                cmap = cmap_peak,
                norm = peak_norm,
                s = 7,
                alpha = 0.8
                )
        if part % 10 == 0:
            print(f"Particle {part} plotted.\n")
        
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
            ax.set_xlabel(r"$\beta_0$", fontsize = axis_label_size)
        
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
            ax.set_ylabel(r"$\beta_1$", fontsize = axis_label_size)
    
    sm = plt.cm.ScalarMappable(norm=peak_norm, cmap=cmap_peak)
    sm.set_array([])       
    
    cbar_ax = fig_beta_vs_peak.add_axes([0.95, 0.15, 0.02, 0.7])
    cbar = fig_beta_vs_peak.colorbar(sm, cax = cbar_ax)
    
    cbar.set_label("Peak day", fontsize = axis_label_size)
    cbar.ax.tick_params(labelsize = tick_label_size)
    
    fig_beta_vs_peak.suptitle("Beta vs peak time", fontsize = 28, y = 0.99)
    fig_beta_vs_peak.tight_layout(rect=[0,0,0.93,0.95])
    fig_beta_vs_peak.savefig(fig_filepath_beta_peak)
    plt.close(fig_beta_vs_peak)
    
    print("Finished plotting beta vs infection incidence peak time.\n")
    
def plot_particles_correlation_enriched(fig_folder_path, df_red):
    
    print("Plotting correlation plot.\n")
    
    fig_filepath = os.path.join(fig_folder_path, "corr_parts_params_scenarios")
    
    df = df_red.drop(columns = ["particle", "old_particle_nr", "scenario", "pct_reduction"])
    
    sample_ranges = {
        "beta_0": [0.037, 0.055],
        "beta_1": [0.05, 0.3],
        "shift": [40, 190],
        "phi_S1": [0.2, 0.8], 
        "phi_V0": [0.2, 0.8], 
        "phi_V1": [0.2, 0.8],
        "frac_S1": [0.2, 0.8],
        "LAIV_phi_V0_pes": [0, 0.8],
        "LAIV_phi_V1_pes": [0, 0.8],
        "LAIV_phi_V0_mid": [0, 0.8],
        "LAIV_phi_V1_mid": [0, 0.8],
        "LAIV_phi_V0_opt": [0, 0.8],
        "LAIV_phi_V1_opt": [0, 0.8]
        }    
    
    phi_subscript_cols = [
        "LAIV_phi_V0_pes",
        "LAIV_phi_V1_pes",
        "LAIV_phi_V0_mid",
        "LAIV_phi_V1_mid",
        "LAIV_phi_V0_opt",
        "LAIV_phi_V1_opt"
    ]
    
    plot_labs = [
        r"$\beta_0$",
        r"$\beta_1$", 
        r"fraction $S^1$", 
        r"$\phi^{S^1}$", 
        r"$\phi^{V^0}$",
        r"$\phi^{V^1}$", 
        "Shift", 
        r"$\phi_{0}^{pes}$", 
        r"$\phi_{1}^{pes}$",
        r"$\phi_{0}^{mid}$", 
        r"$\phi_{1}^{mid}$", 
        r"$\phi_{0}^{opt}$", 
        r"$\phi_{1}^{opt}$"
        ]
    
    n = len(plot_labs)
    
    pct_norm = Normalize(
        vmin=df_red["pct_reduction"].min(),
        vmax=df_red["pct_reduction"].max()
    )
    
    cmap_pct = "coolwarm_r"
    
    # customise visualisation to make points with negative pct_reduction values stand out
    pct_red = df_red["pct_reduction"].to_numpy()
    neg_strength = np.clip(-pct_red, 0, None)
    if neg_strength.max() > 0:
        neg_strength = neg_strength / neg_strength.max()
    dot_sizes = 60 + 80 * neg_strength
    dot_alphas = 1
    
    with plt.rc_context({
       "axes.titlesize": 46,
       "axes.labelsize": 44,
       "xtick.labelsize": 38,
       "ytick.labelsize": 38,
       "legend.fontsize": 40,
       "legend.title_fontsize": 42,
       "figure.titlesize": 60
       }):
    
        fig, axes = plt.subplots(n, n, figsize=(4*n,4*n))        
        
        fmt_1 = FuncFormatter(lambda x, _: f"{x:.1f}")
        fmt_2 = FuncFormatter(lambda x, _: f"{x:.2f}")
        fmt_int = FuncFormatter(lambda x, _: f"{int(x)}")
        
        for i, y in enumerate(df):
            for j, x in enumerate(df):
                ax = axes[i,j]
                if i < j:
                    ax.set_visible(False)
                    continue
                elif i == j:                    
                    bins = np.linspace(min(sample_ranges[x]), max(sample_ranges[x]), 21)
                    ax.hist(
                        df[x],
                        bins=bins,
                        color = "royalblue",
                        alpha=0.7,
                        edgecolor="none",
                        density=True
                    )
                else:
                    ax.scatter(df[x], df[y], c = df_red["pct_reduction"], cmap = cmap_pct, s = dot_sizes, alpha=dot_alphas, edgecolors='none')
                
                ax.set_xlim(min(sample_ranges[x]), max(sample_ranges[x]))
                
                if x in phi_subscript_cols:
                    ax.set_xticks([0.0, 0.4, 0.8])
                
                # axes labels
                if i < n-1: 
                    ax.set_xticklabels([])
                else: 
                    ax.set_xlabel(plot_labs[j], labelpad=10)
                    
                    if i == n-1 and j == 0:
                        # bottom-left: 2 decimals
                        ax.xaxis.set_major_formatter(fmt_2)
                    
                    elif x == "shift":
                        ax.xaxis.set_major_formatter(fmt_int)
                        ax.set_xticks([50, 100, 150])
                
                    elif i == n-1 and j == n-1:
                        # bottom-right: integers only
                        ax.xaxis.set_major_formatter(fmt_1)
                
                    else:
                        # other bottom-row panels
                        ax.xaxis.set_major_formatter(fmt_1)
                
                ax.tick_params(axis = "x", labelrotation = 55)
                
                if j > 0:
                    ax.set_yticklabels([])
                else: 
                    ax.set_ylabel(plot_labs[i], labelpad=10)
                
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
        
        sm = plt.cm.ScalarMappable(
            norm=pct_norm,
            cmap=cmap_pct
        )
        sm.set_array([])
        
        cbar_left = 0.91
        cbar_width = 0.035
        cbar_height = 0.40
        cbar_bottom = (1 - cbar_height) / 2
        # Add a colorbar axis in the empty upper triangle
        # [left, bottom, width, height] in figure coordinates
        cax = fig.add_axes([cbar_left, cbar_bottom, cbar_width, cbar_height])  
        
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label("Percentage reduction in infections", fontsize = 44)
        cbar.ax.tick_params(labelsize = 38)
        
        plt.suptitle("Correlations of parameters — coloured by percentage reduction in infections", x=0.5, y=0.98)
        fig.align_xlabels()
        fig.align_ylabels()
        fig.subplots_adjust(left=0.08, right=0.96, bottom=0.05, top=0.95, hspace=0.35, wspace=0.35)        
        fig.savefig(fig_filepath)
        plt.close(fig)
        
        print("Finished plotting_correlation plot.\n")

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
        ax.set_ylabel(str(age), fontsize=axis_label_size, weight = "bold")
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
    
    nr_runs = df["run_nr"].nunique()
    
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
            ha = "left",
            weight = "bold"
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
    
    fig_age_scen_spag, axes = plt.subplots(nr_age, nr_scens, figsize=(24, 32), sharex=True, sharey = "row")
    
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
            ax.grid(True)
            
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
    
    plot_reproduction_calcs(old_vs_new_plots)
    
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Elapsed time: {elapsed:.3f} seconds\n")
    print("Plotting of special particles done. Yippie-ki-yay!!!\n")