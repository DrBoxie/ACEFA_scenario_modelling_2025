import gzip
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import main_sim as sim

curr_dir = os.getcwd().lower()
files_folder_path = os.path.join(curr_dir, "Forward projection")
filename_AR_pickle = os.path.join(files_folder_path, "multi_run_attack_rates.pkl")

with gzip.open(filename_AR_pickle, "rb") as f:
    AR_totals = pickle.load(f)

filename_particles = os.path.join(files_folder_path, "particles.csv")
filename_pes_phis = os.path.join(files_folder_path, "LAIV_pes_phis.csv")
filename_mid_phis = os.path.join(files_folder_path, "LAIV_mid_phis.csv")
filename_opt_phis = os.path.join(files_folder_path, "LAIV_opt_phis.csv")
filename_reduction_df = os.path.join(files_folder_path, "mean_reduction_df_scenario_F.csv")
particles = pd.read_csv(filename_particles, index_col = 0)
laiv_pes_phis = pd.read_csv(filename_pes_phis, index_col = 0)
laiv_mid_phis = pd.read_csv(filename_mid_phis, index_col = 0)
laiv_opt_phis = pd.read_csv(filename_opt_phis, index_col = 0)
particles_enriched = pd.concat([particles, laiv_pes_phis, laiv_mid_phis, laiv_opt_phis], axis=1)

params = sim.init_params()
scenarios = params["scenarios"]

overview_list = []
nr_runs = len(AR_totals[0]["A"])
for particle, scen_dict in AR_totals.items():
    for run in range(nr_runs):
        baseline_AR = scen_dict["A"][run]
        for scenario, value_list in scen_dict.items():
            if scenario == "A":
                continue
            
            attack_rate = value_list[run]
            
            pct_reduction = (
                100 * (baseline_AR - attack_rate) / baseline_AR
                if baseline_AR != 0 else float("nan")
            )
            
            overview_list.append({
                "particle": particle,
                "scenario": scenario,
                "run_nr": run,
                "baseline_AR": baseline_AR,
                "scenario_AR": attack_rate,
                "pct_reduction": pct_reduction
            })

df_reduction = pd.DataFrame(overview_list)
df_mean_reduction = (df_reduction.groupby(["particle", "scenario"], as_index = False)["pct_reduction"].mean())
df_mean_reduction["rank"] = df_mean_reduction.groupby("scenario")["pct_reduction"].rank(ascending = True, method = "min")

particles_lookup = particles_enriched.reset_index(names = "old_particle_nr")

df_mean_reduction = df_mean_reduction.merge(
    particles_lookup,
    left_on="particle",
    right_index=True,
    how="left"
)

rank_matrix = df_mean_reduction.pivot(index = "particle", columns = "scenario", values = "rank").sort_values(by = "F")
scenario_names = [scenarios[s][3] for s in rank_matrix.columns]

ranking_matrix_fig_path = os.path.join(files_folder_path, "ranking_matrix")

fig, ax = plt.subplots(figsize=(9, 35))

im = ax.imshow(
    rank_matrix.values,
    aspect="auto",
    cmap="RdYlBu_r"
)

ax.set_xticks(np.arange(len(rank_matrix.columns)))
ax.set_xticklabels(scenario_names, rotation = 45, ha = "right", fontsize = 20)

ax.set_yticks(np.arange(len(rank_matrix.index)))
ax.set_yticklabels(rank_matrix.index, fontsize = 10)

ax.set_xlabel("Scenario", fontsize = 25, labelpad = 20)
ax.set_ylabel("Particle", fontsize = 25, labelpad = 20)

# ax.tick_params(axis = "both", labelsize = 18)

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Rank", fontsize = 26, labelpad = 20)
cbar.ax.tick_params(labelsize = 20)

ax.set_title("Particle rankings across scenarios", fontsize = 34, pad = 30)

plt.tight_layout()
fig.savefig(ranking_matrix_fig_path, dpi=300, bbox_inches = "tight")
plt.close(fig)

#########################

ranking_correlation_fig_path = os.path.join(files_folder_path, "ranking_correlation")

corr_matrix = rank_matrix.corr(method="spearman")
scenario_names_corr = [scenarios[s][3] for s in corr_matrix.columns]

fig, ax = plt.subplots(figsize=(7, 6))

im = ax.imshow(
    corr_matrix.values,
    vmin=0,
    vmax=1,
    cmap="RdBu_r"
)

ax.set_xticks(np.arange(len(corr_matrix.columns)))
ax.set_xticklabels(scenario_names_corr, rotation = 45, ha = "right")

ax.set_yticks(np.arange(len(corr_matrix.index)))
ax.set_yticklabels(scenario_names_corr)

for i in range(len(corr_matrix)):
    for j in range(len(corr_matrix)):
        ax.text(
            j,
            i,
            f"{corr_matrix.iloc[i, j]:.2f}",
            ha="center",
            va="center",
            color="white",
            alpha = 0.1
        )

cbar = plt.colorbar(im, ax=ax)
cbar.set_label("Ranking correlation")

ax.set_title("Correlation of particle rankings")

plt.tight_layout()

fig.savefig(ranking_correlation_fig_path, dpi=300, bbox_inches = "tight")
plt.close(fig)
