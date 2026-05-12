
import gzip
import pickle
import os
import matplotlib.pyplot as plt
import pandas as pd

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

particles_lookup = particles_enriched.reset_index(names = "old_particle_nr")

df_mean_reduction = df_mean_reduction.merge(
    particles_lookup,
    left_on="particle",
    right_index=True,
    how="left"
)

# only keep the mid 5-18 scenario, which corresponds to scenario F
df_only_mid_5_18 = df_mean_reduction[df_mean_reduction["scenario"] == "F"].sort_values(by = "pct_reduction")

df_only_mid_5_18.to_csv(filename_reduction_df, index = True)

# select the rows for the 4 particles that will be analysed
most_neg_part = df_only_mid_5_18.loc[df_only_mid_5_18["pct_reduction"].idxmin()]

zero_change_part = df_only_mid_5_18.loc[(df_only_mid_5_18["pct_reduction"].abs()).idxmin()]

n = len(df_only_mid_5_18)

if n % 2 == 1:
    median_part = df_only_mid_5_18.iloc[n // 2]
else:
    # lower of the two middle values if there are an even number of particles in total
    median_part = df_only_mid_5_18.iloc[n // 2 - 1]

max_reduction_part = df_only_mid_5_18.loc[df_only_mid_5_18["pct_reduction"].idxmax()]

selected_rows = [
    ("Most negative reduction", most_neg_part),
    ("No change", zero_change_part),
    ("Median reduction", median_part),
    ("Maximal reduction", max_reduction_part),
]

table_rows = []

for reason, row in selected_rows:

    old_particle_id = row["old_particle_nr"]
    new_particle_id = row["particle"]
    pct_red = round(row["pct_reduction"], 2)

    particle_vals = particles_enriched.loc[old_particle_id]

    out_row = {
        "Reason": reason,
        "Particle ID (new)": new_particle_id,
        "Particle ID (old)": old_particle_id,
        "Infection incidence reduction": f"{pct_red:.2f}",
    }

    # add the particle values (after rounding)
    formatted_particle_vals = {}
    for key, value in particle_vals.to_dict().items():
        formatted_particle_vals[key] = f"{value:.2f}"
    out_row.update(formatted_particle_vals)

    table_rows.append(out_row)

summary_table = pd.DataFrame(table_rows)

# -------------------------------------------------
# Plot as figure-table
# -------------------------------------------------

table_filename = os.path.join(files_folder_path, "overview_table")

col_headers = summary_table["Particle ID (new)"].astype(str).tolist()
summary_table = summary_table.drop(columns = ["Particle ID (new)"]).T

row_order = [
    "Particle ID (old)",
    "Reason",
    "Infection incidence reduction",
    "beta_0",
    "beta_1",
    "shift",
    "frac_S1",
    "phi_S1",
    "phi_V0",
    "LAIV_phi_V0_pes",
    "LAIV_phi_V0_mid",
    "LAIV_phi_V0_opt",
    "phi_V1",
    "LAIV_phi_V1_pes",
    "LAIV_phi_V1_mid",
    "LAIV_phi_V1_opt"
    ]

summary_table = summary_table.loc[row_order]
# summary_table.loc[""] = [""] * summary_table.shape[1]


reduction_row = "Infection incidence reduction"
summary_table.loc[reduction_row] = summary_table.loc[reduction_row].astype(str) + " %"

fig, ax = plt.subplots(figsize=(18, 12))
ax.axis("off")

fig.suptitle("Selected particles: Central scenario 5-18 compared to Baseline", fontsize = 24, fontweight = "bold", x = 0.37, y = 0.80)

tbl = ax.table(
    cellText=summary_table.values,
    rowLabels=summary_table.index,
    colLabels = col_headers,
    cellLoc="center",
    loc="center"
)

tbl.add_cell(
    0, -1,
    width=tbl[(1, -1)].get_width(),
    height=tbl[(0, 0)].get_height(),
    text="Particle ID",
    loc="left"
)

tbl[(0, -1)].set_text_props(weight="bold", fontsize=20)

tbl.auto_set_font_size(False)
tbl.set_fontsize(18)
tbl.scale(1.2, 2.2)


tbl[(0, -1)].set_text_props(weight = "bold", fontsize = 20)
for col in range(len(col_headers)):
    tbl[(0, col)].set_text_props(weight="bold", fontsize=20)

shaded_rows = {0, 1, 2 ,3, 7, 9, 10, 11, 12}

for row in range(0, len(summary_table.index)+1):
    if row in shaded_rows:
        for col in range(-1, len(col_headers)):
            tbl[(row, col)].set_facecolor("#f0f0f0")
            
boundaries = []            

max_row = len(summary_table.index)

for row in range(0, max_row):
    if (row in shaded_rows) != ((row+1) in shaded_rows):
        boundaries.append(row)

fig.canvas.draw()

eps = 0.002

for row in boundaries:
    row_idx = row

    left = tbl[(row_idx, -1)]
    right = tbl[(row_idx, len(col_headers)-1)]

    x0 = left.get_x() + eps
    x1 = right.get_x() + right.get_width() - eps
    y = left.get_y()

    ax.plot(
        [x0, x1],
        [y, y],
        transform=ax.transAxes,
        color="black",
        linewidth=3,
        clip_on = False
    )
    
# plt.tight_layout()
plt.savefig(table_filename, dpi = 300, bbox_inches="tight")
plt.close(fig)


    
