# ACEFA Scenario Modelling 2025

This repository contains the code for the ACEFA scenario-modelling analysis.

## Execution pipeline

The full execution pipeline begins with `abc_code.py`.

Before running the ABC analysis, review and adapt the following variables as required:

| Variable | Description |
|---|---|
| `ABC_output_results` | Activates or suppresses output of data frames and plots from the ABC runs. |
| `seed` | Seed used to generate the seed list supplied to the individual simulation runs. |
| `seed2` | Secondary random seed. |
| `parallel` | Determines whether simulations are run in parallel (`True`) or sequentially (`False`). |
| `nr_cores` | Number of CPU cores used for parallel processing. Ignored when `parallel = False`. |
| `nr_particles` | Number of ABC particles to run. |
| `peak_time_range` | Acceptance range for the timing of the epidemic peak. |
| `AR_range` | Acceptance range for the attack rate (AR). |
| `VE_range` | Acceptance range for vaccine efficacy (VE). |

## ABC parameter sampling

The sampling ranges for the ABC particles are defined in `list_to_sample_params`.

The following parameter ranges can be specified:

- `beta_0_range`
- `beta_1_range`
- `shift`
- `half-life`
- `phi_S1_range`
- `phi_V0_range`
- `phi_V1_range`
- `frac_S1_range`

These ranges determine the parameter space from which the ABC particles are sampled.

## Model parameters

Core model variables are defined in the `init_params` function in `main_sim.py`. If the underlying model structure or assumptions are changed, the corresponding parameters should be updated there.

| Variable | Description |
|---|---|
| `tau` | Size of the tau-leap time step. |
| `t_start`, `t_end` | Start and end days of the simulation. |
| `age_groups` | Age groups included in the model. |
| `compartments` | Epidemiological and vaccination compartments included in the model. |
| `tot_pop` | Total population size. |
| `frac_ages` | Fraction of the population belonging to each age group. |
| `C` | Age-stratified contact matrix. |
| `gamma` | Recovery rate for each age group. |
| `target_cover` | Target vaccination coverage for each age group. |
| `iiv_rates` | Weekly cumulative rollout percentages for inactivated influenza vaccine (IIV). |
| `laiv_rates` | Weekly cumulative rollout percentages for live attenuated influenza vaccine (LAIV). |
| `p_ext` | External force of infection. |

## Recommended workflow

1. Review the ABC settings and acceptance criteria in `abc_code.py`.
2. Review the ABC parameter sampling ranges in `list_to_sample_params`.
3. If the model structure or fixed parameters have changed, update `init_params` in `main_sim.py`.
4. Run `abc_code.py`.
5. Inspect the accepted ABC particles and generated outputs before proceeding to subsequent scenario analyses.
