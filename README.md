# Fire Suppression Crew Routing

This repository contains tools to run empirical experiments for routing fire suppression crews to wildfire incidents. The instructions below explain how to set up the Julia environment and launch a new experiment from raw data.

## 1. Clone the repository and checkout the empirical branch

```bash
git clone https://github.com/jacobwachspress/fire-suppression-crew-routing.git
cd fire-suppression-crew-routing
git checkout empirical-fire-model
```

## 2. Set up the Julia environment

The project uses Julia and the JuMP ecosystem with Gurobi as the solver.

1. Install Julia 1.10 or newer.
2. Ensure Gurobi is installed and a valid license is available. On most systems this means setting `GUROBI_HOME` and `GRB_LICENSE_FILE` appropriately.
3. From the repository root, run

```bash
julia --project=package_dependencies/julia
```

4. In the Julia REPL, instantiate and precompile the dependencies:

```julia
julia> using Pkg
julia> Pkg.instantiate()
julia> Pkg.precompile()
```

5. Exit the Julia session.

## 3. Running an experiment

From the repository root, execute:

```bash
julia --project=package_dependencies/julia EmpiricalMain.jl
```

An optional `--debug` flag can be passed to expose verbose logging:

```bash
julia --project=package_dependencies/julia EmpiricalMain.jl --debug
```

Use `--firefighters-per-crew` to control the number of personnel assigned to each crew (default `70`):

```bash
julia --project=package_dependencies/julia EmpiricalMain.jl --firefighters-per-crew 20
```

The run produces JSON files describing crew and fire arcs in `data/output/` and writes `selected_fires_sorted.csv` to `data/empirical_fire_models/raw/arc_arrays/` for visualization.

## 4. Preparing data for new case studies

All raw data files live in `data/empirical_fire_models/raw/arc_arrays/`. To run a new experiment, replace or augment the following files:

* **Fire and base distances**
  * `fire_fire_distances.csv` – pairwise distances between fires.
  * `base_fire_distances.csv` – distances from each crew base to each fire.
* **Selected fires**
  * `selected_fires.csv` – list of fires to include. After a run, the script produces `selected_fires_sorted.csv` in the same folder.
* **Arc descriptions and costs**
  * `arc_arrays_*.csv` – one file per fire containing the state transition arcs. The filename should follow `arc_arrays_<NIFCID>_<FIRENAME>_day0.csv`.
  * `arc_costs_*.csv` – cost of each arc, matched to the corresponding `arc_arrays_*` file: `arc_costs_<NIFCID>_<FIRENAME>_day0.csv`.

Place each new CSV in `data/empirical_fire_models/raw/arc_arrays/` using the same naming pattern. The program automatically loads every `arc_arrays_*.csv` and `arc_costs_*.csv` present.

## 5. Tweaking experiment parameters

Edit [`EmpiricalMain.jl`](EmpiricalMain.jl) to modify run parameters:

* `num_fires`, `num_crews`, and `num_time_periods` control the size of the case study.
* `travel_speed = 40.0 * 6.0` encodes a 40 mph average speed for 6 hours of travel per day. Change the second factor to adjust allowed daily travel time.
* Pass `--firefighters-per-crew` to set the number of personnel per crew (default `70`).

Save the file and rerun the script to evaluate the new settings.

## 6. Output

Each invocation writes arc information for every fire and crew to JSON files in `data/output/` with filenames of the form `fire_arcs_<fire>_<time>.json`, `fire_arc_costs_<fire>_<time>.json`, `crew_arcs_<crew>_<time>.json`, and `crew_arc_costs_<crew>_<time>.json`.

## 7. Repeating experiments

To run another case study, swap in a new set of CSV inputs (keeping the naming conventions), adjust parameters as needed, and repeat the run command in Section 3.

