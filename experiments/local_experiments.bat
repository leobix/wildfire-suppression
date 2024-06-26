julia --threads 12 --project=package_dependencies/julia experiments/profile.jl
julia --threads 12 --project=package_dependencies/julia experiments/cuts_at_root_node.jl
julia --threads 12 --project=package_dependencies/julia experiments/cglp_lazy_constraints.jl
julia --threads 12 --project=package_dependencies/julia experiments/branch_price_and_cut.jl
julia --threads 12 --project=package_dependencies/julia experiments/network_flow_direct.jl
julia --threads 12 --project=package_dependencies/julia experiments/triage_then_route.jl
python -m venv package_dependencies\python & package_dependencies\python\Scripts\Activate & pip install -r package_dependencies\python\requirements.txt & python experiments/process_outputs.py & deactivate