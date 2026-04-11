using Pkg
Pkg.activate(".")

include("../src/types.jl")
include("../src/states.jl")
include("../src/conformers.jl")
include("../src/orca_manager.jl")
include("../src/sigma_parser.jl")
include("../src/sigma_profile.jl")
include("../src/boltzmann.jl")
include("../src/mixture.jl")

Pkg.instantiate()
orca_path = "/home/analima/Programs/orca-6.1.0-f.0_linux_x86-64/bin/orca"

state = MolecularState("octane", "CCCCCCCC", 0, 1)
conf = generate_conformers(state; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)

result = run_cosmo_workflow(
    conf, state;
    base_dir  = "noctane_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result.conformers]

energies = [c.energy for c in result.conformers]

profile = average_sigma_profiles(profiles, energies)

using Plots

plot(profile.sigma_grid, profile.area_distribution, label="Sigma Profile", xlabel="Sigma (e/Angs²)", ylabel="Area Distribution (Angs²)", title="Sigma Profile of Octane Conformers")

savefig("octane_sigma_profile_5conf.png")

state = MolecularState("sorbitol", "OC[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)CO", 0, 1)
conf = generate_conformers(state; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
conf
result = run_cosmo_workflow(
    conf, state;
    base_dir  = "sorbitol_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result.conformers]
energies = [c.energy for c in result.conformers]
profile = average_sigma_profiles(profiles, energies)

plot(profile.sigma_grid, profile.area_distribution, label="Sigma Profile", xlabel="Sigma (e/Angs²)", ylabel="Area Distribution (Angs²)", title="Sigma Profile of Sorbitol Conformers")
savefig("sorbitol_sigma_profile_19conf.png")


states = generate_states("glycine")
ensembles = generate_conformers(states; nconfs=50)

ensembles = run_cosmo_workflow(ensembles;
    base_dir  = "glycine_cosmo",
    nprocs    = 5,
    maxcore   = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

total_profile = average_over_states(ensembles)

plot(total_profile.sigma_grid, total_profile.area_distribution, label="Sigma Profile", xlabel="Sigma (e/Angs²)", ylabel="Area Distribution (Angs²)", title="Sigma Profile of Glycine Conformers")
savefig("glycine_sigma_profile_50conf.png")



states = generate_states("his")

conf = generate_conformers(states[3]; nconfs=50, rms_threshold=0.5, energy_window_kJmol=6.0)

result = run_cosmo_workflow(
    conf, states[3];
    base_dir  = "his_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result.conformers]
energies = [c.energy for c in result.conformers]
profile = average_sigma_profiles(profiles, energies)    


profile_1 = parse_sigma_surface(result.conformers[1].sigma_surface_file)
profile_2 = parse_sigma_surface(result.conformers[2].sigma_surface_file)
profile_3 = parse_sigma_surface(result.conformers[3].sigma_surface_file)


surface_1 = compute_sigma_profile(profile_1)
surface_2 = compute_sigma_profile(profile_2)
surface_3 = compute_sigma_profile(profile_3)

plot(surface_1.sigma_grid, surface_1.area_distribution, label="Conformer 1", xlabel="Sigma (e/Angs²)", ylabel="Area Distribution (Angs²)", title="Sigma Profiles of His Conformers")
plot!(surface_2.sigma_grid, surface_2.area_distribution, label="Conformer 2")
plot!(surface_3.sigma_grid, surface_3.area_distribution, label="Conformer 3")
plot!(profile.sigma_grid, profile.area_distribution, label="Boltzmann Averaged Profile", xlabel="Sigma (e/Angs²)", ylabel="Area Distribution (Angs²)", title="Sigma Profiles of His Conformers", linewidth=3, color=:black)

savefig("his_sigma_profile_50conf_zwitterion.png")

conf = generate_conformers(states[1:2]; nconfs=50, rms_threshold=0.5, energy_window_kJmol=6.0)

result_2 = run_cosmo_workflow(
    conf; 
    base_dir  = "his_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

ensemble = average_over_states(result_2)

plot(ensemble.sigma_grid, ensemble.area_distribution, label="Sigma Profile", xlabel="Sigma (e/Angs²)", 
    ylabel="Area Distribution (Angs²)", title="Sigma Profile of His neutral states")
savefig("his_sigma_profile_50conf_all_states.png")

profile = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_2.conformers]
energies = [c.energy for c in result_2.conformers]
surface = average_sigma_profiles(profile, energies)



states = generate_states("tryptophan")

conf = generate_conformers(states[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)

results_3 = run_cosmo_workflow(
    conf, states[2];
    base_dir  = "trp_cosmo_zwitterion",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in results_3.conformers]
energies = [c.energy for c in results_3.conformers]

states = generate_states("diglycine")
conf = generate_conformers(states[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)

results_4 = run_cosmo_workflow(
    conf, states[2];
    base_dir  = "diglycine_cosmo_zwitterion",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)