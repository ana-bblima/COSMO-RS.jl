using Pkg
Pkg.activate(".")
using Plots

include("../src/types.jl")
include("../src/states.jl")
include("../src/conformers.jl")
include("../src/orca_manager.jl")
include("../src/sigma_parser.jl")
include("../src/sigma_profile.jl")
include("../src/boltzmann.jl")
include("../src/mixture.jl")


orca_path = "/home/analima/Programs/orca-6.1.0-f.0_linux_x86-64/bin/orca"

st_dgly = generate_states("diglycine")
cf_dgly = generate_conformers(st_dgly[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
#
result_dgly = run_cosmo_workflow(
    cf_dgly, st_dgly[2];
    base_dir  = "dgly_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

dgly_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_dgly.conformers]
dgly_energies = [c.energy for c in result_dgly.conformers]
d_gly_final_sigma_profile = average_sigma_profiles(dgly_profiles, dgly_energies)

plot(d_gly_final_sigma_profile.sigma_grid, 
    d_gly_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Diglycine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)
savefig("dgly_sigma_profile_3conf.png")

st_tgly = generate_states("triglycine")
cf_tgly = generate_conformers(st_tgly[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
result_tgly = run_cosmo_workflow(
    cf_tgly, st_tgly[2];
    base_dir  = "tgly_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

tgly_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_tgly.conformers]
tgly_energies = [c.energy for c in result_tgly.conformers]
t_gly_final_sigma_profile = average_sigma_profiles(tgly_profiles, tgly_energies)

plot(t_gly_final_sigma_profile.sigma_grid, 
    t_gly_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Triglycine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)
savefig("tgly_sigma_profile_3conf.png")

st_ala = generate_states("alanine")
cf_ala = generate_conformers(st_ala[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
result_ala = run_cosmo_workflow(
    cf_ala, st_ala[2];
    base_dir  = "ala_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

ala_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_ala.conformers]
ala_energies = [c.energy for c in result_ala.conformers]
ala_final_sigma_profile = average_sigma_profiles(ala_profiles, ala_energies)

plot(ala_final_sigma_profile.sigma_grid, 
    ala_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Alanine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)
savefig("ala_sigma_profile_3conf.png")

st_val = generate_states("valine")
cf_val = generate_conformers(st_val[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
result_val = run_cosmo_workflow(
    cf_val, st_val[2];
    base_dir  = "val_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

val_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_val.conformers]
val_energies = [c.energy for c in result_val.conformers]
val_final_sigma_profile = average_sigma_profiles(val_profiles, val_energies)

plot(val_final_sigma_profile.sigma_grid, 
    val_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Valine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)
savefig("val_sigma_profile_3conf.png")


st_ser = generate_states("serine")
cf_ser = generate_conformers(st_ser[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
result_ser = run_cosmo_workflow(
    cf_ser, st_ser[2];
    base_dir  = "ser_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

ser_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_ser.conformers]
ser_energies = [c.energy for c in result_ser.conformers]
ser_final_sigma_profile = average_sigma_profiles(ser_profiles, ser_energies)


plot(ser_final_sigma_profile.sigma_grid, 
    ser_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Serine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)

savefig("ser_sigma_profile_3conf.png")


st_thr = generate_states("threonine")
cf_thr = generate_conformers(st_thr[2]; nconfs=100, rms_threshold=0.5, energy_window_kJmol=6.0)
result_thr = run_cosmo_workflow(
    cf_thr, st_thr[2];
    base_dir  = "thr_cosmo",
    nprocs = 5,
    maxcore = 2000,
    orca_path = orca_path,
    max_cpcm_conformers_after_fast = 3
)

thr_profiles = [compute_sigma_profile(parse_sigma_surface(c.sigma_surface_file)) for c in result_thr.conformers]
thr_energies = [c.energy for c in result_thr.conformers]
thr_final_sigma_profile = average_sigma_profiles(thr_profiles, thr_energies)

plot(thr_final_sigma_profile.sigma_grid, 
    thr_final_sigma_profile.area_distribution, 
    label="Sigma Profile", xlabel="Sigma (e/Angs²)",
     ylabel="Area Distribution (Angs²)", title="Sigma Profile of Threonine with 3 Conformers",
     legend=false,
     linewidth=2,
     color=:blue)

savefig("thr_sigma_profile_3conf.png")