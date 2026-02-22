using Printf
using Dates

struct OrcaSettings
    method::String
    basis::String
    solvent::String
    nprocs::Int
    maxcore::Int
end


# Input file generator
function write_orca_input(
        filepath::String,
        coords::Matrix{Float64},
        atomic_numbers::Vector{Int},
        charge::Int,
        multiplicity::Int,
        settings::OrcaSettings)

    open(filepath, "w") do io

        println(io, "%maxcore $(settings.maxcore)")

        if settings.nprocs > 1
            println(io, "%pal nprocs $(settings.nprocs) end")
        end

        println(io)
        println(io,
            "! $(settings.method) $(settings.basis) CPCM($(settings.solvent)) TightSCF")

        # ---- THIS IS CRITICAL ----
        println(io)
        println(io, "%cpcm")
        println(io, "   writemos true")
        println(io, "end")

        println(io)
        println(io, "* xyz $charge $multiplicity")

        for i in eachindex(atomic_numbers)
            Z = atomic_numbers[i]
            x, y, z = coords[i, :]
            @printf(io, "%-3d %14.8f %14.8f %14.8f\n",
                Z, x, y, z)
        end

        println(io, "*")
    end
end

# Run orca
function run_orca_job(input_file::String)

    path_to_orca = "/home/analima/Downloads/orca-6.1.0-f.0_linux_x86-64_openmpi41/orca-6.1.0-f.0_linux_x86-64/bin/orca"

    job_dir = dirname(input_file)
    input_name = basename(input_file)

    output_name = replace(input_name, ".inp" => ".out")

    base_cmd = Cmd([path_to_orca, input_name])
    cmd = Cmd(base_cmd; dir=job_dir)

    open(joinpath(job_dir, output_name), "w") do io
        run(pipeline(cmd, stdout=io, stderr=io))
    end
end

# Energy Extraction
function parse_orca_energy(output_file::String)

    energy = nothing

    for line in eachline(output_file)
        if occursin("FINAL SINGLE POINT ENERGY", line)
            parts = split(line)
            energy = parse(Float64, parts[end])
        end
    end

    energy === nothing && error("Energy not found")

    return energy  # Hartree
end


function find_orcacosmo_file(directory::String)
    files = readdir(directory)
    for f in files
        if endswith(f, ".orcacosmo")
            return joinpath(directory, f)
        end
    end
    error("No .orcacosmo file found")
end

function run_orca_for_conformers(
        conformers::Vector{Conformer},
        state::MolecularState,
        settings::OrcaSettings;
        base_dir::String = "orca_jobs")

    mkpath(base_dir)

    updated = Conformer[]

    for (i, conf) in enumerate(conformers)

        job_dir = joinpath(base_dir,
                           "$(state.name)_conf_$i")

        mkpath(job_dir)

        input_path = joinpath(job_dir, "input.inp")

        write_orca_input(
            input_path,
            conf.coordinates,
            conf.atomic_numbers,
            state.charge,
            state.multiplicity,
            settings)

        run_orca_job(input_path)

        output_file = replace(input_path, ".inp" => ".out")

        energy_hartree =
            parse_orca_energy(output_file)

        # Convert Hartree → J/mol
        energy =
            energy_hartree * 2625.49962 * 1000.0

        sigma_file =
            find_orcacosmo_file(job_dir)

        push!(updated,
            Conformer(
                conf.coordinates,
                conf.atomic_numbers,
                energy,
                sigma_file
            )
        )
    end

    return updated
end


