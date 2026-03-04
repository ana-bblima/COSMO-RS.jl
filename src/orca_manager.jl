using Printf
using Dates

# Known ORCA CPCM solvent names (lowercase). Unknown solvents produce a warning.
const ORCA_CPCM_SOLVENTS = Set([
    "water", "acetonitrile", "acetone", "ammonia", "benzene",
    "ccl4", "ch2cl2", "chloroform", "cyclohexane", "dmf",
    "dmso", "ethanol", "ether", "hexane", "methanol",
    "octanol", "pyridine", "thf", "toluene", "dichloroethane",
])

"""
    OrcaSettings(method, basis, solvent, nprocs, maxcore; orca_path=nothing)

ORCA calculation settings. If `orca_path` is not provided, it is resolved
from the system PATH via `Sys.which("orca")`.
"""
struct OrcaSettings
    method::String
    basis::String
    solvent::String
    nprocs::Int
    maxcore::Int
    orca_path::String
end

function OrcaSettings(method::String, basis::String, solvent::String,
                      nprocs::Int, maxcore::Int;
                      orca_path::Union{String,Nothing} = nothing)
    if orca_path === nothing
        resolved = Sys.which("orca")
        resolved === nothing && error(
            "ORCA executable not found in PATH. " *
            "Either add it to PATH or pass orca_path explicitly to OrcaSettings.")
        orca_path = resolved
    end
    return OrcaSettings(method, basis, solvent, nprocs, maxcore, orca_path)
end


"""
    validate_orca_input(coords, atomic_numbers, charge, multiplicity, settings;
                        overlap_threshold=0.1)

Pre-flight validation of ORCA input parameters. Throws `OrcaInputError` on failure.

Checks:
1. Coordinate matrix dimensions match atom count
2. All coordinates are finite (no NaN/Inf)
3. No overlapping atoms (distance < `overlap_threshold` angstroms)
4. Charge/multiplicity electron count consistency
5. Solvent keyword against known ORCA CPCM solvents (warning only)
6. Settings sanity (nprocs >= 1, maxcore >= 100)
"""
function validate_orca_input(
        coords::Matrix{Float64},
        atomic_numbers::Vector{Int},
        charge::Int,
        multiplicity::Int,
        settings::OrcaSettings;
        overlap_threshold::Float64 = 0.1)

    natoms = length(atomic_numbers)

    # 1. Dimension consistency
    if size(coords) != (natoms, 3)
        throw(OrcaInputError(
            "Coordinate matrix has size $(size(coords)), expected ($natoms, 3)"))
    end

    # 2. Finite coordinates (no NaN/Inf)
    for i in 1:natoms
        for j in 1:3
            if !isfinite(coords[i, j])
                throw(OrcaInputError(
                    "Non-finite coordinate at atom $i, dimension $j: $(coords[i, j])"))
            end
        end
    end

    # 3. Overlapping atoms
    for i in 1:natoms
        for j in (i+1):natoms
            dx = coords[i, 1] - coords[j, 1]
            dy = coords[i, 2] - coords[j, 2]
            dz = coords[i, 3] - coords[j, 3]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            if dist < overlap_threshold
                throw(OrcaInputError(
                    "Atoms $i (Z=$(atomic_numbers[i])) and $j (Z=$(atomic_numbers[j])) " *
                    "overlap: distance = $(round(dist; digits=4)) A < threshold $overlap_threshold A"))
            end
        end
    end

    # 4. Charge/multiplicity electron count consistency
    total_electrons = sum(atomic_numbers) - charge
    if total_electrons <= 0
        throw(OrcaInputError(
            "Total electron count is $total_electrons " *
            "(sum(Z)=$(sum(atomic_numbers)), charge=$charge). This is non-physical."))
    end
    electrons_even = (total_electrons % 2 == 0)
    mult_even = (multiplicity % 2 == 0)
    # Even electrons -> odd multiplicity (1, 3, 5, ...)
    # Odd electrons  -> even multiplicity (2, 4, 6, ...)
    if electrons_even == mult_even
        throw(OrcaInputError(
            "Charge/multiplicity inconsistency: $total_electrons electrons " *
            "($(electrons_even ? "even" : "odd")) is incompatible with " *
            "multiplicity $multiplicity ($(mult_even ? "even" : "odd")). " *
            "For ground-state organics: even electrons -> mult 1, " *
            "odd electrons -> mult 2."))
    end

    # 5. Solvent keyword validation (warning, not error)
    if !isempty(settings.solvent)
        if !(lowercase(settings.solvent) in ORCA_CPCM_SOLVENTS)
            @warn "Solvent '$(settings.solvent)' not in known ORCA CPCM solvent list. " *
                  "ORCA may reject this. Known solvents: " *
                  join(sort(collect(ORCA_CPCM_SOLVENTS)), ", ")
        end
    end

    # 6. Settings sanity
    if settings.nprocs < 1
        throw(OrcaInputError("nprocs must be >= 1, got $(settings.nprocs)"))
    end
    if settings.maxcore < 100
        throw(OrcaInputError("maxcore must be >= 100 MB, got $(settings.maxcore)"))
    end

    return nothing
end


# Input file generator
function write_orca_input(
        filepath::String,
        coords::Matrix{Float64},
        atomic_numbers::Vector{Int},
        charge::Int,
        multiplicity::Int,
        settings::OrcaSettings;
        optimize_keyword::String = "")

    # Pre-flight validation
    validate_orca_input(coords, atomic_numbers, charge, multiplicity, settings)

    open(filepath, "w") do io

        println(io, "%maxcore $(settings.maxcore)")

        if settings.nprocs > 1
            println(io, "%pal nprocs $(settings.nprocs) end")
        end

        # Set base name for output files (.cpcm, .gbw, etc.)
        basename_no_ext = replace(basename(filepath), ".inp" => "")
        println(io, "%base \"$basename_no_ext\"")

        println(io)
        opt_str = isempty(optimize_keyword) ? "" : " $optimize_keyword"
        println(io,
            "! $(settings.method) $(settings.basis) CPCM($(settings.solvent)) TightSCF$opt_str")

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

    # Validate the generated file has balanced ORCA blocks
    validate_orca_blocks(filepath)
end

# Run orca
function run_orca_job(input_file::String, settings::OrcaSettings)
    job_dir = dirname(input_file)
    input_name = basename(input_file)
    output_name = replace(input_name, ".inp" => ".out")
    base_cmd = Cmd([settings.orca_path, input_name])
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
        base_dir::String = "orca_jobs",
        enable_copt_fallback::Bool = true)

    mkpath(base_dir)

    # Element symbols lookup for XYZ file generation
    element_symbols = Dict(
        1 => "H", 2 => "He", 3 => "Li", 4 => "Be", 5 => "B",
        6 => "C", 7 => "N", 8 => "O", 9 => "F", 10 => "Ne",
        11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si", 15 => "P",
        16 => "S", 17 => "Cl", 18 => "Ar", 19 => "K", 20 => "Ca",
        26 => "Fe", 29 => "Cu", 30 => "Zn", 35 => "Br", 53 => "I"
    )

    updated = Conformer[]

    for (i, conf) in enumerate(conformers)

        job_dir = joinpath(base_dir,
                           "$(state.name)_conf_$i")

        mkpath(job_dir)

        structname = "input"
        input_path = joinpath(job_dir, "$structname.inp")

        # Write XYZ file for .orcacosmo assembly (before ORCA run)
        xyzfile = joinpath(job_dir, "$structname.xyz")
        open(xyzfile, "w") do io
            natoms = length(conf.atomic_numbers)
            println(io, natoms)
            println(io, "$(state.name) conformer $i")
            for j in eachindex(conf.atomic_numbers)
                sym = get(element_symbols, conf.atomic_numbers[j],
                          string(conf.atomic_numbers[j]))
                x, y, z = conf.coordinates[j, :]
                @printf(io, "%-3s %14.8f %14.8f %14.8f\n", sym, x, y, z)
            end
        end

        # Run with retry logic; validation happens inside write_orca_input
        result = run_orca_with_retry(
            input_path, conf.coordinates, conf.atomic_numbers,
            state.charge, state.multiplicity, settings;
            enable_copt_fallback=enable_copt_fallback)

        # Build the .orcacosmo file from ORCA outputs
        method_str = "$(settings.method) $(settings.basis) CPCM($(settings.solvent))"
        sigma_file = build_orcacosmo(
            job_dir, structname, method_str,
            result.energy, result.dipole)

        # Convert Hartree -> kJ/mol
        energy = result.energy * KJMOL_PER_HARTREE

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

"""
    check_orca_success(logfile) -> OrcaResult

Parse ORCA output and return a structured `OrcaResult`.
Throws `OrcaTerminationError` if ORCA did not terminate normally.
Throws `OrcaConvergenceError` if SCF or geometry optimization did not converge.
"""
function check_orca_success(logfile::String)

    terminated_normally = false
    scf_not_converged = false
    opt_not_converged = false
    final_energy = nothing
    dipole = nothing

    for line in eachline(logfile)

        if occursin("****ORCA TERMINATED NORMALLY****", line)
            terminated_normally = true
        end

        if occursin("The optimization did not converge but reached", line)
            opt_not_converged = true
        end

        if occursin("SCF NOT CONVERGED", line)
            scf_not_converged = true
        end

        if occursin("FINAL SINGLE POINT ENERGY", line)
            parts = split(strip(line))
            final_energy = parse(Float64, parts[end])
        end

        if occursin("x,y,z [Debye]:", line)
            parts = split(strip(line))
            dipole = (
                parse(Float64, parts[end-2]),
                parse(Float64, parts[end-1]),
                parse(Float64, parts[end])
            )
        end
    end

    if !terminated_normally
        throw(OrcaTerminationError(
            "ORCA did not terminate normally", logfile))
    end

    if scf_not_converged
        throw(OrcaConvergenceError(
            "SCF did not converge", logfile))
    end

    if opt_not_converged
        throw(OrcaConvergenceError(
            "Geometry optimization did not converge", logfile))
    end

    if final_energy === nothing
        throw(OrcaTerminationError(
            "ORCA terminated but no 'FINAL SINGLE POINT ENERGY' found in output",
            logfile))
    end

    return OrcaResult(
        final_energy,
        dipole,
        terminated_normally,
        !scf_not_converged && !opt_not_converged,
        logfile
    )
end

function build_orcacosmo(
    job_dir::String,
    structname::String,
    method::String,
    energy::Float64,
    dipole
)

    outfile = joinpath(job_dir, "$structname.orcacosmo")
    xyzfile = joinpath(job_dir, "$structname.xyz")
    cpcmfile = joinpath(job_dir, "$structname.cpcm")

    open(outfile, "w") do io

        println(io, "$structname : $method")
        println(io)
        println(io, "#"^50)
        println(io, "#ENERGY")
        println(io, "FINAL SINGLE POINT ENERGY  $energy")

        if dipole !== nothing
            println(io)
            println(io, "#"^50)
            println(io, "#DIPOLE MOMENT (Debye)")
            println(io, "$(dipole[1]) $(dipole[2]) $(dipole[3])")
        end

        println(io)
        println(io, "#"^50)
        println(io, "#XYZ_FILE")

        for line in eachline(xyzfile)
            println(io, line)
        end

        if isfile(cpcmfile)
            println(io)
            println(io, "#"^50)
            println(io, "#COSMO")

            for line in eachline(cpcmfile)
                println(io, line)
            end
        end

        cpcm_corr_file = joinpath(job_dir, "$structname.cpcm_corr")
        if isfile(cpcm_corr_file)
            println(io)
            println(io, "#"^50)
            println(io, "#COSMO_corrected")

            for line in eachline(cpcm_corr_file)
                println(io, line)
            end
        end
    end

    return outfile
end


"""
    validate_orca_blocks(filepath)

Validate that `%block...end` pairs are balanced in an ORCA input file.
Throws `OrcaInputError` if unclosed blocks are found.

Single-line blocks (e.g., `%maxcore 2000` or `%pal nprocs 4 end`) are allowed.
Multi-line blocks (e.g., `%cpcm\\n...\\nend`) must have a matching `end`.
"""
function validate_orca_blocks(filepath::String)
    open_blocks = String[]

    for line in eachline(filepath)
        stripped = strip(line)

        if startswith(stripped, '%')
            block_name = split(stripped)[1]
            # Skip single-line directives: %maxcore, %base
            if block_name in ["%maxcore", "%base"]
                continue
            end
            # If the line contains "end", it's a single-line block (e.g., %pal nprocs 4 end)
            if !occursin("end", stripped)
                push!(open_blocks, block_name)
            end
        elseif lowercase(stripped) == "end" && !isempty(open_blocks)
            pop!(open_blocks)
        end
    end

    if !isempty(open_blocks)
        throw(OrcaInputError(
            "Unclosed ORCA blocks in $filepath: $(join(open_blocks, ", "))"))
    end

    return nothing
end


"""
    run_orca_with_retry(input_path, coords, atomic_numbers, charge, multiplicity,
                        settings; enable_copt_fallback=true) -> OrcaResult

Run ORCA job with optional COPT fallback on geometry optimization convergence failure.
Mirrors the Python `ORCA.execute()` retry pattern: if optimization fails with standard
settings, retry with Cartesian optimization (COPT).
"""
function run_orca_with_retry(
        input_path::String,
        coords::Matrix{Float64},
        atomic_numbers::Vector{Int},
        charge::Int,
        multiplicity::Int,
        settings::OrcaSettings;
        enable_copt_fallback::Bool = true)

    # First attempt: standard run
    write_orca_input(input_path, coords, atomic_numbers, charge, multiplicity, settings)
    run_orca_job(input_path, settings)

    output_file = replace(input_path, ".inp" => ".out")

    try
        return check_orca_success(output_file)
    catch e
        if e isa OrcaConvergenceError && enable_copt_fallback
            @warn "ORCA convergence failed, retrying with COPT (Cartesian optimization)..." logfile=e.logfile
            write_orca_input(input_path, coords, atomic_numbers,
                            charge, multiplicity, settings;
                            optimize_keyword="COPT")
            run_orca_job(input_path, settings)
            return check_orca_success(output_file)
        else
            rethrow()
        end
    end
end