using Printf
using Dates

# Known ORCA CPCM solvent names (lowercase). Unknown solvents produce a warning.
const ORCA_CPCM_SOLVENTS = Set([
    "water", "acetonitrile", "acetone", "ammonia", "benzene",
    "ccl4", "ch2cl2", "chloroform", "cyclohexane", "dmf",
    "dmso", "ethanol", "ether", "hexane", "methanol",
    "octanol", "pyridine", "thf", "toluene", "dichloroethane",
])

# Default CPCM radii for COSMO-RS calculations (atomic number => radius in bohr)
const DEFAULT_CPCM_RADII = Dict{Int,Float64}(
    1  => 1.3,    # H
    2  => 1.601,  # He
    5  => 2.048,  # B
    6  => 2.0,    # C
    7  => 1.83,   # N
    8  => 1.72,   # O
    9  => 1.72,   # F
    10 => 1.771,  # Ne
    13 => 2.152,  # Al
    14 => 2.2,    # Si
    15 => 2.106,  # P
    16 => 2.16,   # S
    17 => 2.05,   # Cl
    18 => 2.184,  # Ar
    22 => 2.261,  # Ti
    26 => 2.195,  # Fe
    27 => 2.153,  # Co
    31 => 2.172,  # Ga
    32 => 2.7,    # Ge
    33 => 2.148,  # As
    34 => 2.2,    # Se
    35 => 2.16,   # Br
    36 => 2.354,  # Kr
    49 => 2.245,  # In
    50 => 2.55,   # Sn
    51 => 2.402,  # Sb
    53 => 2.32,   # I
    54 => 2.524,  # Xe
    80 => 1.978,  # Hg
    82 => 2.36,   # Pb
    86 => 2.573,  # Rn
)

# Element symbols lookup for XYZ file generation
const ELEMENT_SYMBOLS = Dict(
    1 => "H", 2 => "He", 3 => "Li", 4 => "Be", 5 => "B",
    6 => "C", 7 => "N", 8 => "O", 9 => "F", 10 => "Ne",
    11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si", 15 => "P",
    16 => "S", 17 => "Cl", 18 => "Ar", 19 => "K", 20 => "Ca",
    26 => "Fe", 29 => "Cu", 30 => "Zn", 35 => "Br", 53 => "I"
)

"""
    write_xyz_file(filepath, coords, atomic_numbers, comment="")

Write an XYZ file from coordinates and atomic numbers.
"""
function write_xyz_file(filepath::String, coords::Matrix{Float64},
                        atomic_numbers::Vector{Int}, comment::String = "")
    open(filepath, "w") do io
        natoms = length(atomic_numbers)
        println(io, natoms)
        println(io, comment)
        for j in eachindex(atomic_numbers)
            sym = get(ELEMENT_SYMBOLS, atomic_numbers[j], string(atomic_numbers[j]))
            x, y, z = coords[j, :]
            @printf(io, "%-3s %14.8f %14.8f %14.8f\n", sym, x, y, z)
        end
    end
end

"""
    parse_xyz_file(filepath) -> (coords, elements)

Parse an XYZ file and return coordinates matrix (natoms × 3) and element symbols.
"""
function parse_xyz_file(filepath::String)
    lines = readlines(filepath)
    natoms = parse(Int, strip(lines[1]))
    coords = zeros(Float64, natoms, 3)
    elements = String[]
    for i in 1:natoms
        parts = split(strip(lines[i + 2]))
        push!(elements, parts[1])
        coords[i, 1] = parse(Float64, parts[2])
        coords[i, 2] = parse(Float64, parts[3])
        coords[i, 3] = parse(Float64, parts[4])
    end
    return coords, elements
end

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
function run_orca_job(input_file::String, orca_path::String)
    job_dir = dirname(input_file)
    input_name = basename(input_file)
    output_name = replace(input_name, ".inp" => ".out")
    base_cmd = Cmd([orca_path, input_name])
    cmd = Cmd(base_cmd; dir=job_dir)
    open(joinpath(job_dir, output_name), "w") do io
        run(pipeline(cmd, stdout=io, stderr=io))
    end
end

run_orca_job(input_file::String, settings::OrcaSettings) =
    run_orca_job(input_file, settings.orca_path)


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

    updated = Conformer[]

    for (i, conf) in enumerate(conformers)

        job_dir = joinpath(base_dir,
                           "$(state.name)_conf_$i")

        mkpath(job_dir)

        structname = "input"
        input_path = joinpath(job_dir, "$structname.inp")

        # Write XYZ file for .orcacosmo assembly (before ORCA run)
        xyzfile = joinpath(job_dir, "$structname.xyz")
        write_xyz_file(xyzfile, conf.coordinates, conf.atomic_numbers,
                       "$(state.name) conformer $i")

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
            block_name = lowercase(split(stripped)[1])
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


# ================================================================
# GAS-PHASE & CPCM WORKFLOW
# ================================================================
#
# The functions below implement the full COSMO-RS quantum-chemistry
# workflow mirroring the Python ConformerGenerator pipeline:
#
#   Gas phase:
#     1. ORCA_DFT_FAST    — OPT BP86/def2-TZVP(-f) (gas)
#     2. ORCA_DFT_FINAL   — OPT BP86/def2-TZVP (gas) + SP BP86/def2-TZVPD
#
#   CPCM (solvation):
#     3. ORCA_XTB2_ALPB   — XTB2 OPT ALPB(water) pre-screening
#     4. ORCA_DFT_CPCM_FAST  — OPT BP86/def2-TZVP(-f) + CPCM
#     5. ORCA_DFT_CPCM_FINAL — OPT BP86/def2-TZVP + CPCM, then
#                               SP  BP86/def2-TZVPD + CPCM
# ================================================================

# ----------------------------------------------------------------
# Gas-phase input writers
# ----------------------------------------------------------------

"""
    write_dft_fast_input(filepath, xyz_filename, charge, multiplicity,
                         nprocs, maxcore; optimize_kw="OPT")

Write ORCA input for gas-phase geometry optimisation at BP86/def2-TZVP(-f).
Equivalent to Python ORCA_DFT_FAST.  Output base name: `geo_opt`.
"""
function write_dft_fast_input(filepath::String, xyz_filename::String,
                              charge::Int, multiplicity::Int,
                              nprocs::Int, maxcore::Int;
                              optimize_kw::String = "OPT")
    pal_str = nprocs > 1 ? " PAL$nprocs" : ""
    open(filepath, "w") do io
        println(io, "%MaxCore $maxcore")
        println(io)
        println(io, "! DFT $optimize_kw BP86 def2-TZVP(-f)$pal_str")
        println(io)
        println(io, "%base \"geo_opt\"")
        println(io)
        println(io, "* xyzfile $charge $multiplicity $xyz_filename")
        println(io)
    end
    validate_orca_blocks(filepath)
end

"""
    write_dft_final_input(filepath, xyz_filename, charge, multiplicity,
                          nprocs, maxcore; optimize_kw="OPT", do_opt=true)

Write ORCA input for gas-phase DFT final calculation (two-job compound):
  Job 1  — BP86/def2-TZVP geometry optimisation (if `do_opt`)
  Job 2  — BP86/def2-TZVPD single point (gives E_gas)
Equivalent to Python ORCA_DFT_FINAL.
"""
function write_dft_final_input(filepath::String, xyz_filename::String,
                               charge::Int, multiplicity::Int,
                               nprocs::Int, maxcore::Int;
                               optimize_kw::String = "OPT",
                               do_opt::Bool = true)
    pal_str = nprocs > 1 ? " PAL$nprocs" : ""
    open(filepath, "w") do io
        println(io, "%MaxCore $maxcore")
        println(io)
        if do_opt
            println(io, "! DFT $optimize_kw BP86 def2-TZVP$pal_str")
            println(io)
            println(io, "%base \"geo_opt_tzvp\"")
            println(io)
            println(io, "* xyzfile $charge $multiplicity $xyz_filename")
            println(io)
            println(io, "\$new_job")
            println(io)
        end
        println(io, "! def2-TZVPD SP$pal_str")
        println(io)
        println(io, "%base \"single_point_tzvpd\"")
        println(io)
        sp_xyz = do_opt ? "geo_opt_tzvp.xyz" : xyz_filename
        println(io, "* xyzfile $charge $multiplicity $sp_xyz")
        println(io)
    end
    validate_orca_blocks(filepath)
end

# ----------------------------------------------------------------
# CPCM / solvation input writers
# ----------------------------------------------------------------

"""
    format_cpcm_block(cpcm_radii; cut_area_bohr2=nothing) -> Vector{String}

Return the lines for an ORCA `%cpcm ... end` block with custom radii.
"""
function format_cpcm_block(cpcm_radii::Dict{Int,Float64};
                           cut_area_bohr2::Union{Nothing,Float64} = nothing)
    lines = String["%cpcm"]
    for (Z, r) in sort(collect(cpcm_radii))
        push!(lines, "radius[$Z]  $r")
    end
    if cut_area_bohr2 !== nothing
        push!(lines, "cut_area $cut_area_bohr2")
    end
    push!(lines, "end")
    return lines
end

"""
    write_xtb2_alpb_input(filepath, xyz_filename, charge, multiplicity,
                          maxcore; solvent="water", optimize_kw="OPT")

Write ORCA input for XTB2 + ALPB solvation pre-screening.
Equivalent to Python ORCA_XTB2_ALPB. Output base name: `geo_opt`.
"""
function write_xtb2_alpb_input(filepath::String, xyz_filename::String,
                               charge::Int, multiplicity::Int,
                               maxcore::Int;
                               solvent::String = "water",
                               optimize_kw::String = "OPT")
    open(filepath, "w") do io
        println(io, "%MaxCore $maxcore")
        println(io)
        println(io, "! XTB2 $optimize_kw ALPB($solvent)")
        println(io)
        println(io, "%base \"geo_opt\"")
        println(io)
        println(io, "* xyzfile $charge $multiplicity $xyz_filename")
        println(io)
    end
    validate_orca_blocks(filepath)
end

"""
    write_cpcm_fast_input(filepath, xyz_filename, charge, multiplicity,
                          nprocs, maxcore, cpcm_radii; optimize_kw="OPT")

Write ORCA input for CPCM geometry optimisation at BP86/def2-TZVP(-f).
Equivalent to Python ORCA_DFT_CPCM_FAST.  Output base name: `geo_opt`.
"""
function write_cpcm_fast_input(filepath::String, xyz_filename::String,
                               charge::Int, multiplicity::Int,
                               nprocs::Int, maxcore::Int,
                               cpcm_radii::Dict{Int,Float64};
                               optimize_kw::String = "OPT")
    pal_str = nprocs > 1 ? " PAL$nprocs" : ""
    open(filepath, "w") do io
        println(io, "%MaxCore $maxcore")
        println(io)
        println(io, "%base \"geo_opt\"")
        for line in format_cpcm_block(cpcm_radii)
            println(io, line)
        end
        println(io)
        println(io, "! DFT $optimize_kw CPCM BP86 def2-TZVP(-f)$pal_str")
        println(io)
        println(io, "* xyzfile $charge $multiplicity $xyz_filename")
        println(io)
    end
    validate_orca_blocks(filepath)
end

"""
    write_cpcm_final_input(filepath, xyz_filename, charge, multiplicity,
                           nprocs, maxcore, cpcm_radii;
                           optimize_kw="OPT", do_opt=true,
                           cut_area_angstrom2=nothing,
                           calculate_polarizabilities=true)

Write ORCA input for CPCM final calculation (two-job compound):
  Job 1  — OPT BP86/def2-TZVP + CPCM  (if `do_opt`)
  Job 2  — SP  BP86/def2-TZVPD + CPCM  (with optional atomic polarizabilities)
Equivalent to Python ORCA_DFT_CPCM_FINAL.

The `%cpcm` block (custom radii, cut_area) is placed before the first `!` line
and persists across `\$new_job`.
"""
function write_cpcm_final_input(filepath::String, xyz_filename::String,
                                charge::Int, multiplicity::Int,
                                nprocs::Int, maxcore::Int,
                                cpcm_radii::Dict{Int,Float64};
                                optimize_kw::String = "OPT",
                                do_opt::Bool = true,
                                cut_area_angstrom2::Union{Nothing,Float64} = nothing,
                                calculate_polarizabilities::Bool = true)
    pal_str = nprocs > 1 ? " PAL$nprocs" : ""

    cut_area_bohr2 = nothing
    if cut_area_angstrom2 !== nothing
        cut_area_bohr2 = cut_area_angstrom2 / (ANGSTROM_PER_BOHR^2)
    end

    open(filepath, "w") do io
        println(io, "%MaxCore $maxcore")
        println(io)

        # Base name for first job (or only job when skipping optimisation)
        if do_opt
            println(io, "%base \"geo_opt_tzvp\"")
        else
            println(io, "%base \"single_point_tzvpd\"")
        end
        println(io)

        # %cpcm block — persists across $new_job
        for line in format_cpcm_block(cpcm_radii; cut_area_bohr2)
            println(io, line)
        end
        println(io)

        if do_opt
            println(io, "! $optimize_kw CPCM BP86 def2-TZVP$pal_str")
            println(io)
            println(io, "* xyzfile $charge $multiplicity $xyz_filename")
            println(io)
            println(io, "\$new_job")
            println(io)
        end

        # SP job — ORCA 6.0 bug: polarizability must run single-core
        sp_pal = calculate_polarizabilities ? "" : pal_str
        println(io, "! CPCM BP86 def2-TZVPD SP$sp_pal")
        println(io)

        if calculate_polarizabilities
            println(io, "%elprop")
            println(io, "        Polar 1")
            println(io, "        Polaratom 1")
            println(io, "end")
            println(io)
        end

        println(io, "%base \"single_point_tzvpd\"")
        println(io)

        sp_xyz = do_opt ? "geo_opt_tzvp.xyz" : xyz_filename
        println(io, "* xyzfile $charge $multiplicity $sp_xyz")
        println(io)
    end
    validate_orca_blocks(filepath)
end

# ----------------------------------------------------------------
# Enhanced .orcacosmo builder for multi-job workflows
# ----------------------------------------------------------------

"""
    build_orcacosmo_from_files(; outfile, xyzfile, method_str, energy, dipole,
                                 cpcmfile=nothing, cpcm_corr_file=nothing)

Assemble a `.orcacosmo` file from explicit file paths.
Unlike `build_orcacosmo`, this accepts pre-parsed energy/dipole and lets the
caller specify arbitrary paths for the XYZ / CPCM / CPCM-corrected files
(needed when Job 1 and Job 2 use different `%base` names).
"""
function build_orcacosmo_from_files(;
    outfile::String,
    xyzfile::String,
    method_str::String,
    energy::Float64,
    dipole::Union{Nothing, NTuple{3,Float64}} = nothing,
    cpcmfile::Union{String,Nothing} = nothing,
    cpcm_corr_file::Union{String,Nothing} = nothing)

    structname = replace(basename(outfile), ".orcacosmo" => "")

    open(outfile, "w") do io
        println(io, "$structname : $method_str")

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

        if cpcmfile !== nothing && isfile(cpcmfile)
            println(io)
            println(io, "#"^50)
            println(io, "#COSMO")
            for line in eachline(cpcmfile)
                println(io, line)
            end
        end

        if cpcm_corr_file !== nothing && isfile(cpcm_corr_file)
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

# ----------------------------------------------------------------
# CPCM radii loader
# ----------------------------------------------------------------

"""
    load_cpcm_radii(filepath) -> Dict{Int,Float64}

Load CPCM radii from a tab-separated file (`atomic_number<TAB>radius`).
"""
function load_cpcm_radii(filepath::String)
    radii = Dict{Int,Float64}()
    for line in eachline(filepath)
        parts = split(strip(line), '\t')
        if length(parts) == 2
            Z = parse(Int, parts[1])
            r = parse(Float64, parts[2])
            radii[Z] = r
        end
    end
    return radii
end

# ----------------------------------------------------------------
# Step runner with OPT → COPT retry
# ----------------------------------------------------------------

"""
    run_step_with_retry(write_fn, input_path, orca_path;
                        enable_copt_fallback=true) -> OrcaResult

Run an ORCA job with automatic COPT fallback on geometry-optimisation failure.

`write_fn(optimize_kw::String)` is a callable that (re)writes the input file
using the given optimisation keyword ("OPT" or "COPT").
"""
function run_step_with_retry(write_fn::Function, input_path::String,
                             orca_path::String;
                             enable_copt_fallback::Bool = true)
    # First attempt — standard OPT
    write_fn("OPT")
    run_orca_job(input_path, orca_path)
    output_file = replace(input_path, ".inp" => ".out")

    try
        return check_orca_success(output_file)
    catch e
        if e isa OrcaConvergenceError && enable_copt_fallback
            @warn "Convergence failed, retrying with COPT…" logfile=e.logfile
            write_fn("COPT")
            run_orca_job(input_path, orca_path)
            return check_orca_success(output_file)
        else
            rethrow()
        end
    end
end

# ----------------------------------------------------------------
# Main workflow
# ----------------------------------------------------------------

"""
    run_cosmo_workflow(conformers, state; kwargs...)
        -> (gas_energy, conformers)

Run the full COSMO-RS quantum-chemistry workflow for a set of conformers.

## Gas-phase pipeline (produces E_gas reference energy)
1. **DFT_FAST** — geometry optimisation of every conformer at
   BP86/def2-TZVP(-f) in the gas phase.
2. Sort by energy, keep only the lowest-energy conformer.
3. **DFT_FINAL** — re-optimise at BP86/def2-TZVP (gas), then
   single-point at BP86/def2-TZVPD.  The SP energy is `E_gas`.

## CPCM pipeline (produces `.orcacosmo` sigma-surface files)
4. **XTB2_ALPB** — fast solvation pre-screening of the *original*
   conformers with XTB2/ALPB(water).  Filter to ≤ 3.
5. **DFT_CPCM_FAST** — CPCM optimisation at BP86/def2-TZVP(-f).
   Filter to ≤ 1.
6. **DFT_CPCM_FINAL** — CPCM optimisation at BP86/def2-TZVP, then
   single-point at BP86/def2-TZVPD + CPCM.  Produces `.orcacosmo`.

## Returns
`NamedTuple` with fields:
- `gas_energy` — gas-phase energy at BP86/def2-TZVPD (Hartree)
- `conformers` — `Vector{Conformer}` with CPCM energies and `.orcacosmo` paths

## Keyword arguments
- `base_dir`       — root directory for calculation files (default `"orca_jobs"`)
- `nprocs`         — number of parallel ORCA cores (default `4`)
- `maxcore`        — max memory per core in MB (default `2000`)
- `cpcm_radii`     — `Dict{Int,Float64}` mapping atomic number → CPCM radius
- `energy_window_kJmol` — energy window for conformer filtering (default `25.104` ≈ 6 kcal/mol)
- `rms_threshold`  — RMSD threshold for deduplication in Å (default `1.0`)
- `orca_path`      — path to ORCA executable (default: auto-detect)
- `do_geometry_optimization` — run geometry optimisations (default `true`; set `false` for single atoms)
- `max_cpcm_conformers_after_xtb`  — max conformers after XTB pre-screen (default `3`)
- `max_cpcm_conformers_after_fast` — max conformers after CPCM fast DFT (default `1`)
- `calculate_polarizabilities` — include atomic polarizabilities in CPCM final SP (default `true`)
"""
function run_cosmo_workflow(
    conformers::Vector{Conformer},
    state::MolecularState;
    base_dir::String = "orca_jobs",
    nprocs::Int = 4,
    maxcore::Int = 2000,
    cpcm_radii::Dict{Int,Float64} = copy(DEFAULT_CPCM_RADII),
    energy_window_kJmol::Float64 = 25.104,
    rms_threshold::Float64 = 1.0,
    orca_path::Union{String,Nothing} = nothing,
    do_geometry_optimization::Bool = true,
    max_cpcm_conformers_after_xtb::Int = 3,
    max_cpcm_conformers_after_fast::Int = 1,
    calculate_polarizabilities::Bool = true)

    # Resolve ORCA executable path
    if orca_path === nothing
        resolved = Sys.which("orca")
        resolved === nothing && error(
            "ORCA executable not found in PATH. " *
            "Pass orca_path explicitly to run_cosmo_workflow.")
        orca_path = resolved
    end

    # Validate CPCM radii coverage
    if !isempty(cpcm_radii)
        all_Z = Set{Int}()
        for conf in conformers
            union!(all_Z, conf.atomic_numbers)
        end
        missing_Z = setdiff(all_Z, keys(cpcm_radii))
        if !isempty(missing_Z)
            error("No CPCM radii specified for atomic numbers: $(sort(collect(missing_Z))). " *
                  "Provide them in the cpcm_radii dict or use load_cpcm_radii().")
        end
    end

    mol_dir = joinpath(base_dir, state.name)
    mkpath(mol_dir)

    # Keep original conformers for the CPCM pipeline
    original_conformers = copy(conformers)

    # ============================================================
    # GAS-PHASE PIPELINE
    # ============================================================

    gas_energy_hartree = NaN

    if do_geometry_optimization
        # --- Step 1: DFT_FAST — BP86/def2-TZVP(-f) gas-phase OPT ---
        @info "[$(state.name)] Gas DFT_FAST — $(length(conformers)) conformers"

        gas_fast_dir = joinpath(mol_dir, "01_gas_dft_fast")
        gas_fast_conformers = Conformer[]

        for (i, conf) in enumerate(conformers)
            job_dir = joinpath(gas_fast_dir, "conf_$i")
            mkpath(job_dir)

            xyz_fn = "input.xyz"
            write_xyz_file(joinpath(job_dir, xyz_fn), conf.coordinates,
                           conf.atomic_numbers, "$(state.name) conformer $i")

            inp = joinpath(job_dir, "input.inp")
            result = run_step_with_retry(inp, orca_path) do opt_kw
                write_dft_fast_input(inp, xyz_fn, state.charge, state.multiplicity,
                                     nprocs, maxcore; optimize_kw=opt_kw)
            end

            opt_coords, _ = parse_xyz_file(joinpath(job_dir, "geo_opt.xyz"))
            push!(gas_fast_conformers,
                  Conformer(opt_coords, conf.atomic_numbers,
                            result.energy * KJMOL_PER_HARTREE, nothing))
        end

        # Sort by energy, keep only the lowest
        sort!(gas_fast_conformers, by = c -> c.energy)
        best = gas_fast_conformers[1]
        @info "[$(state.name)] Gas DFT_FAST best energy = $(round(best.energy; digits=2)) kJ/mol"

        # --- Step 2: DFT_FINAL — BP86/def2-TZVP OPT + BP86/def2-TZVPD SP ---
        @info "[$(state.name)] Gas DFT_FINAL — optimise + single point"

        job_dir = joinpath(mol_dir, "02_gas_dft_final", "conf_best")
        mkpath(job_dir)

        xyz_fn = "input.xyz"
        write_xyz_file(joinpath(job_dir, xyz_fn), best.coordinates,
                       best.atomic_numbers, "$(state.name) best conformer")

        inp = joinpath(job_dir, "input.inp")
        result = run_step_with_retry(inp, orca_path) do opt_kw
            write_dft_final_input(inp, xyz_fn, state.charge, state.multiplicity,
                                  nprocs, maxcore; optimize_kw=opt_kw, do_opt=true)
        end

        gas_energy_hartree = result.energy
    else
        # No optimisation — single-point only
        @info "[$(state.name)] Gas DFT_FINAL — single point (no optimisation)"

        job_dir = joinpath(mol_dir, "02_gas_dft_final", "conf_1")
        mkpath(job_dir)

        conf = conformers[1]
        xyz_fn = "input.xyz"
        write_xyz_file(joinpath(job_dir, xyz_fn), conf.coordinates,
                       conf.atomic_numbers, state.name)

        inp = joinpath(job_dir, "input.inp")
        result = run_step_with_retry(inp, orca_path) do opt_kw
            write_dft_final_input(inp, xyz_fn, state.charge, state.multiplicity,
                                  nprocs, maxcore; optimize_kw=opt_kw, do_opt=false)
        end

        gas_energy_hartree = result.energy
    end

    @info "[$(state.name)] E_gas (def2-TZVPD) = $gas_energy_hartree Hartree"

    # ============================================================
    # CPCM PIPELINE
    # ============================================================

    cpcm_conformers = original_conformers

    if do_geometry_optimization
        # --- Step 3: XTB2_ALPB — fast solvation pre-screening ---
        @info "[$(state.name)] CPCM XTB2_ALPB — $(length(cpcm_conformers)) conformers"

        xtb_dir = joinpath(mol_dir, "03_cpcm_xtb2_alpb")
        xtb_results = Conformer[]

        for (i, conf) in enumerate(cpcm_conformers)
            job_dir = joinpath(xtb_dir, "conf_$i")
            mkpath(job_dir)

            xyz_fn = "input.xyz"
            write_xyz_file(joinpath(job_dir, xyz_fn), conf.coordinates,
                           conf.atomic_numbers, "$(state.name) conformer $i")

            inp = joinpath(job_dir, "input.inp")
            result = run_step_with_retry(inp, orca_path) do opt_kw
                write_xtb2_alpb_input(inp, xyz_fn, state.charge, state.multiplicity,
                                      maxcore; optimize_kw=opt_kw)
            end

            opt_coords, _ = parse_xyz_file(joinpath(job_dir, "geo_opt.xyz"))
            push!(xtb_results,
                  Conformer(opt_coords, conf.atomic_numbers,
                            result.energy * KJMOL_PER_HARTREE, nothing))
        end

        # Filter: energy window → RMS → cap
        sort!(xtb_results, by = c -> c.energy)
        if length(xtb_results) > 1
            e_min = xtb_results[1].energy
            xtb_results = filter(c -> (c.energy - e_min) <= energy_window_kJmol, xtb_results)
        end
        if rms_threshold > 0 && length(xtb_results) > 1
            xtb_results = prune_by_rms(xtb_results; threshold=rms_threshold)
        end
        if length(xtb_results) > max_cpcm_conformers_after_xtb
            xtb_results = xtb_results[1:max_cpcm_conformers_after_xtb]
        end
        @info "[$(state.name)] After XTB2_ALPB filtering — $(length(xtb_results)) conformers"

        cpcm_conformers = xtb_results

        # --- Step 4: DFT_CPCM_FAST — BP86/def2-TZVP(-f) + CPCM ---
        @info "[$(state.name)] CPCM DFT_CPCM_FAST — $(length(cpcm_conformers)) conformers"

        cpcm_fast_dir = joinpath(mol_dir, "04_cpcm_dft_fast")
        cpcm_fast_results = Conformer[]

        for (i, conf) in enumerate(cpcm_conformers)
            job_dir = joinpath(cpcm_fast_dir, "conf_$i")
            mkpath(job_dir)

            xyz_fn = "input.xyz"
            write_xyz_file(joinpath(job_dir, xyz_fn), conf.coordinates,
                           conf.atomic_numbers, "$(state.name) conformer $i")

            inp = joinpath(job_dir, "input.inp")
            result = run_step_with_retry(inp, orca_path) do opt_kw
                write_cpcm_fast_input(inp, xyz_fn, state.charge, state.multiplicity,
                                      nprocs, maxcore, cpcm_radii; optimize_kw=opt_kw)
            end

            opt_coords, _ = parse_xyz_file(joinpath(job_dir, "geo_opt.xyz"))
            push!(cpcm_fast_results,
                  Conformer(opt_coords, conf.atomic_numbers,
                            result.energy * KJMOL_PER_HARTREE, nothing))
        end

        # Filter: energy window → RMS → cap
        sort!(cpcm_fast_results, by = c -> c.energy)
        if length(cpcm_fast_results) > 1
            e_min = cpcm_fast_results[1].energy
            cpcm_fast_results = filter(c -> (c.energy - e_min) <= energy_window_kJmol,
                                       cpcm_fast_results)
        end
        if rms_threshold > 0 && length(cpcm_fast_results) > 1
            cpcm_fast_results = prune_by_rms(cpcm_fast_results; threshold=rms_threshold)
        end
        if length(cpcm_fast_results) > max_cpcm_conformers_after_fast
            cpcm_fast_results = cpcm_fast_results[1:max_cpcm_conformers_after_fast]
        end
        @info "[$(state.name)] After DFT_CPCM_FAST filtering — $(length(cpcm_fast_results)) conformers"

        cpcm_conformers = cpcm_fast_results
    end

    # --- Step 5: DFT_CPCM_FINAL — BP86/def2-TZVP + CPCM OPT, then
    #                                BP86/def2-TZVPD + CPCM SP ---
    @info "[$(state.name)] CPCM DFT_CPCM_FINAL — $(length(cpcm_conformers)) conformers"

    cpcm_final_dir = joinpath(mol_dir, "05_cpcm_dft_final")
    final_conformers = Conformer[]

    for (i, conf) in enumerate(cpcm_conformers)
        job_dir = joinpath(cpcm_final_dir, "conf_$i")
        mkpath(job_dir)

        xyz_fn = "input.xyz"
        write_xyz_file(joinpath(job_dir, xyz_fn), conf.coordinates,
                       conf.atomic_numbers, "$(state.name) conformer $i")

        inp = joinpath(job_dir, "input.inp")
        result = run_step_with_retry(inp, orca_path) do opt_kw
            write_cpcm_final_input(inp, xyz_fn, state.charge, state.multiplicity,
                                   nprocs, maxcore, cpcm_radii;
                                   optimize_kw=opt_kw,
                                   do_opt=do_geometry_optimization,
                                   cut_area_angstrom2=0.01,
                                   calculate_polarizabilities=calculate_polarizabilities)
        end

        # Coordinates: from OPT job if optimised, else from input
        if do_geometry_optimization
            opt_coords, _ = parse_xyz_file(joinpath(job_dir, "geo_opt_tzvp.xyz"))
        else
            opt_coords = conf.coordinates
        end

        # Assemble .orcacosmo
        structname = "$(state.name)_conf_$i"
        orcacosmo_path = joinpath(job_dir, "$structname.orcacosmo")

        build_orcacosmo_from_files(
            outfile        = orcacosmo_path,
            xyzfile        = do_geometry_optimization ?
                             joinpath(job_dir, "geo_opt_tzvp.xyz") :
                             joinpath(job_dir, xyz_fn),
            cpcmfile       = joinpath(job_dir, "single_point_tzvpd.cpcm"),
            cpcm_corr_file = joinpath(job_dir, "single_point_tzvpd.cpcm_corr"),
            method_str     = "DFT_CPCM_BP86_def2-TZVP+def2-TZVPD_SP",
            energy         = result.energy,
            dipole         = result.dipole)

        push!(final_conformers,
              Conformer(opt_coords, conf.atomic_numbers,
                        result.energy * KJMOL_PER_HARTREE, orcacosmo_path))
    end

    @info "[$(state.name)] COSMO workflow complete — " *
          "E_gas=$gas_energy_hartree Ha, $(length(final_conformers)) CPCM conformer(s)"

    return (gas_energy = gas_energy_hartree, conformers = final_conformers)
end