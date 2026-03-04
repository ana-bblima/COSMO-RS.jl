using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test

# Load the COSMORS module
include(joinpath(@__DIR__, "..", "src", "COSMORS.jl"))
using .COSMORS

# Import non-exported internals needed for testing
import .COSMORS: rmsd, select_atoms, prune_by_rms,
                 validate_orca_blocks, build_orcacosmo,
                 R_gas, filter_energy_window,
                 SigmaSurface, SigmaProfile, AminoAcidPKa,
                 ORCA_CPCM_SOLVENTS,
                 gaussian_sigma_averaging!,
                 ANGSTROM_PER_BOHR, KJMOL_PER_HARTREE,
                 polyprotic_fractions

# Path to the existing ORCA test output (methane single-point)
const ORCA_TEST_DIR = joinpath(@__DIR__, "..", "test_orca_run", "methane_conf_1")
const ORCA_TEST_OUTPUT = joinpath(ORCA_TEST_DIR, "input.out")


# ══════════════════════════════════════════════════════════════
#  1. STATE GENERATION
# ══════════════════════════════════════════════════════════════

@testset "State Generation" begin

    @testset "All 20 amino acids are defined" begin
        @test length(AMINO_ACID_NAMES) == 20
        expected = sort([
            "glycine", "alanine", "valine", "leucine", "isoleucine",
            "proline", "phenylalanine", "tryptophan", "methionine",
            "serine", "threonine", "cysteine", "tyrosine",
            "asparagine", "glutamine", "aspartic_acid", "glutamic_acid",
            "lysine", "arginine", "histidine"
        ])
        @test AMINO_ACID_NAMES == expected
    end

    @testset "Every amino acid returns non-empty states" begin
        for name in AMINO_ACID_NAMES
            states = generate_states(name)
            @test length(states) >= 4  # at least neutral, zwitterion, protonated, deprotonated
            @test all(s -> !isempty(s.smiles), states)
            @test all(s -> !isempty(s.name), states)
        end
    end

    @testset "Standard states present for simple amino acids" begin
        simple = ["glycine", "alanine", "valine", "leucine", "isoleucine",
                  "proline", "phenylalanine", "tryptophan", "methionine",
                  "serine", "threonine", "asparagine", "glutamine"]
        for name in simple
            states = generate_states(name)
            state_names = [s.name for s in states]
            @test "neutral" in state_names || any(n -> startswith(n, "neutral"), state_names)
            @test "zwitterion" in state_names
            @test "protonated" in state_names
            @test "deprotonated" in state_names
        end
    end

    @testset "Ionizable side chain amino acids have extra states" begin
        # Cysteine: thiolate
        cys = generate_states("cysteine")
        @test any(s -> s.name == "thiolate", cys)

        # Tyrosine: phenolate
        tyr = generate_states("tyrosine")
        @test any(s -> s.name == "phenolate", tyr)

        # Aspartic acid: doubly_deprotonated
        asp = generate_states("aspartic_acid")
        @test any(s -> s.name == "doubly_deprotonated", asp)

        # Glutamic acid: doubly_deprotonated
        glu = generate_states("glutamic_acid")
        @test any(s -> s.name == "doubly_deprotonated", glu)

        # Lysine: doubly_protonated
        lys = generate_states("lysine")
        @test any(s -> s.name == "doubly_protonated", lys)

        # Arginine: guanidinium
        arg = generate_states("arginine")
        @test any(s -> s.name == "guanidinium", arg)

        # Histidine: imidazolium + two tautomers
        his = generate_states("histidine")
        @test any(s -> s.name == "imidazolium", his)
        @test any(s -> s.name == "neutral_epsilon", his)
        @test any(s -> s.name == "neutral_delta", his)
    end

    @testset "Charge and multiplicity consistency" begin
        for name in AMINO_ACID_NAMES
            states = generate_states(name)
            for s in states
                # Multiplicity must be positive
                @test s.multiplicity >= 1
                # All amino acid states defined here are closed-shell singlets
                @test s.multiplicity == 1

                # Check charge labels are consistent with names
                if s.name == "neutral" || startswith(s.name, "neutral_")
                    @test s.charge == 0
                end
                if s.name == "zwitterion"
                    @test s.charge == 0
                end
            end
        end
    end

    @testset "Lookup by three-letter code" begin
        @test generate_states("gly") == generate_states("glycine")
        @test generate_states("Ala") == generate_states("alanine")
        @test generate_states("TRP") == generate_states("tryptophan")
        @test generate_states("asp") == generate_states("aspartic_acid")
        @test generate_states("his") == generate_states("histidine")
    end

    @testset "Lookup by one-letter code" begin
        @test generate_states("G") == generate_states("glycine")
        @test generate_states("a") == generate_states("alanine")
        @test generate_states("W") == generate_states("tryptophan")
        @test generate_states("D") == generate_states("aspartic_acid")
        @test generate_states("H") == generate_states("histidine")
    end

    @testset "Unknown amino acid throws error" begin
        @test_throws ErrorException generate_states("unknown_molecule")
        @test_throws ErrorException generate_states("xyz")
        @test_throws ErrorException generate_states("")
    end

    @testset "Glycine backward compatibility" begin
        states = generate_states("glycine")
        @test states[1].name == "neutral"
        @test states[1].smiles == "NCC(=O)O"
        @test states[1].charge == 0
        @test states[1].multiplicity == 1
        @test states[2].name == "zwitterion"
        @test states[2].smiles == "[NH3+]CC(=O)[O-]"
        @test states[2].charge == 0
    end
end


# ══════════════════════════════════════════════════════════════
#  2. CONFORMER GENERATION
# ══════════════════════════════════════════════════════════════

@testset "Conformer Generation" begin

    @testset "Basic conformer generation from glycine neutral" begin
        state = MolecularState("neutral", "NCC(=O)O", 0, 1)
        confs = generate_conformers(state; nconfs=5, rms_threshold=0.0)
        @test length(confs) >= 1
        @test all(c -> size(c.coordinates, 2) == 3, confs)
        @test all(c -> length(c.atomic_numbers) == size(c.coordinates, 1), confs)
        @test all(c -> isfinite(c.energy), confs)
        @test all(c -> c.sigma_surface_file === nothing, confs)
    end

    @testset "Conformer generation from glycine zwitterion" begin
        state = MolecularState("zwitterion", "[NH3+]CC(=O)[O-]", 0, 1)
        confs = generate_conformers(state; nconfs=5, rms_threshold=0.0)
        @test length(confs) >= 1
        @test all(c -> size(c.coordinates, 2) == 3, confs)
    end

    @testset "Conformer generation from alanine" begin
        state = generate_states("alanine")[1]  # neutral
        confs = generate_conformers(state; nconfs=5)
        @test length(confs) >= 1
        @test all(c -> length(c.atomic_numbers) > 0, confs)
    end

    @testset "Atomic numbers are sensible" begin
        state = MolecularState("neutral", "NCC(=O)O", 0, 1)
        confs = generate_conformers(state; nconfs=3, rms_threshold=0.0)
        @test length(confs) >= 1
        # Glycine: H, C, N, O atoms only (Z = 1, 6, 7, 8)
        for c in confs
            @test all(z -> z in [1, 6, 7, 8], c.atomic_numbers)
        end
    end

    @testset "Energy window filtering" begin
        state = MolecularState("test", "CCCCCCCC", 0, 1)  # octane, many conformers
        confs = generate_conformers(state; nconfs=30, rms_threshold=0.0)
        n_before = length(confs)

        confs_filtered = generate_conformers(state; nconfs=30,
                                              rms_threshold=0.0,
                                              energy_window_kJmol=5.0)
        @test length(confs_filtered) <= n_before

        if length(confs_filtered) > 1
            energies = [c.energy for c in confs_filtered]
            Emin = minimum(energies)
            @test all(e -> (e - Emin) <= 5.0 * 1000.0 + 1e-6, energies)
        end
    end

    @testset "RMS pruning reduces conformers" begin
        state = MolecularState("test", "CCCCCCCC", 0, 1)

        # Generate without pruning, then prune the same set
        confs_no_prune = generate_conformers(state; nconfs=20, rms_threshold=0.0)
        confs_pruned = prune_by_rms(confs_no_prune; threshold=1.0)

        @test length(confs_pruned) <= length(confs_no_prune)
        @test length(confs_pruned) >= 1
    end

    @testset "RMSD function" begin
        P = [1.0 0.0 0.0;
             0.0 1.0 0.0;
             0.0 0.0 1.0]
        @test rmsd(P, P) < 1e-10

        Q = P .+ 5.0  # translation only
        @test rmsd(P, Q) < 1e-10

        R = [2.0 0.0 0.0;
             0.0 2.0 0.0;
             0.0 0.0 2.0]
        @test rmsd(P, R) > 0.0
    end

    @testset "select_atoms heavy only" begin
        coords = [0.0 0.0 0.0;
                  1.0 0.0 0.0;
                  2.0 0.0 0.0]
        atomic_numbers = [6, 1, 8]  # C, H, O
        conf = Conformer(coords, atomic_numbers, 0.0, nothing)

        heavy = select_atoms(conf; heavy_only=true)
        @test size(heavy) == (2, 3)  # C and O only

        all_atoms = select_atoms(conf; heavy_only=false)
        @test size(all_atoms) == (3, 3)
    end
end


# ══════════════════════════════════════════════════════════════
#  3. ORCA INPUT VALIDATION
# ══════════════════════════════════════════════════════════════

@testset "ORCA Input Validation" begin

    # Construct settings bypassing Sys.which by using the 6-arg inner constructor
    test_settings = OrcaSettings("BP86", "def2-SVP", "Water", 1, 1000, "/dev/null")

    # Water molecule: H2O
    water_coords = [0.0 0.0 0.0;
                    0.757 0.586 0.0;
                   -0.757 0.586 0.0]
    water_Z = [8, 1, 1]

    @testset "Valid input passes" begin
        @test validate_orca_input(water_coords, water_Z, 0, 1, test_settings) === nothing
    end

    @testset "Dimension mismatch" begin
        bad_coords = [0.0 0.0 0.0;
                      1.0 0.0 0.0]  # 2 rows but 3 atoms
        @test_throws OrcaInputError validate_orca_input(bad_coords, water_Z, 0, 1, test_settings)
    end

    @testset "NaN coordinates" begin
        nan_coords = copy(water_coords)
        nan_coords[1, 1] = NaN
        @test_throws OrcaInputError validate_orca_input(nan_coords, water_Z, 0, 1, test_settings)
    end

    @testset "Inf coordinates" begin
        inf_coords = copy(water_coords)
        inf_coords[2, 3] = Inf
        @test_throws OrcaInputError validate_orca_input(inf_coords, water_Z, 0, 1, test_settings)
    end

    @testset "Overlapping atoms" begin
        overlap_coords = [0.0 0.0 0.0;
                          0.05 0.0 0.0;  # too close (< 0.1 A)
                          2.0 0.0 0.0]
        @test_throws OrcaInputError validate_orca_input(overlap_coords, water_Z, 0, 1, test_settings)
    end

    @testset "Charge/multiplicity inconsistency" begin
        # Water: 10 electrons (charge 0) -> even -> multiplicity must be odd (1,3,5...)
        # Multiplicity 2 (even) is inconsistent with even electron count
        @test_throws OrcaInputError validate_orca_input(water_coords, water_Z, 0, 2, test_settings)

        # Charge +1: 9 electrons -> odd -> multiplicity must be even (2,4,6...)
        # Multiplicity 1 is inconsistent
        @test_throws OrcaInputError validate_orca_input(water_coords, water_Z, 1, 1, test_settings)
    end

    @testset "Negative electron count" begin
        # Absurd charge that yields non-physical electron count
        @test_throws OrcaInputError validate_orca_input(water_coords, water_Z, 100, 1, test_settings)
    end

    @testset "Settings sanity: nprocs < 1" begin
        bad_settings = OrcaSettings("BP86", "def2-SVP", "Water", 0, 1000, "/dev/null")
        @test_throws OrcaInputError validate_orca_input(water_coords, water_Z, 0, 1, bad_settings)
    end

    @testset "Settings sanity: maxcore < 100" begin
        bad_settings = OrcaSettings("BP86", "def2-SVP", "Water", 1, 50, "/dev/null")
        @test_throws OrcaInputError validate_orca_input(water_coords, water_Z, 0, 1, bad_settings)
    end

    @testset "Unknown solvent produces warning (not error)" begin
        warn_settings = OrcaSettings("BP86", "def2-SVP", "kryptonite", 1, 1000, "/dev/null")
        @test_logs (:warn, r"not in known ORCA CPCM solvent") validate_orca_input(
            water_coords, water_Z, 0, 1, warn_settings)
    end
end


# ══════════════════════════════════════════════════════════════
#  4. ORCA INPUT FILE WRITING
# ══════════════════════════════════════════════════════════════

@testset "ORCA Input File Writing" begin

    test_settings = OrcaSettings("BP86", "def2-SVP", "Water", 2, 2000, "/dev/null")

    water_coords = [0.0 0.0 0.0;
                    0.757 0.586 0.0;
                   -0.757 0.586 0.0]
    water_Z = [8, 1, 1]

    @testset "write_orca_input creates valid file" begin
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "test_water.inp")

        write_orca_input(filepath, water_coords, water_Z, 0, 1, test_settings)
        @test isfile(filepath)

        content = read(filepath, String)

        # Check key components
        @test occursin("%maxcore 2000", content)
        @test occursin("%pal nprocs 2 end", content)
        @test occursin("BP86", content)
        @test occursin("def2-SVP", content)
        @test occursin("CPCM(Water)", content)
        @test occursin("TightSCF", content)
        @test occursin("* xyz 0 1", content)
        @test occursin("*", content)

        rm(tmpdir; recursive=true)
    end

    @testset "write_orca_input with optimize keyword" begin
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "test_opt.inp")

        write_orca_input(filepath, water_coords, water_Z, 0, 1, test_settings;
                         optimize_keyword="OPT")
        content = read(filepath, String)

        @test occursin("OPT", content)

        rm(tmpdir; recursive=true)
    end

    @testset "write_orca_input with COPT keyword" begin
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "test_copt.inp")

        write_orca_input(filepath, water_coords, water_Z, 0, 1, test_settings;
                         optimize_keyword="COPT")
        content = read(filepath, String)

        @test occursin("COPT", content)

        rm(tmpdir; recursive=true)
    end

    @testset "Single-proc input omits %pal block" begin
        single_settings = OrcaSettings("BP86", "def2-SVP", "Water", 1, 1000, "/dev/null")
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "test_single.inp")

        write_orca_input(filepath, water_coords, water_Z, 0, 1, single_settings)
        content = read(filepath, String)

        @test !occursin("%pal", content)

        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  5. ORCA BLOCK VALIDATION
# ══════════════════════════════════════════════════════════════

@testset "ORCA Block Validation" begin

    @testset "Balanced blocks pass" begin
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "balanced.inp")
        write(filepath, """
        %maxcore 2000
        %pal nprocs 4 end
        %cpcm
        end
        ! BP86 def2-SVP
        * xyz 0 1
        8 0.0 0.0 0.0
        *
        """)
        @test validate_orca_blocks(filepath) === nothing
        rm(tmpdir; recursive=true)
    end

    @testset "Unclosed block throws error" begin
        tmpdir = mktempdir()
        filepath = joinpath(tmpdir, "unclosed.inp")
        write(filepath, """
        %maxcore 2000
        %cpcm
        ! BP86 def2-SVP
        """)
        @test_throws OrcaInputError validate_orca_blocks(filepath)
        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  6. ORCA OUTPUT PARSING & CONVERGENCE CHECKING
# ══════════════════════════════════════════════════════════════

@testset "ORCA Output Parsing" begin

    # Tests using the real methane output if available
    if isfile(ORCA_TEST_OUTPUT)

        @testset "check_orca_success with real output" begin
            result = check_orca_success(ORCA_TEST_OUTPUT)
            @test result.terminated_normally == true
            @test result.converged == true
            @test isfinite(result.energy)
            @test result.energy < 0.0  # Hartree energy is negative
            @test result.energy ≈ -40.481771719006 atol=1e-6
            @test result.logfile == ORCA_TEST_OUTPUT
        end

        @testset "parse_orca_energy with real output" begin
            energy = parse_orca_energy(ORCA_TEST_OUTPUT)
            @test energy ≈ -40.481771719006 atol=1e-6
        end

        @testset "Dipole moment parsed from real output" begin
            result = check_orca_success(ORCA_TEST_OUTPUT)
            @test result.dipole !== nothing
            dx, dy, dz = result.dipole
            @test isfinite(dx) && isfinite(dy) && isfinite(dz)
        end

    else
        @warn "Skipping real ORCA output tests: $ORCA_TEST_OUTPUT not found"
    end

    @testset "check_orca_success: missing TERMINATED NORMALLY" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "no_termination.out")
        write(logfile, """
        Some ORCA output
        FINAL SINGLE POINT ENERGY     -76.123456789
        """)
        @test_throws OrcaTerminationError check_orca_success(logfile)
        rm(tmpdir; recursive=true)
    end

    @testset "check_orca_success: SCF not converged" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "scf_fail.out")
        write(logfile, """
        SCF NOT CONVERGED AFTER 128 CYCLES
        FINAL SINGLE POINT ENERGY     -76.123456789
        ****ORCA TERMINATED NORMALLY****
        """)
        @test_throws OrcaConvergenceError check_orca_success(logfile)
        rm(tmpdir; recursive=true)
    end

    @testset "check_orca_success: geometry optimization not converged" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "opt_fail.out")
        write(logfile, """
        The optimization did not converge but reached the maximum number of steps
        FINAL SINGLE POINT ENERGY     -76.123456789
        ****ORCA TERMINATED NORMALLY****
        """)
        @test_throws OrcaConvergenceError check_orca_success(logfile)
        rm(tmpdir; recursive=true)
    end

    @testset "check_orca_success: terminated but no energy" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "no_energy.out")
        write(logfile, """
        ****ORCA TERMINATED NORMALLY****
        """)
        @test_throws OrcaTerminationError check_orca_success(logfile)
        rm(tmpdir; recursive=true)
    end

    @testset "check_orca_success: successful minimal output" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "success.out")
        write(logfile, """
        FINAL SINGLE POINT ENERGY     -76.123456789000
        x,y,z [Debye]:    1.234    -0.567    0.890
        ****ORCA TERMINATED NORMALLY****
        """)
        result = check_orca_success(logfile)
        @test result.terminated_normally == true
        @test result.converged == true
        @test result.energy ≈ -76.123456789 atol=1e-9
        @test result.dipole !== nothing
        @test result.dipole[1] ≈ 1.234 atol=1e-3
        @test result.dipole[2] ≈ -0.567 atol=1e-3
        @test result.dipole[3] ≈ 0.890 atol=1e-3
        rm(tmpdir; recursive=true)
    end

    @testset "parse_orca_energy: no energy found" begin
        tmpdir = mktempdir()
        logfile = joinpath(tmpdir, "empty.out")
        write(logfile, "no energy here\n")
        @test_throws ErrorException parse_orca_energy(logfile)
        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  7. BUILD ORCACOSMO FILE
# ══════════════════════════════════════════════════════════════

@testset "Build orcacosmo" begin

    @testset "orcacosmo file construction" begin
        tmpdir = mktempdir()

        # Create a minimal xyz file
        xyzfile = joinpath(tmpdir, "mol.xyz")
        write(xyzfile, "3\ntest molecule\nO   0.0  0.0  0.0\nH   0.757  0.586  0.0\nH  -0.757  0.586  0.0\n")

        outfile = build_orcacosmo(tmpdir, "mol", "BP86 def2-SVP CPCM(Water)",
                                   -76.123, (1.0, -0.5, 0.3))

        @test isfile(outfile)
        content = read(outfile, String)

        @test occursin("#ENERGY", content)
        @test occursin("-76.123", content)
        @test occursin("#DIPOLE", content)
        @test occursin("#XYZ_FILE", content)
        @test occursin("O   0.0", content)

        rm(tmpdir; recursive=true)
    end

    @testset "orcacosmo without dipole" begin
        tmpdir = mktempdir()

        xyzfile = joinpath(tmpdir, "mol.xyz")
        write(xyzfile, "1\ntest\nC   0.0  0.0  0.0\n")

        outfile = build_orcacosmo(tmpdir, "mol", "method", -40.0, nothing)
        content = read(outfile, String)

        @test occursin("#ENERGY", content)
        @test !occursin("#DIPOLE", content)

        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  8. BOLTZMANN WEIGHTS & THERMODYNAMICS
# ══════════════════════════════════════════════════════════════

@testset "Boltzmann Weights" begin

    @testset "Weights sum to 1" begin
        energies = [0.0, 1.0, 5.0]  # kJ/mol
        w = boltzmann_weights(energies)
        @test sum(w) ≈ 1.0
    end

    @testset "Single energy gives weight 1" begin
        w = boltzmann_weights([42.0])
        @test w ≈ [1.0]
    end

    @testset "Degenerate energies give equal weights" begin
        w = boltzmann_weights([0.0, 0.0, 0.0])
        @test all(w .≈ 1/3)
    end

    @testset "Lowest energy gets highest weight" begin
        w = boltzmann_weights([0.0, 10.0, 50.0])  # kJ/mol
        @test w[1] > w[2] > w[3]
    end

    @testset "Large energy gap gives nearly-zero weight for high-energy state" begin
        w = boltzmann_weights([0.0, 100.0])  # 100 kJ/mol gap
        @test w[1] ≈ 1.0 atol=1e-6
        @test w[2] < 1e-6
    end

    @testset "Temperature dependence" begin
        energies = [0.0, 5.0]  # kJ/mol

        w_low = boltzmann_weights(energies; T=100.0)
        w_high = boltzmann_weights(energies; T=1000.0)

        # At higher T, populations are more equal
        @test abs(w_high[1] - w_high[2]) < abs(w_low[1] - w_low[2])
    end
end

@testset "State Free Energy" begin

    @testset "Single conformer returns its energy" begin
        conf = Conformer(zeros(1, 3), [6], 1000.0, nothing)
        G = state_free_energy([conf])
        @test G ≈ 1000.0 atol=1e-6
    end

    @testset "Degenerate conformers lower free energy" begin
        conf1 = Conformer(zeros(1, 3), [6], 0.0, nothing)
        conf2 = Conformer(ones(1, 3), [6], 0.0, nothing)

        G_single = state_free_energy([conf1])
        G_double = state_free_energy([conf1, conf2])

        # Adding a degenerate state lowers free energy by -RT*ln(2)
        @test G_double < G_single
        expected_diff = -R_gas * 298.15 * log(2)
        @test (G_double - G_single) ≈ expected_diff atol=1.0
    end
end

@testset "State Population Weights" begin

    @testset "Single state gets weight 1" begin
        s = MolecularState("test", "C", 0, 1)
        c = Conformer(zeros(1, 3), [6], 0.0, nothing)
        ensemble = StateEnsemble(s, [c])

        w = state_population_weights([ensemble])
        @test w ≈ [1.0]
    end

    @testset "Equal-energy states get equal weights" begin
        s1 = MolecularState("s1", "C", 0, 1)
        s2 = MolecularState("s2", "C", 0, 1)
        c1 = Conformer(zeros(1, 3), [6], 0.0, nothing)
        c2 = Conformer(zeros(1, 3), [6], 0.0, nothing)

        e1 = StateEnsemble(s1, [c1])
        e2 = StateEnsemble(s2, [c2])

        w = state_population_weights([e1, e2])
        @test w[1] ≈ 0.5 atol=1e-10
        @test w[2] ≈ 0.5 atol=1e-10
    end

    @testset "Lower-energy state dominates" begin
        s1 = MolecularState("low", "C", 0, 1)
        s2 = MolecularState("high", "C", 0, 1)
        c_low = Conformer(zeros(1, 3), [6], 0.0, nothing)
        c_high = Conformer(zeros(1, 3), [6], 100.0, nothing)  # 100 kJ/mol higher

        e1 = StateEnsemble(s1, [c_low])
        e2 = StateEnsemble(s2, [c_high])

        w = state_population_weights([e1, e2])
        @test w[1] > w[2]
        @test w[1] ≈ 1.0 atol=1e-6
    end
end


# ══════════════════════════════════════════════════════════════
#  9. PARSE SIGMA SURFACE
# ══════════════════════════════════════════════════════════════

const ORCA_TEST_ORCACOSMO = joinpath(ORCA_TEST_DIR, "input.orcacosmo")

@testset "Parse sigma surface" begin

    if isfile(ORCA_TEST_ORCACOSMO)

        @testset "Basic parsing of methane orcacosmo" begin
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)

            # Should have 510 surface segments
            @test length(surface.seg_area) == 510
            @test length(surface.seg_sigma_raw) == 510
            @test size(surface.seg_pos) == (510, 3)

            # Energy should be finite and negative (converted to kJ/mol)
            @test isfinite(surface.energy)
            @test surface.energy < 0.0

            # Verify energy conversion: -40.481771719006 Hartree * 2625.499639479
            expected_energy = -40.481771719006 * KJMOL_PER_HARTREE
            @test surface.energy ≈ expected_energy atol=1.0  # within 1 kJ/mol

            # All areas should be non-negative
            @test all(a -> a >= 0.0, surface.seg_area)

            # Positions should be in Angstrom range (not bohr)
            # methane is ~1A radius, surface points should be < ~5A from origin
            max_dist = maximum(sqrt.(sum(surface.seg_pos.^2, dims=2)))
            @test max_dist < 5.0  # Angstrom
            @test max_dist > 0.5  # Should be nonzero

            # seg_sigma_avg should be nothing (not yet averaged)
            @test surface.seg_sigma_avg === nothing
        end

        @testset "First segment spot-check values" begin
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)

            # First surface point from orcacosmo file:
            # X=3.866161596 Y=0.018609683 Z=-0.020114566 (bohr)
            # area=0.139721426 (bohr^2)
            expected_x1 = 3.866161596 * ANGSTROM_PER_BOHR
            expected_y1 = 0.018609683 * ANGSTROM_PER_BOHR
            expected_z1 = -0.020114566 * ANGSTROM_PER_BOHR
            expected_area1 = 0.139721426 * ANGSTROM_PER_BOHR^2

            @test surface.seg_pos[1, 1] ≈ expected_x1 atol=1e-6
            @test surface.seg_pos[1, 2] ≈ expected_y1 atol=1e-6
            @test surface.seg_pos[1, 3] ≈ expected_z1 atol=1e-6
            @test surface.seg_area[1] ≈ expected_area1 atol=1e-10
        end

        @testset "Sigma averaging works on parsed surface" begin
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)

            gaussian_sigma_averaging!(surface)

            @test surface.seg_sigma_avg !== nothing
            @test length(surface.seg_sigma_avg) == 510
            @test all(isfinite, surface.seg_sigma_avg)
        end

        @testset "Full pipeline: parse -> sigma profile" begin
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)
            profile = compute_sigma_profile(surface)

            @test profile isa SigmaProfile
            @test length(profile.sigma_grid) > 0
            @test length(profile.area_distribution) == length(profile.sigma_grid)
            @test sum(profile.area_distribution) > 0.0
        end

    else
        @warn "Skipping sigma surface parsing tests: $ORCA_TEST_ORCACOSMO not found"
    end
end


# ══════════════════════════════════════════════════════════════
#  9b. SIGMA PROFILE (linear interpolation binning)
# ══════════════════════════════════════════════════════════════

@testset "Sigma profile linear interpolation" begin

    @testset "Area is conserved after binning" begin
        # Create a small synthetic surface with known values
        n = 5
        seg_area = [0.1, 0.2, 0.15, 0.05, 0.3]
        seg_sigma = [0.0, 0.005, -0.01, 0.0125, -0.005]
        seg_pos = zeros(n, 3)
        surface = SigmaSurface(seg_area, seg_sigma, seg_pos, -100.0)

        profile = compute_sigma_profile(surface; use_averaged=false)

        # Total area must be conserved
        @test sum(profile.area_distribution) ≈ sum(seg_area) atol=1e-10
    end

    @testset "Exact grid-point sigma goes to single bin" begin
        # sigma = 0.005 should land exactly on a grid point (step=0.001)
        seg_area = [1.0]
        seg_sigma = [0.005]
        seg_pos = zeros(1, 3)
        surface = SigmaSurface(seg_area, seg_sigma, seg_pos, 0.0)

        grid = -0.03:0.001:0.03
        profile = compute_sigma_profile(surface; sigma_range=grid, use_averaged=false)

        # Find the bin at sigma=0.005
        idx = findfirst(x -> x ≈ 0.005, profile.sigma_grid)
        @test profile.area_distribution[idx] ≈ 1.0
        @test sum(profile.area_distribution) ≈ 1.0
    end

    @testset "Mid-point sigma splits 50/50" begin
        # sigma = 0.0005 is exactly between 0.0 and 0.001
        seg_area = [2.0]
        seg_sigma = [0.0005]
        seg_pos = zeros(1, 3)
        surface = SigmaSurface(seg_area, seg_sigma, seg_pos, 0.0)

        grid = -0.03:0.001:0.03
        profile = compute_sigma_profile(surface; sigma_range=grid, use_averaged=false)

        idx_left = findfirst(x -> x ≈ 0.0, profile.sigma_grid)
        idx_right = findfirst(x -> x ≈ 0.001, profile.sigma_grid)

        @test profile.area_distribution[idx_left] ≈ 1.0 atol=1e-10
        @test profile.area_distribution[idx_right] ≈ 1.0 atol=1e-10
        @test sum(profile.area_distribution) ≈ 2.0
    end

    @testset "Quarter-point sigma splits 75/25" begin
        # sigma = 0.00025 is 25% between 0.0 and 0.001
        seg_area = [4.0]
        seg_sigma = [0.00025]
        seg_pos = zeros(1, 3)
        surface = SigmaSurface(seg_area, seg_sigma, seg_pos, 0.0)

        grid = -0.03:0.001:0.03
        profile = compute_sigma_profile(surface; sigma_range=grid, use_averaged=false)

        idx_left = findfirst(x -> x ≈ 0.0, profile.sigma_grid)
        idx_right = findfirst(x -> x ≈ 0.001, profile.sigma_grid)

        # 75% goes to left, 25% to right
        @test profile.area_distribution[idx_left] ≈ 3.0 atol=1e-10
        @test profile.area_distribution[idx_right] ≈ 1.0 atol=1e-10
    end

    @testset "Real data: area conservation from orcacosmo" begin
        if isfile(ORCA_TEST_ORCACOSMO)
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)
            profile = compute_sigma_profile(surface)

            total_seg_area = sum(surface.seg_area)
            total_profile_area = sum(profile.area_distribution)

            @test total_profile_area ≈ total_seg_area atol=1e-6
        end
    end
end


# ══════════════════════════════════════════════════════════════
#  9c. SIGMA MOMENTS
# ══════════════════════════════════════════════════════════════

@testset "Sigma moments" begin

    @testset "Basic sigma moments from real data" begin
        if isfile(ORCA_TEST_ORCACOSMO)
            surface = parse_sigma_surface(ORCA_TEST_ORCACOSMO)
            moments = calculate_sigma_moments(surface)

            @test length(moments.sigma_moments) == 7
            @test length(moments.hb_acceptor_moments) == 7
            @test length(moments.hb_donor_moments) == 7

            # 0th moment = total surface area
            @test moments.sigma_moments[1] ≈ sum(surface.seg_area) atol=1e-6

            # 1st moment = total charge (sum of sigma * area)
            gaussian_sigma_averaging!(surface)
            expected_m1 = sum(surface.seg_sigma_avg .* surface.seg_area)
            @test moments.sigma_moments[2] ≈ expected_m1 atol=1e-10

            # All moments should be finite
            @test all(isfinite, moments.sigma_moments)
            @test all(isfinite, moments.hb_acceptor_moments)
            @test all(isfinite, moments.hb_donor_moments)

            # HB moments at indices 1 and 2 should be zero (only computed for i=2,3,4)
            @test moments.hb_acceptor_moments[1] == 0.0
            @test moments.hb_acceptor_moments[2] == 0.0
            @test moments.hb_donor_moments[1] == 0.0
            @test moments.hb_donor_moments[2] == 0.0
        end
    end

    @testset "Sigma moments with uniform charge" begin
        # Uniform sigma=0 -> all moments except 0th should be ~0
        n = 10
        seg_area = fill(0.1, n)
        seg_sigma = zeros(n)
        seg_pos = randn(n, 3)
        surface = SigmaSurface(seg_area, seg_sigma, seg_pos, 0.0)

        moments = calculate_sigma_moments(surface; use_averaged=false)

        @test moments.sigma_moments[1] ≈ 1.0  # total area = 10 * 0.1
        for i in 2:7
            @test moments.sigma_moments[i] ≈ 0.0 atol=1e-15
        end
    end
end
# ══════════════════════════════════════════════════════════════
#  10. BUILD ORCACOSMO WITH COSMO_CORRECTED
# ══════════════════════════════════════════════════════════════

@testset "Build orcacosmo with COSMO_corrected" begin

    @testset "cpcm_corr file is included" begin
        tmpdir = mktempdir()

        # Create minimal xyz file
        xyzfile = joinpath(tmpdir, "mol.xyz")
        write(xyzfile, "1\ntest\nC   0.0  0.0  0.0\n")

        # Create minimal cpcm file
        cpcmfile = joinpath(tmpdir, "mol.cpcm")
        write(cpcmfile, "# fake cpcm data\n1  # Number of atoms\n")

        # Create cpcm_corr file
        corrfile = joinpath(tmpdir, "mol.cpcm_corr")
        write(corrfile, """Corrected dielectric energy   =     -0.001234000
Total C-PCM charge            =     -0.000000000
C-PCM corrected charges:
    -0.000050000
     0.000050000
""")

        outfile = build_orcacosmo(tmpdir, "mol", "BP86 def2-SVP CPCM(Water)",
                                   -76.123, (1.0, -0.5, 0.3))

        content = read(outfile, String)

        @test occursin("#COSMO_corrected", content)
        @test occursin("Corrected dielectric energy", content)
        @test occursin("-0.001234000", content)
        @test occursin("C-PCM corrected charges:", content)

        rm(tmpdir; recursive=true)
    end

    @testset "no cpcm_corr file -> no corrected section" begin
        tmpdir = mktempdir()

        xyzfile = joinpath(tmpdir, "mol.xyz")
        write(xyzfile, "1\ntest\nC   0.0  0.0  0.0\n")

        outfile = build_orcacosmo(tmpdir, "mol", "method", -40.0, nothing)
        content = read(outfile, String)

        @test !occursin("#COSMO_corrected", content)

        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  11. ROUND-TRIP: BUILD THEN PARSE
# ══════════════════════════════════════════════════════════════

@testset "Round-trip: build then parse" begin

    if isdir(ORCA_TEST_DIR)

        @testset "Round-trip with real ORCA output" begin
            tmpdir = mktempdir()

            # Copy necessary files to tmpdir
            for fname in ["input.xyz", "input.cpcm", "input.cpcm_corr"]
                src = joinpath(ORCA_TEST_DIR, fname)
                if isfile(src)
                    cp(src, joinpath(tmpdir, fname))
                end
            end

            # Build orcacosmo from the test data
            energy_hartree = -40.481771719006
            dipole = (-0.005633, -0.010018, -0.012657)

            outfile = build_orcacosmo(tmpdir, "input",
                                       "BP86 def2-SVP CPCM(Water)",
                                       energy_hartree, dipole)

            @test isfile(outfile)

            # Verify the corrected section is present
            content = read(outfile, String)
            @test occursin("#COSMO_corrected", content)

            # Parse it back
            surface = parse_sigma_surface(outfile)

            @test length(surface.seg_area) == 510
            @test isfinite(surface.energy)
            @test size(surface.seg_pos) == (510, 3)

            # Energy should be adjusted with corrected dielectric
            @test surface.energy < 0.0

            rm(tmpdir; recursive=true)
        end

    else
        @warn "Skipping round-trip test: $ORCA_TEST_DIR not found"
    end

    @testset "Parse without corrected section (fallback)" begin
        tmpdir = mktempdir()

        xyzfile = joinpath(tmpdir, "mol.xyz")
        write(xyzfile, "1\ntest\nC   0.0  0.0  0.0\n")

        # Build without cpcm or cpcm_corr -> no COSMO sections
        outfile = build_orcacosmo(tmpdir, "mol", "method", -40.0, nothing)

        # Parser should still work (returns empty segments)
        surface = parse_sigma_surface(outfile)
        @test length(surface.seg_area) == 0
        @test surface.energy ≈ -40.0 * KJMOL_PER_HARTREE atol=1.0

        rm(tmpdir; recursive=true)
    end
end


# ══════════════════════════════════════════════════════════════
#  12. pKa DATA & pH-DEPENDENT POPULATION WEIGHTING
# ══════════════════════════════════════════════════════════════

@testset "pKa Data" begin

    @testset "All amino acids have pKa data" begin
        for name in AMINO_ACID_NAMES
            @test haskey(AMINO_ACID_PKA, name)
            pKa = AMINO_ACID_PKA[name]
            @test length(pKa.pKa_values) >= 2
            @test issorted(pKa.pKa_values)
        end
    end

    @testset "pKa values are physically reasonable" begin
        for (name, pKa) in AMINO_ACID_PKA
            for val in pKa.pKa_values
                @test 0.0 < val < 15.0
            end
            @test pKa.charge_fully_deprotonated ∈ [-2, -1]
        end
    end
end

@testset "Polyprotic Fractions" begin

    @testset "Monoprotic: pH = pKa gives 50/50" begin
        f = polyprotic_fractions([4.0], 4.0)
        @test length(f) == 2
        @test f[1] ≈ 0.5 atol=1e-10  # protonated
        @test f[2] ≈ 0.5 atol=1e-10  # deprotonated
    end

    @testset "Monoprotic: low pH gives mostly protonated" begin
        f = polyprotic_fractions([4.0], 1.0)  # 3 units below pKa
        @test f[1] > 0.999
        @test f[2] < 0.001
    end

    @testset "Monoprotic: high pH gives mostly deprotonated" begin
        f = polyprotic_fractions([4.0], 7.0)  # 3 units above pKa
        @test f[1] < 0.001
        @test f[2] > 0.999
    end

    @testset "Fractions sum to 1" begin
        for pKa_vals in [[2.34, 9.60], [1.88, 3.65, 9.60], [1.82, 6.00, 9.17]]
            for pH in [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 14.0]
                f = polyprotic_fractions(pKa_vals, pH)
                @test sum(f) ≈ 1.0 atol=1e-10
                @test all(f .>= 0.0)
            end
        end
    end

    @testset "Diprotic glycine at characteristic pH values" begin
        pKa = [2.34, 9.60]

        # At pH 1 (well below pKa1): fully protonated dominates
        f = polyprotic_fractions(pKa, 1.0)
        @test f[1] > 0.95   # NH3+-COOH

        # At pH 6 (between pKa1 and pKa2): singly deprotonated dominates
        f = polyprotic_fractions(pKa, 6.0)
        @test f[2] > 0.99   # zwitterion + neutral level

        # At pH 12 (well above pKa2): fully deprotonated dominates
        f = polyprotic_fractions(pKa, 12.0)
        @test f[3] > 0.99   # NH2-COO-
    end

    @testset "Triprotic histidine at pH 7" begin
        pKa = [1.82, 6.00, 9.17]
        f = polyprotic_fractions(pKa, 7.0)

        # At pH 7: 1 proton removed from fully protonated is dominant
        # (between pKa2=6.0 and pKa3=9.17)
        @test f[3] > 0.5    # 2 protons removed (zwitterion level)
        @test f[2] > f[1]   # 1 removed > 0 removed
        @test f[3] > f[4]   # 2 removed > 3 removed
    end
end


@testset "pH-Dependent State Population Weights" begin

    @testset "Glycine at pH 6: zwitterion + neutral level dominates" begin
        states = generate_states("glycine")
        # Create fake ensembles with energy differences
        # Zwitterion much lower energy than neutral (CPCM stabilization)
        ensembles = StateEnsemble[]
        for s in states
            # Zwitterion 20 kJ/mol lower than neutral; others at 0
            E = 0.0
            if s.name == "zwitterion"
                E = -20.0  # more stable in solution
            end
            c = Conformer(zeros(1, 3), [6], E, nothing)
            push!(ensembles, StateEnsemble(s, [c]))
        end

        w = state_population_weights_pH(ensembles, "glycine", 6.0)

        @test length(w) == 4
        @test sum(w) ≈ 1.0 atol=1e-10

        # Find indices
        idx_zwit = findfirst(e -> e.state.name == "zwitterion", ensembles)
        idx_neut = findfirst(e -> e.state.name == "neutral", ensembles)
        idx_prot = findfirst(e -> e.state.name == "protonated", ensembles)
        idx_dprt = findfirst(e -> e.state.name == "deprotonated", ensembles)

        # At pH 6: singly deprotonated level (zwitterion + neutral) should dominate
        @test w[idx_zwit] + w[idx_neut] > 0.99

        # Zwitterion should dominate over neutral due to lower energy
        @test w[idx_zwit] > w[idx_neut]
    end

    @testset "Glycine at pH 1: protonated dominates" begin
        states = generate_states("glycine")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w = state_population_weights_pH(ensembles, "glycine", 1.0)

        idx_prot = findfirst(e -> e.state.name == "protonated", ensembles)
        @test w[idx_prot] > 0.9
    end

    @testset "Glycine at pH 12: deprotonated dominates" begin
        states = generate_states("glycine")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w = state_population_weights_pH(ensembles, "glycine", 12.0)

        idx_dprt = findfirst(e -> e.state.name == "deprotonated", ensembles)
        @test w[idx_dprt] > 0.99
    end

    @testset "Weights sum to 1 for all amino acids at pH 7" begin
        for name in AMINO_ACID_NAMES
            states = generate_states(name)
            ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                         for s in states]

            w = state_population_weights_pH(ensembles, name, 7.0)
            @test sum(w) ≈ 1.0 atol=1e-10
            @test all(w .>= 0.0)
        end
    end

    @testset "Three-letter and one-letter code lookup" begin
        states = generate_states("glycine")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w_full = state_population_weights_pH(ensembles, "glycine", 7.0)
        w_three = state_population_weights_pH(ensembles, "gly", 7.0)
        w_one = state_population_weights_pH(ensembles, "g", 7.0)

        @test w_full ≈ w_three
        @test w_full ≈ w_one
    end

    @testset "Direct pKa/charge interface matches amino acid lookup" begin
        states = generate_states("glycine")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w_lookup = state_population_weights_pH(ensembles, "glycine", 7.0)
        w_direct = state_population_weights_pH(
            ensembles, [2.34, 9.60], 7.0;
            charge_fully_deprotonated=-1)

        @test w_lookup ≈ w_direct
    end

    @testset "Histidine at pH 7: imidazolium vs protonated tautomers" begin
        states = generate_states("histidine")
        # At pH 7, level with 1 proton removed from fully protonated
        # (protonated + imidazolium, both charge +1) should have some weight
        # but the 2-protons-removed level (zwitterion + neutrals) dominates
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w = state_population_weights_pH(ensembles, "histidine", 7.0)

        # Charge 0 states (zwitterion + neutrals) should dominate at pH 7
        charge0_weight = sum(w[i] for i in eachindex(ensembles)
                            if ensembles[i].state.charge == 0)
        @test charge0_weight > 0.5
    end

    @testset "Aspartic acid at pH 7: doubly_deprotonated level" begin
        states = generate_states("aspartic_acid")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w = state_population_weights_pH(ensembles, "aspartic_acid", 7.0)

        # At pH 7 (between pKa2=3.65 and pKa3=9.60):
        # 2-protons-removed level should dominate
        # That's deprotonated + doubly_deprotonated (both charge -1)
        charge_m1_weight = sum(w[i] for i in eachindex(ensembles)
                              if ensembles[i].state.charge == -1)
        @test charge_m1_weight > 0.9
    end

    @testset "Energy-only weights differ from pH-dependent weights" begin
        states = generate_states("glycine")
        ensembles = [StateEnsemble(s, [Conformer(zeros(1,3), [6], 0.0, nothing)])
                     for s in states]

        w_energy = state_population_weights(ensembles)
        w_pH = state_population_weights_pH(ensembles, "glycine", 7.0)

        # With equal energies, energy-only gives equal weights
        @test all(w_energy .≈ 0.25)

        # pH-dependent weights are NOT equal (pH selects protonation level)
        @test !all(w_pH .≈ 0.25)
    end
end


# ══════════════════════════════════════════════════════════════
#  13. SMILES VALIDITY (via MoleculeFlow)
# ══════════════════════════════════════════════════════════════

@testset "SMILES Validity" begin
    using MoleculeFlow

    @testset "All amino acid SMILES are valid in MoleculeFlow" begin
        for name in AMINO_ACID_NAMES
            states = generate_states(name)
            for s in states
                mol = mol_from_smiles(s.smiles)
                @test mol.valid
            end
        end
    end
end
