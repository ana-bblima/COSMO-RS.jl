using MoleculeFlow
using LinearAlgebra
using Statistics

function generate_conformers(state::MolecularState;
                             nconfs::Int = 100,
                             optimize::Bool = true,
                             force_field::Symbol = :mmff,
                             random_seed::Union{Nothing,Int} = nothing,
                             rms_threshold::Float64 = 0.5,
                             energy_window_kJmol::Union{Nothing,Float64} = nothing)

    mol = mol_from_smiles(state.smiles)

    @assert mol.valid "Invalid SMILES for state $(state.name)"

    conformers = generate_3d_conformers(
        mol,
        nconfs;
        optimize = optimize,
        force_field = force_field,
        random_seed = random_seed
    )

    converted = Conformer[]

    for conf in conformers

        coords_full = conf.molecule.props[:coordinates_3d]

        # Get heavy atoms only (since get_atoms excludes hydrogens)
        atoms = get_atoms(conf.molecule)

        heavy_indices = Int[]
        atomic_numbers = Int[]

        for (i, atom) in enumerate(atoms)
            Z = get_atomic_number(atom)
            if Z != 1
                push!(heavy_indices, i)
                push!(atomic_numbers, Z)
            end
        end

        # Extract heavy atom coordinates
        coords = coords_full[heavy_indices, :]

        # Convert kcal/mol → J/mol
        energy_kcal = conf.conformer_result.energy
        energy = energy_kcal * 4184.0

        push!(converted,
            Conformer(
                coords,
                atomic_numbers,
                energy,
                nothing
            )
        )
    end

        # --- Energy window filter (optional) ---
    if energy_window_kJmol !== nothing
        energies = [c.energy for c in converted]
        Emin = minimum(energies)
        window_J = energy_window_kJmol * 1000.0

        converted = [
            c for c in converted
            if (c.energy - Emin) <= window_J
        ]
    end

    # --- RMS pruning ---
    if rms_threshold > 0
        converted = prune_by_rms(
            converted;
            threshold = rms_threshold
        )
    end

    return converted
end

function filter_energy_window(conformers;
                              window_kJmol=10.0)

    energies = [c.energy for c in conformers]
    Emin = minimum(energies)

    return [c for c in conformers
            if (c.energy - Emin) <= window_kJmol]
end


"""
    rmsd(P, Q)

Compute optimal RMSD after Kabsch alignment.
P and Q are (n_atoms × 3).
"""
function rmsd(P::Matrix{Float64}, Q::Matrix{Float64})
    @assert size(P) == size(Q)

    Pc = P .- mean(P, dims=1)
    Qc = Q .- mean(Q, dims=1)

    C = Pc' * Qc
    U, _, V = svd(C)

    d = sign(det(V * U'))
    D = Diagonal([1.0, 1.0, d])

    R = V * D * U'
    P_rot = Pc * R

    diff = P_rot - Qc

    return sqrt(sum(diff.^2) / size(P,1))
end

"""
    select_atoms(conf; heavy_only=true)

Return coordinate matrix filtered by atom type.
"""
function select_atoms(conf::Conformer; heavy_only=true)

    if !heavy_only
        return conf.coordinates
    end

    mask = conf.atomic_numbers .!= 1

    return conf.coordinates[mask, :]
end

"""
    prune_by_rms(conformers;
                 threshold=0.5,
                 heavy_only=true)

Remove redundant conformers using RMS threshold (Å).
Returns pruned conformer vector.
"""
function prune_by_rms(conformers::Vector{Conformer};
                      threshold::Float64 = 0.5)

    sorted = sort(conformers, by = c -> c.energy)

    kept = Conformer[]

    for conf in sorted
        keep_flag = true

        for kept_conf in kept
            r = rmsd(conf.coordinates, kept_conf.coordinates)

            if r < threshold
                keep_flag = false
                break
            end
        end

        keep_flag && push!(kept, conf)
    end

    return kept
end