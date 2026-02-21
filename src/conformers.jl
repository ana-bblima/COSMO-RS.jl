using MoleculeFlow
using LinearAlgebra


function generate_conformers(state::MolecularState;
                             nconfs=100)

    mol = MoleculeFlow.from_smiles(state.smiles)
    mol3d = MoleculeFlow.embed(mol)

    conformers = MoleculeFlow.generate_conformers(
        mol3d;
        n_conformers=nconfs,
        method=:etkdg
    )

    return conformers
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
function rmsd(P::Matrix{Float64},
              Q::Matrix{Float64})

    @assert size(P) == size(Q)

    # Center coordinates
    Pc = P .- mean(P, dims=1)
    Qc = Q .- mean(Q, dims=1)

    # Covariance matrix
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
                      threshold::Float64 = 0.5,
                      heavy_only::Bool = true)

    # Sort by energy
    sorted = sort(conformers, by=c->c.energy)

    kept = Conformer[]

    for conf in sorted

        coords_i = select_atoms(conf; heavy_only=heavy_only)

        keep_flag = true

        for kept_conf in kept

            coords_j = select_atoms(kept_conf;
                                    heavy_only=heavy_only)

            r = rmsd(coords_i, coords_j)

            if r < threshold
                keep_flag = false
                break
            end
        end

        keep_flag && push!(kept, conf)
    end

    return kept
end