############
# TYPES
############

struct MolecularState
    name::String
    smiles::String
    charge::Int
    multiplicity::Int
end

struct Conformer
    coordinates::Matrix{Float64}   # (n_atoms × 3)
    atomic_numbers::Vector{Int}
    energy::Float64
    sigma_surface_file::Union{Nothing,String}
end

struct SigmaSurface
    seg_area::Vector{Float64}
    seg_sigma_raw::Vector{Float64}
    seg_pos::Matrix{Float64}
    energy::Float64
    seg_sigma_avg::Union{Nothing, Vector{Float64}}
end

SigmaSurface(area, sigma, pos, energy) = SigmaSurface(area, sigma, pos, energy, nothing)

struct SigmaProfile
    sigma_grid::Vector{Float64}
    area_distribution::Vector{Float64}
end

struct MoleculeSystem
    name::String
    states::Vector{MolecularState}
end

struct StateEnsemble
    state::MolecularState
    conformers::Vector{Conformer}
end