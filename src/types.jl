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
    coordinates::Matrix{Float64}
    energy::Float64
    sigma_surface_file::String
end

struct SigmaSurface
    seg_area::Vector{Float64}
    seg_sigma_raw::Vector{Float64}
    seg_pos::Matrix{Float64}
    energy::Float64
end

struct SigmaProfile
    sigma_grid::Vector{Float64}
    area_distribution::Vector{Float64}
end

struct MoleculeSystem
    name::String
    states::Vector{MolecularState}
end