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

mutable struct SigmaSurface
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

struct AminoAcidPKa
    pKa_values::Vector{Float64}          # sorted ascending, one per ionizable group
    charge_fully_deprotonated::Int       # charge of the fully deprotonated form
end

struct StateEnsemble
    state::MolecularState
    conformers::Vector{Conformer}
end

############
# ORCA ERROR TYPES
############

struct OrcaInputError <: Exception
    msg::String
end
Base.showerror(io::IO, e::OrcaInputError) = print(io, "OrcaInputError: ", e.msg)

struct OrcaTerminationError <: Exception
    msg::String
    logfile::String
end
Base.showerror(io::IO, e::OrcaTerminationError) =
    print(io, "OrcaTerminationError: ", e.msg, " (see ", e.logfile, ")")

struct OrcaConvergenceError <: Exception
    msg::String
    logfile::String
end
Base.showerror(io::IO, e::OrcaConvergenceError) =
    print(io, "OrcaConvergenceError: ", e.msg, " (see ", e.logfile, ")")

############
# ORCA RESULT
############

struct OrcaResult
    energy::Float64                              # Hartree
    dipole::Union{Nothing, NTuple{3, Float64}}
    terminated_normally::Bool
    converged::Bool
    logfile::String
end