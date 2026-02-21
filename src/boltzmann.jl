const R_gas = 8.314  # J/mol/K

function boltzmann_weights(energies;
                           T=298.15)

    β = 1/(R_gas*T)

    shifted = energies .- minimum(energies)
    w = exp.(-β .* shifted)

    return w ./ sum(w)
end


function average_sigma_profiles(profiles::Vector{SigmaProfile},
                                energies)

    w = boltzmann_weights(energies)

    σgrid = profiles[1].sigma_grid
    avg = zeros(length(σgrid))

    for i in eachindex(profiles)
        avg .+= w[i] .* profiles[i].area_distribution
    end

    return SigmaProfile(σgrid, avg)
end


"""
    state_free_energy(conformers; T=298.15)

Compute free energy of a molecular state
from its conformer ensemble.
"""
function state_free_energy(conformers::Vector{Conformer};
                           T::Float64 = 298.15)

    energies = [c.energy for c in conformers]

    Emin = minimum(energies)
    β = 1 / (R_gas * T)

    Z = sum(exp.(-β .* (energies .- Emin)))

    return Emin - (1/β) * log(Z)
end

"""
    state_population_weights(ensembles; T=298.15)

Return normalized state populations.
"""
function state_population_weights(
        ensembles::Vector{StateEnsemble};
        T::Float64 = 298.15)

    G = [state_free_energy(e.conformers; T=T)
         for e in ensembles]

    Gmin = minimum(G)
    β = 1 / (R_gas * T)

    weights = exp.(-β .* (G .- Gmin))

    return weights ./ sum(weights)
end

"""
    average_over_states(ensembles; T=298.15)

Return total sigma profile weighted
over states and conformers.
"""
function average_over_states(
        ensembles::Vector{StateEnsemble};
        T::Float64 = 298.15)

    # Compute state weights
    W = state_population_weights(ensembles; T=T)

    total_profile = nothing

    for (idx, ensemble) in enumerate(ensembles)

        # Conformer energies
        energies = [c.energy for c in ensemble.conformers]

        # Conformer weights
        w_conf = boltzmann_weights(energies; T=T)

        # Build state-averaged profile
        profiles = [
            compute_sigma_profile(
                parse_sigma_surface(c.sigma_surface_file)
            )
            for c in ensemble.conformers
        ]

        state_profile = average_sigma_profiles(
            profiles, energies)

        if total_profile === nothing
            total_profile = SigmaProfile(
                state_profile.sigma_grid,
                zeros(length(state_profile.sigma_grid))
            )
        end

        total_profile.area_distribution .+=
            W[idx] .* state_profile.area_distribution
    end

    return total_profile
end