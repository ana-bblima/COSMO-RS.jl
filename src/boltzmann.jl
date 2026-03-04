const R_gas = 8.314e-3  # kJ/(mol·K)

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


############
# pH-DEPENDENT POPULATION WEIGHTING
############

"""
    polyprotic_fractions(pKa_values, pH)

Compute macroscopic protonation-level fractions for a polyprotic acid.

Given `n` sorted (ascending) pKa values, returns a vector of `n+1` fractions:
- `fractions[1]`   = fully protonated (0 protons removed)
- `fractions[j+1]` = `j` protons removed
- `fractions[n+1]` = fully deprotonated

Uses the standard polyprotic acid dissociation formalism:

    H_n A  ⇌  H_{n-1}A⁻ + H⁺      Ka₁
    H_{n-1}A⁻  ⇌  H_{n-2}A²⁻ + H⁺  Ka₂
    ...

The fraction at level j (j protons removed) is:

    α_j = (Ka₁·Ka₂·...·Ka_j / [H⁺]^j) / Z

where Z = Σ_j α_j (unnormalized).
"""
function polyprotic_fractions(pKa_values::Vector{Float64}, pH::Float64)
    n = length(pKa_values)
    H = 10.0^(-pH)

    # terms[j+1] = Π_{k=1}^{j} Ka_k / H^j
    terms = zeros(n + 1)
    terms[1] = 1.0  # fully protonated reference

    cumulative_Ka = 1.0
    for j in 1:n
        cumulative_Ka *= 10.0^(-pKa_values[j])
        terms[j + 1] = cumulative_Ka / H^j
    end

    Z = sum(terms)
    return terms ./ Z
end


"""
    state_population_weights_pH(ensembles, pKa_values, pH;
                                charge_fully_deprotonated, T=298.15)

Compute pH-dependent state population weights.

Combines polyprotic acid equilibria (for macroscopic protonation levels)
with Boltzmann weighting of tautomers/conformers within each level.

Each state is assigned to a protonation level based on:

    n_protons = state.charge - charge_fully_deprotonated
    protons_removed = n_ionizable - n_protons

States at the same protonation level are tautomers; their relative
populations are determined by their conformer-averaged free energies.

Returns normalized weights (one per ensemble).
"""
function state_population_weights_pH(
        ensembles::Vector{StateEnsemble},
        pKa_values::Vector{Float64},
        pH::Float64;
        charge_fully_deprotonated::Int,
        T::Float64 = 298.15)

    n_ionizable = length(pKa_values)

    # Macrostate fractions from pKa
    fractions = polyprotic_fractions(pKa_values, pH)

    n_states = length(ensembles)
    weights = zeros(n_states)

    # Group states by protonation level and assign weights
    for j_removed in 0:n_ionizable

        # Find ensembles at this protonation level
        level_indices = Int[]
        for (idx, e) in enumerate(ensembles)
            n_protons = e.state.charge - charge_fully_deprotonated
            if n_ionizable - n_protons == j_removed
                push!(level_indices, idx)
            end
        end

        isempty(level_indices) && continue

        # Boltzmann weight tautomers within this level by their free energies
        G_tautomers = [state_free_energy(ensembles[idx].conformers; T=T)
                       for idx in level_indices]
        tautomer_weights = boltzmann_weights(G_tautomers; T=T)

        # Distribute the macrostate fraction among tautomers
        for (k, idx) in enumerate(level_indices)
            weights[idx] = fractions[j_removed + 1] * tautomer_weights[k]
        end
    end

    # Renormalize (accounts for protonation levels with no states)
    total = sum(weights)
    if total > 0
        weights ./= total
    end

    return weights
end


"""
    state_population_weights_pH(ensembles, amino_acid_name, pH; T=298.15)

Convenience method using built-in amino acid pKa data from `AMINO_ACID_PKA`.

Accepts full names ("glycine"), three-letter codes ("gly"),
or one-letter codes ("G"). Lookup is case-insensitive.
"""
function state_population_weights_pH(
        ensembles::Vector{StateEnsemble},
        amino_acid_name::String,
        pH::Float64;
        T::Float64 = 298.15)

    key = lowercase(strip(amino_acid_name))

    if haskey(AMINO_ACID_THREE_TO_FULL, key)
        key = AMINO_ACID_THREE_TO_FULL[key]
    elseif haskey(AMINO_ACID_ONE_TO_FULL, key)
        key = AMINO_ACID_ONE_TO_FULL[key]
    end

    if !haskey(AMINO_ACID_PKA, key)
        error("pKa data not found for '$amino_acid_name'. " *
              "Available: $(join(sort(collect(keys(AMINO_ACID_PKA))), ", "))")
    end

    pKa_data = AMINO_ACID_PKA[key]
    return state_population_weights_pH(
        ensembles, pKa_data.pKa_values, pH;
        charge_fully_deprotonated=pKa_data.charge_fully_deprotonated, T=T)
end