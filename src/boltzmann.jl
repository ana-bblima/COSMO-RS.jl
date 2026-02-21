function boltzmann_weights(energies;
                           T=298.15)

    R = 8.314
    β = 1/(R*T)

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