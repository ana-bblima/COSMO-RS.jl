using LinearAlgebra

"""
    gaussian_sigma_averaging!(surface; r_av=0.5)

Compute Gaussian-averaged sigmas in-place.
r_av in Angstrom.
"""
function gaussian_sigma_averaging!(surface::SigmaSurface;
                                   r_av::Float64 = 0.5)

    A = surface.seg_area
    σ = surface.seg_sigma_raw
    pos = surface.seg_pos

    n = length(σ)

    @assert size(pos,1) == n
    @assert size(pos,2) == 3

    σ_avg = zeros(Float64, n)

    # Precompute squared radii from segment area
    r2 = A ./ π
    r_av2 = r_av^2

    for i in 1:n

        numerator = 0.0
        denominator = 0.0

        xi = pos[i,1]
        yi = pos[i,2]
        zi = pos[i,3]

        for j in 1:n

            dx = pos[j,1] - xi
            dy = pos[j,2] - yi
            dz = pos[j,3] - zi

            d2 = dx*dx + dy*dy + dz*dz

            denom_r = r2[j] + r_av2

            weight =
                (r2[j] * r_av2 / denom_r) *
                exp(-d2 / denom_r)

            numerator += σ[j] * weight
            denominator += weight
        end

        σ_avg[i] = numerator / denominator
    end

    surface.seg_sigma_avg = σ_avg

    return surface
end

function gaussian_sigma_averaging_fast!(surface::SigmaSurface;
                                        r_av::Float64 = 0.5)

    A = surface.seg_area
    σ = surface.seg_sigma_raw
    pos = surface.seg_pos

    n = length(σ)
    σ_avg = zeros(Float64, n)

    r2 = A ./ π
    r_av2 = r_av^2

    for i in 1:n

        d2 = sum((pos .- pos[i,:]').^2, dims=2)[:]

        denom_r = r2 .+ r_av2

        weights =
            (r2 .* r_av2 ./ denom_r) .*
            exp.(-d2 ./ denom_r)

        σ_avg[i] =
            sum(σ .* weights) /
            sum(weights)
    end

    surface.seg_sigma_avg = σ_avg

    return surface
end


function compute_sigma_profile(surface::SigmaSurface;
                               sigma_range=-0.03:0.001:0.03,
                               use_averaged=true)

    if use_averaged
        surface.seg_sigma_avg === nothing &&
            gaussian_sigma_averaging_fast!(surface)

        σ_values = surface.seg_sigma_avg
    else
        σ_values = surface.seg_sigma_raw
    end

    areas = zeros(length(sigma_range))

    for (σ, A) in zip(σ_values, surface.seg_area)
        idx = findmin(abs.(sigma_range .- σ))[2]
        areas[idx] += A
    end

    return SigmaProfile(collect(sigma_range), areas)
end