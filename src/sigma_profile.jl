function compute_sigma_profile(surface::SigmaSurface;
                               sigma_range=-0.03:0.001:0.03)

    areas = zeros(length(sigma_range))

    for (σ, A) in zip(surface.seg_sigma_raw,
                      surface.seg_area)

        idx = findmin(abs.(sigma_range .- σ))[2]
        areas[idx] += A
    end

    return SigmaProfile(collect(sigma_range), areas)
end