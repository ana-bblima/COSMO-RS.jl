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


"""
    compute_sigma_profile(surface; sigma_range=-0.03:0.001:0.03, use_averaged=true)

Build a discretised sigma profile from a `SigmaSurface`.

Each segment's area is distributed to the two neighbouring grid points
by linear interpolation (matching the Python `cluster_segments_into_segmenttypes`
implementation). A segment whose sigma falls exactly on a grid point
contributes its full area to that bin.

If `use_averaged` is true (default), Gaussian sigma averaging is applied
first (if not already done).
"""
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

    grid = collect(sigma_range)
    ngrid = length(grid)
    areas = zeros(ngrid)

    σ_lo = first(grid)
    σ_hi = last(grid)
    dσ = step(sigma_range)

    for (σ, A) in zip(σ_values, surface.seg_area)

        # Clamp sigma to grid range
        σ_clamped = clamp(σ, σ_lo, σ_hi)

        # Index of the left grid point (<= σ_clamped)
        idx_left = floor(Int, (σ_clamped - σ_lo) / dσ) + 1
        idx_left = clamp(idx_left, 1, ngrid)

        left_val = grid[idx_left]

        if idx_left == ngrid || left_val ≈ σ_clamped
            # Exactly on a grid point or at right boundary
            areas[idx_left] += A
        else
            # Linear interpolation between left and right grid points
            idx_right = idx_left + 1
            right_val = grid[idx_right]
            Δ = right_val - left_val

            areas[idx_left]  += A * (right_val - σ_clamped) / Δ
            areas[idx_right] += A * (σ_clamped - left_val) / Δ
        end
    end

    return SigmaProfile(grid, areas)
end


"""
    calculate_sigma_moments(surface; use_averaged=true, n_moments=7)

Compute sigma moments and hydrogen-bond donor/acceptor moments
from a `SigmaSurface`, following the Python reference implementation.

Returns a `NamedTuple` with:
- `sigma_moments`    : Vector of length `n_moments` (0th through 6th)
- `hb_acceptor_moments` : Vector of length `n_moments`
- `hb_donor_moments`    : Vector of length `n_moments`

Moments of order >= 2 are scaled by `100^order` (following the convention
used in COSMO-RS sigma-moment descriptors).
"""
function calculate_sigma_moments(surface::SigmaSurface;
                                 use_averaged::Bool = true,
                                 n_moments::Int = 7)

    if use_averaged
        surface.seg_sigma_avg === nothing &&
            gaussian_sigma_averaging_fast!(surface)
        sigmas = surface.seg_sigma_avg
    else
        sigmas = surface.seg_sigma_raw
    end

    areas = surface.seg_area

    # ── Sigma moments ──────────────────────────────────────────
    sigma_moments = zeros(n_moments)
    for i in 0:(n_moments - 1)
        sigma_moments[i + 1] = sum(sigmas .^ i .* areas)
        if i > 1
            sigma_moments[i + 1] *= 100^i
        end
    end

    # ── HB acceptor / donor moments ───────────────────────────
    hb_acceptor = zeros(n_moments)
    hb_donor    = zeros(n_moments)

    for i in [2, 3, 4]
        current_hb_threshold = 0.006 + i * 0.002

        hb_acceptor[i + 1] = sum(max.(sigmas .- current_hb_threshold, 0.0) .* areas)
        hb_donor[i + 1]    = sum(max.(-sigmas .- current_hb_threshold, 0.0) .* areas)
    end

    hb_acceptor .= 100.0 .* abs.(hb_acceptor)
    hb_donor    .= 100.0 .* abs.(hb_donor)

    return (sigma_moments       = sigma_moments,
            hb_acceptor_moments = hb_acceptor,
            hb_donor_moments    = hb_donor)
end
