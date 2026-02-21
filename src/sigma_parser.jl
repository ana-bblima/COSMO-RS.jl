# Add later all the other terms


function parse_sigma_surface(filepath::String)

    seg_area = Float64[]
    seg_sigma = Float64[]
    seg_pos = Float64[]
    energy = 0.0

    lines = readlines(filepath)

    for line in lines
        if occursin("FINAL SINGLE POINT ENERGY", line)
            energy = parse(Float64, split(line)[end])
        end
    end

    # Expand this to parse:
    # - segment areas
    # - charges
    # - coordinates

    return SigmaSurface(seg_area, seg_sigma, seg_pos, energy)
end