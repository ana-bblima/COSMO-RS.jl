# Unit conversion constants (matching Python input_parsers.py)
const ANGSTROM_PER_BOHR = 0.52917721092       # Angstrom / bohr
const KJMOL_PER_HARTREE = 2625.499639479      # (kJ/mol) / Hartree

"""
    parse_sigma_surface(filepath::String) -> SigmaSurface

Parse an `.orcacosmo` file produced by `build_orcacosmo`.

Returns a `SigmaSurface` with:
- `seg_area`      : segment areas (Angstrom^2)
- `seg_sigma_raw` : charge density sigma = charge / area (e/Angstrom^2)
- `seg_pos`       : segment positions (Angstrom), n_segments x 3
- `energy`        : total electronic energy (kJ/mol)

If a `#COSMO_corrected` section is present, corrected charges and
dielectric energy are used. Otherwise uncorrected values are used.
"""
function parse_sigma_surface(filepath::String)

    lines = readlines(filepath)
    nlines = length(lines)

    # Energies (kJ/mol after conversion)
    energy_tot = 0.0
    energy_dielectric_uncorr = 0.0

    # Raw segment data (in atomic units, converted later)
    seg_x_bohr   = Float64[]
    seg_y_bohr   = Float64[]
    seg_z_bohr   = Float64[]
    seg_area_bohr2 = Float64[]
    seg_charge_uncorr = Float64[]

    # Corrected data (populated only if #COSMO_corrected exists)
    corrected_dielectric = nothing      # Float64 kJ/mol
    corrected_charges    = nothing      # Vector{Float64}

    i = 1
    while i <= nlines
        line = strip(lines[i])

        # ── Energy ──────────────────────────────────────────────
        if occursin("FINAL SINGLE POINT ENERGY", line)
            energy_tot = parse(Float64, split(line)[end]) * KJMOL_PER_HARTREE
        end

        # ── Uncorrected dielectric energy (inside #COSMO section) ──
        if occursin("# CPCM dielectric energy", line)
            energy_dielectric_uncorr = parse(Float64, split(line)[1]) * KJMOL_PER_HARTREE
        end

        # ── Surface points table (inside #COSMO section) ────────
        if occursin("SURFACE POINTS (A.U.)", line)
            i += 1  # skip the #--- separator line
            i += 1  # skip the column-header line (X Y Z area potential ...)
            i += 1  # first data line
            while i <= nlines
                data_line = strip(lines[i])
                if isempty(data_line) || startswith(data_line, "#")
                    break
                end
                parts = split(data_line)
                if length(parts) >= 6
                    push!(seg_x_bohr,        parse(Float64, parts[1]))
                    push!(seg_y_bohr,        parse(Float64, parts[2]))
                    push!(seg_z_bohr,        parse(Float64, parts[3]))
                    push!(seg_area_bohr2,    parse(Float64, parts[4]))
                    push!(seg_charge_uncorr, parse(Float64, parts[6]))
                end
                i += 1
            end
            continue   # i already at next section boundary
        end

        # ── COSMO corrected section ─────────────────────────────
        if occursin("#COSMO_corrected", line) || occursin("#CPCM_corrected", line)
            i += 1  # "Corrected dielectric energy = ..."
            if i <= nlines
                corr_line = strip(lines[i])
                m = match(r"Corrected\s+dielectric\s+energy\s+=\s*([0-9+\-.eE]+)", corr_line)
                if m !== nothing
                    corrected_dielectric = parse(Float64, m.captures[1]) * KJMOL_PER_HARTREE
                end
            end
            i += 1  # "Total C-PCM charge = ..."
            i += 1  # "C-PCM corrected charges:"
            i += 1  # first corrected charge value
            corrected_charges = Float64[]
            while i <= nlines
                charge_line = strip(lines[i])
                if isempty(charge_line) || startswith(charge_line, "#")
                    break
                end
                push!(corrected_charges, parse(Float64, charge_line))
                i += 1
            end
            continue
        end

        i += 1
    end

    # ── Unit conversions ────────────────────────────────────────
    n = length(seg_area_bohr2)

    seg_area = seg_area_bohr2 .* (ANGSTROM_PER_BOHR^2)

    seg_pos = Matrix{Float64}(undef, n, 3)
    seg_pos[:, 1] = seg_x_bohr .* ANGSTROM_PER_BOHR
    seg_pos[:, 2] = seg_y_bohr .* ANGSTROM_PER_BOHR
    seg_pos[:, 3] = seg_z_bohr .* ANGSTROM_PER_BOHR

    # ── Choose corrected or uncorrected charges ─────────────────
    if corrected_charges !== nothing && length(corrected_charges) == n
        seg_sigma_raw = corrected_charges ./ seg_area
        # Adjust total energy with corrected dielectric
        if corrected_dielectric !== nothing
            energy_tot = energy_tot - energy_dielectric_uncorr + corrected_dielectric
        end
    else
        seg_sigma_raw = seg_charge_uncorr ./ seg_area
    end

    return SigmaSurface(seg_area, seg_sigma_raw, seg_pos, energy_tot)
end
