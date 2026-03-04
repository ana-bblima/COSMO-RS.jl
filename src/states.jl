############
# AMINO ACID MOLECULAR STATES
############
#
# Each amino acid has at minimum four backbone states:
#   neutral:      NH2-CHR-COOH         (charge  0, multiplicity 1)
#   zwitterion:   NH3+-CHR-COO-        (charge  0, multiplicity 1)
#   protonated:   NH3+-CHR-COOH        (charge +1, multiplicity 1)
#   deprotonated: NH2-CHR-COO-         (charge -1, multiplicity 1)
#
# Amino acids with ionizable side chains carry additional states
# corresponding to side-chain (de)protonation.
############

# Three-letter code -> full name
const AMINO_ACID_THREE_TO_FULL = Dict{String,String}(
    "gly" => "glycine",
    "ala" => "alanine",
    "val" => "valine",
    "leu" => "leucine",
    "ile" => "isoleucine",
    "pro" => "proline",
    "phe" => "phenylalanine",
    "trp" => "tryptophan",
    "met" => "methionine",
    "ser" => "serine",
    "thr" => "threonine",
    "cys" => "cysteine",
    "tyr" => "tyrosine",
    "asn" => "asparagine",
    "gln" => "glutamine",
    "asp" => "aspartic_acid",
    "glu" => "glutamic_acid",
    "lys" => "lysine",
    "arg" => "arginine",
    "his" => "histidine",
)

# One-letter code -> full name
const AMINO_ACID_ONE_TO_FULL = Dict{String,String}(
    "g" => "glycine",
    "a" => "alanine",
    "v" => "valine",
    "l" => "leucine",
    "i" => "isoleucine",
    "p" => "proline",
    "f" => "phenylalanine",
    "w" => "tryptophan",
    "m" => "methionine",
    "s" => "serine",
    "t" => "threonine",
    "c" => "cysteine",
    "y" => "tyrosine",
    "n" => "asparagine",
    "q" => "glutamine",
    "d" => "aspartic_acid",
    "e" => "glutamic_acid",
    "k" => "lysine",
    "r" => "arginine",
    "h" => "histidine",
)

const AMINO_ACID_STATES = Dict{String, Vector{MolecularState}}(

    # ──────────────────────────────────────────────────────────
    # Non-polar, aliphatic side chains
    # ──────────────────────────────────────────────────────────

    "glycine" => [
        MolecularState("neutral",      "NCC(=O)O",             0, 1),
        MolecularState("zwitterion",   "[NH3+]CC(=O)[O-]",     0, 1),
        MolecularState("protonated",   "[NH3+]CC(=O)O",        1, 1),
        MolecularState("deprotonated", "NCC(=O)[O-]",         -1, 1),
    ],

    "alanine" => [
        MolecularState("neutral",      "CC(N)C(=O)O",                 0, 1),
        MolecularState("zwitterion",   "CC([NH3+])C(=O)[O-]",         0, 1),
        MolecularState("protonated",   "CC([NH3+])C(=O)O",            1, 1),
        MolecularState("deprotonated", "CC(N)C(=O)[O-]",             -1, 1),
    ],

    "valine" => [
        MolecularState("neutral",      "CC(C)C(N)C(=O)O",            0, 1),
        MolecularState("zwitterion",   "CC(C)C([NH3+])C(=O)[O-]",    0, 1),
        MolecularState("protonated",   "CC(C)C([NH3+])C(=O)O",       1, 1),
        MolecularState("deprotonated", "CC(C)C(N)C(=O)[O-]",        -1, 1),
    ],

    "leucine" => [
        MolecularState("neutral",      "CC(C)CC(N)C(=O)O",           0, 1),
        MolecularState("zwitterion",   "CC(C)CC([NH3+])C(=O)[O-]",   0, 1),
        MolecularState("protonated",   "CC(C)CC([NH3+])C(=O)O",      1, 1),
        MolecularState("deprotonated", "CC(C)CC(N)C(=O)[O-]",       -1, 1),
    ],

    "isoleucine" => [
        MolecularState("neutral",      "CCC(C)C(N)C(=O)O",           0, 1),
        MolecularState("zwitterion",   "CCC(C)C([NH3+])C(=O)[O-]",   0, 1),
        MolecularState("protonated",   "CCC(C)C([NH3+])C(=O)O",      1, 1),
        MolecularState("deprotonated", "CCC(C)C(N)C(=O)[O-]",       -1, 1),
    ],

    "proline" => [
        # Proline has a secondary amine (pyrrolidine ring)
        MolecularState("neutral",      "O=C(O)C1CCCN1",              0, 1),
        MolecularState("zwitterion",   "O=C([O-])C1CCC[NH2+]1",      0, 1),
        MolecularState("protonated",   "O=C(O)C1CCC[NH2+]1",         1, 1),
        MolecularState("deprotonated", "O=C([O-])C1CCCN1",          -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Aromatic side chains
    # ──────────────────────────────────────────────────────────

    "phenylalanine" => [
        MolecularState("neutral",      "NC(Cc1ccccc1)C(=O)O",             0, 1),
        MolecularState("zwitterion",   "[NH3+]C(Cc1ccccc1)C(=O)[O-]",     0, 1),
        MolecularState("protonated",   "[NH3+]C(Cc1ccccc1)C(=O)O",        1, 1),
        MolecularState("deprotonated", "NC(Cc1ccccc1)C(=O)[O-]",         -1, 1),
    ],

    "tryptophan" => [
        MolecularState("neutral",      "NC(Cc1c[nH]c2ccccc12)C(=O)O",          0, 1),
        MolecularState("zwitterion",   "[NH3+]C(Cc1c[nH]c2ccccc12)C(=O)[O-]",  0, 1),
        MolecularState("protonated",   "[NH3+]C(Cc1c[nH]c2ccccc12)C(=O)O",     1, 1),
        MolecularState("deprotonated", "NC(Cc1c[nH]c2ccccc12)C(=O)[O-]",      -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Sulfur-containing side chains
    # ──────────────────────────────────────────────────────────

    "methionine" => [
        MolecularState("neutral",      "CSCCC(N)C(=O)O",              0, 1),
        MolecularState("zwitterion",   "CSCCC([NH3+])C(=O)[O-]",      0, 1),
        MolecularState("protonated",   "CSCCC([NH3+])C(=O)O",         1, 1),
        MolecularState("deprotonated", "CSCCC(N)C(=O)[O-]",          -1, 1),
    ],

    "cysteine" => [
        # Side chain: -CH2-SH   (pKa ~8.3 for thiol)
        MolecularState("neutral",      "NC(CS)C(=O)O",                0, 1),
        MolecularState("zwitterion",   "[NH3+]C(CS)C(=O)[O-]",        0, 1),
        MolecularState("protonated",   "[NH3+]C(CS)C(=O)O",           1, 1),
        MolecularState("deprotonated", "NC(CS)C(=O)[O-]",            -1, 1),
        # Thiolate: zwitterion backbone + deprotonated thiol (relevant above pH ~8)
        MolecularState("thiolate",     "[NH3+]C(C[S-])C(=O)[O-]",    -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Hydroxyl side chains
    # ──────────────────────────────────────────────────────────

    "serine" => [
        MolecularState("neutral",      "NC(CO)C(=O)O",                0, 1),
        MolecularState("zwitterion",   "[NH3+]C(CO)C(=O)[O-]",        0, 1),
        MolecularState("protonated",   "[NH3+]C(CO)C(=O)O",           1, 1),
        MolecularState("deprotonated", "NC(CO)C(=O)[O-]",            -1, 1),
    ],

    "threonine" => [
        MolecularState("neutral",      "CC(O)C(N)C(=O)O",            0, 1),
        MolecularState("zwitterion",   "CC(O)C([NH3+])C(=O)[O-]",    0, 1),
        MolecularState("protonated",   "CC(O)C([NH3+])C(=O)O",       1, 1),
        MolecularState("deprotonated", "CC(O)C(N)C(=O)[O-]",        -1, 1),
    ],

    "tyrosine" => [
        # Side chain: -CH2-C6H4-OH   (pKa ~10.1 for phenol)
        MolecularState("neutral",      "NC(Cc1ccc(O)cc1)C(=O)O",              0, 1),
        MolecularState("zwitterion",   "[NH3+]C(Cc1ccc(O)cc1)C(=O)[O-]",      0, 1),
        MolecularState("protonated",   "[NH3+]C(Cc1ccc(O)cc1)C(=O)O",         1, 1),
        MolecularState("deprotonated", "NC(Cc1ccc(O)cc1)C(=O)[O-]",          -1, 1),
        # Phenolate: zwitterion backbone + deprotonated phenol (relevant above pH ~10)
        MolecularState("phenolate",    "[NH3+]C(Cc1ccc([O-])cc1)C(=O)[O-]",  -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Amide side chains
    # ──────────────────────────────────────────────────────────

    "asparagine" => [
        MolecularState("neutral",      "NC(CC(=O)N)C(=O)O",               0, 1),
        MolecularState("zwitterion",   "[NH3+]C(CC(=O)N)C(=O)[O-]",       0, 1),
        MolecularState("protonated",   "[NH3+]C(CC(=O)N)C(=O)O",          1, 1),
        MolecularState("deprotonated", "NC(CC(=O)N)C(=O)[O-]",           -1, 1),
    ],

    "glutamine" => [
        MolecularState("neutral",      "NC(CCC(=O)N)C(=O)O",              0, 1),
        MolecularState("zwitterion",   "[NH3+]C(CCC(=O)N)C(=O)[O-]",      0, 1),
        MolecularState("protonated",   "[NH3+]C(CCC(=O)N)C(=O)O",         1, 1),
        MolecularState("deprotonated", "NC(CCC(=O)N)C(=O)[O-]",          -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Acidic side chains
    # ──────────────────────────────────────────────────────────

    "aspartic_acid" => [
        # Side chain: -CH2-COOH   (pKa ~3.7)
        MolecularState("neutral",                "NC(CC(=O)O)C(=O)O",                  0, 1),
        MolecularState("zwitterion",             "[NH3+]C(CC(=O)O)C(=O)[O-]",          0, 1),
        MolecularState("protonated",             "[NH3+]C(CC(=O)O)C(=O)O",             1, 1),
        MolecularState("deprotonated",           "NC(CC(=O)O)C(=O)[O-]",              -1, 1),
        # Both carboxyls deprotonated (physiological pH ~7.4)
        MolecularState("doubly_deprotonated",    "[NH3+]C(CC(=O)[O-])C(=O)[O-]",      -1, 1),
    ],

    "glutamic_acid" => [
        # Side chain: -CH2CH2-COOH   (pKa ~4.1)
        MolecularState("neutral",                "NC(CCC(=O)O)C(=O)O",                 0, 1),
        MolecularState("zwitterion",             "[NH3+]C(CCC(=O)O)C(=O)[O-]",         0, 1),
        MolecularState("protonated",             "[NH3+]C(CCC(=O)O)C(=O)O",            1, 1),
        MolecularState("deprotonated",           "NC(CCC(=O)O)C(=O)[O-]",             -1, 1),
        # Both carboxyls deprotonated (physiological pH ~7.4)
        MolecularState("doubly_deprotonated",    "[NH3+]C(CCC(=O)[O-])C(=O)[O-]",     -1, 1),
    ],

    # ──────────────────────────────────────────────────────────
    # Basic side chains
    # ──────────────────────────────────────────────────────────

    "lysine" => [
        # Side chain: -(CH2)4-NH2   (pKa ~10.5 for side-chain amine)
        MolecularState("neutral",             "NC(CCCCN)C(=O)O",                    0, 1),
        MolecularState("zwitterion",          "[NH3+]C(CCCCN)C(=O)[O-]",            0, 1),
        MolecularState("protonated",          "[NH3+]C(CCCCN)C(=O)O",               1, 1),
        MolecularState("deprotonated",        "NC(CCCCN)C(=O)[O-]",                -1, 1),
        # Both amines protonated (physiological pH ~7.4)
        MolecularState("doubly_protonated",   "[NH3+]C(CCCC[NH3+])C(=O)[O-]",       1, 1),
    ],

    "arginine" => [
        # Side chain: -(CH2)3-NHC(=NH)NH2   (pKa ~12.5 for guanidinium)
        MolecularState("neutral",       "NC(CCCNC(=N)N)C(=O)O",                    0, 1),
        MolecularState("zwitterion",    "[NH3+]C(CCCNC(=N)N)C(=O)[O-]",            0, 1),
        MolecularState("protonated",    "[NH3+]C(CCCNC(=N)N)C(=O)O",               1, 1),
        MolecularState("deprotonated",  "NC(CCCNC(=N)N)C(=O)[O-]",                -1, 1),
        # Guanidinium protonated (physiological pH ~7.4)
        MolecularState("guanidinium",   "[NH3+]C(CCCNC(=[NH2+])N)C(=O)[O-]",       1, 1),
    ],

    "histidine" => [
        # Side chain: -CH2-(imidazol-4-yl)  (pKa ~6.0 for imidazole)
        # Two neutral tautomers: Nε-H and Nδ-H
        MolecularState("neutral_epsilon",  "NC(Cc1c[nH]cn1)C(=O)O",                   0, 1),
        MolecularState("neutral_delta",    "NC(Cc1cnc[nH]1)C(=O)O",                   0, 1),
        MolecularState("zwitterion",       "[NH3+]C(Cc1c[nH]cn1)C(=O)[O-]",           0, 1),
        MolecularState("protonated",       "[NH3+]C(Cc1c[nH]cn1)C(=O)O",              1, 1),
        MolecularState("deprotonated",     "NC(Cc1c[nH]cn1)C(=O)[O-]",               -1, 1),
        # Imidazolium: both ring nitrogens protonated (relevant near pH 6)
        MolecularState("imidazolium",      "[NH3+]C(Cc1c[nH]c[nH+]1)C(=O)[O-]",      1, 1),
    ],
)

"""
    AMINO_ACID_NAMES

Sorted list of all supported amino acid names (full lowercase names).
"""
const AMINO_ACID_NAMES = sort(collect(keys(AMINO_ACID_STATES)))


############
# EXPERIMENTAL pKa DATA
############
#
# pKa values (sorted ascending) and charge of the fully deprotonated form.
# Sources: CRC Handbook of Chemistry and Physics, Lehninger Principles of Biochemistry.
#
# For amino acids with 2 ionizable groups (backbone only):
#   pKa1 = α-COOH,  pKa2 = α-NH3⁺
#   charge_fully_deprotonated = -1  (NH2-CHR-COO⁻)
#
# For amino acids with 3 ionizable groups (backbone + side chain):
#   pKa values sorted ascending across all groups.
#   charge_fully_deprotonated depends on side-chain type:
#     acidic side chain (Asp, Glu) or Cys/Tyr:  -2
#     basic side chain (Lys, Arg, His):          -1
############

const AMINO_ACID_PKA = Dict{String, AminoAcidPKa}(

    # ── Non-polar, aliphatic ──────────────────────────────────
    "glycine"       => AminoAcidPKa([2.34, 9.60],        -1),
    "alanine"       => AminoAcidPKa([2.34, 9.69],        -1),
    "valine"        => AminoAcidPKa([2.32, 9.62],        -1),
    "leucine"       => AminoAcidPKa([2.36, 9.60],        -1),
    "isoleucine"    => AminoAcidPKa([2.36, 9.68],        -1),
    "proline"       => AminoAcidPKa([1.99, 10.60],       -1),

    # ── Aromatic ──────────────────────────────────────────────
    "phenylalanine" => AminoAcidPKa([1.83, 9.13],        -1),
    "tryptophan"    => AminoAcidPKa([2.83, 9.39],        -1),

    # ── Sulfur-containing ─────────────────────────────────────
    "methionine"    => AminoAcidPKa([2.28, 9.21],        -1),
    "cysteine"      => AminoAcidPKa([1.96, 8.18, 10.28], -2),

    # ── Hydroxyl ──────────────────────────────────────────────
    "serine"        => AminoAcidPKa([2.21, 9.15],        -1),
    "threonine"     => AminoAcidPKa([2.09, 9.10],        -1),
    "tyrosine"      => AminoAcidPKa([2.20, 9.11, 10.07], -2),

    # ── Amide ─────────────────────────────────────────────────
    "asparagine"    => AminoAcidPKa([2.02, 8.80],        -1),
    "glutamine"     => AminoAcidPKa([2.17, 9.13],        -1),

    # ── Acidic side chains ────────────────────────────────────
    "aspartic_acid" => AminoAcidPKa([1.88, 3.65, 9.60],  -2),
    "glutamic_acid" => AminoAcidPKa([2.19, 4.25, 9.67],  -2),

    # ── Basic side chains ─────────────────────────────────────
    "lysine"        => AminoAcidPKa([2.18, 8.95, 10.53], -1),
    "arginine"      => AminoAcidPKa([2.17, 9.04, 12.48], -1),
    "histidine"     => AminoAcidPKa([1.82, 6.00, 9.17],  -1),
)


"""
    generate_states(name::String) -> Vector{MolecularState}

Return the molecular states for the given amino acid.
Accepts full names ("glycine"), three-letter codes ("gly"),
or one-letter codes ("G"). Lookup is case-insensitive.
"""
function generate_states(name::String)

    key = lowercase(strip(name))

    # Direct lookup by full name
    if haskey(AMINO_ACID_STATES, key)
        return AMINO_ACID_STATES[key]
    end

    # Three-letter code
    if haskey(AMINO_ACID_THREE_TO_FULL, key)
        return AMINO_ACID_STATES[AMINO_ACID_THREE_TO_FULL[key]]
    end

    # One-letter code
    if haskey(AMINO_ACID_ONE_TO_FULL, key)
        return AMINO_ACID_STATES[AMINO_ACID_ONE_TO_FULL[key]]
    end

    error("State definitions not found for '$name'. " *
          "Available amino acids: $(join(AMINO_ACID_NAMES, ", "))")
end
