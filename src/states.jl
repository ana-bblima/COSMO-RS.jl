# I will be extending for all amino acids soon 

function generate_states(name::String)

    if lowercase(name) == "glycine"
        return [
            MolecularState("neutral",
                "NCC(=O)O",
                0,
                1),

            MolecularState("zwitterion",
                "[NH3+]CC(=O)[O-]",
                0,
                1)
        ]
    end

    error("State definitions not implemented for $name")
end