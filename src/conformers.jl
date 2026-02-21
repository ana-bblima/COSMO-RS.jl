using MoleculeFlow

function generate_conformers(state::MolecularState;
                             nconfs=100)

    mol = MoleculeFlow.from_smiles(state.smiles)
    mol3d = MoleculeFlow.embed(mol)

    conformers = MoleculeFlow.generate_conformers(
        mol3d;
        n_conformers=nconfs,
        method=:etkdg
    )

    return conformers
end

# Need to add the RMS filtering and energy pruning