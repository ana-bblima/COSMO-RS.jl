using Pkg
Pkg.activate(".")
using MoleculeFlow

include("../src/types.jl")
include("../src/states.jl")

gly = MolecularState("zwitterion",   "[NH3+]CC(=O)[O-]",     0, 1)
ala = MolecularState("zwitterion",   "CC([NH3+])C(=O)[O-]",  0, 1)

# dl-alpha-Amino n butyric acid
aaba = MolecularState("zwitterion",  "CC[C@H](C(=O)[O-])[NH3+]", 0, 1)

# dl-alpha-Amino-n valeric acid or dl-norvalina
anva = MolecularState("zwitterion",  "CCC[C@@H](C(=O)[O-])[NH3+]", 0, 1)

# a-Aminoisobutyric acid.
aaib = MolecularState("zwitterion",  "CC(C)(C(=O)[O-])[NH3+]", 0, 1)

val = MolecularState("zwitterion",   "CC(C)C([NH3+])C(=O)[O-]", 0, 1)

# beta-alanine
bala = MolecularState("zwitterion",  "C(C[NH3+])C(=O)[O-]", 0, 1)

# dl-beta-Aminobutyric acid
baba = MolecularState("zwitterion",  "CC(CC(=O)[O-])[NH3+]", 0, 1)

# dl-beta-Amino Valeric acid
bava = MolecularState("zwitterion",  "CCC(CC(=O)[O-])[NH3+]", 0, 1)

# y-Aminobutyric acid
gaba = MolecularState("zwitterion",  "C(CC(=O)[O-])C[NH3+]", 0, 1)

# dl-y-Amino Valeric acid
yava = MolecularState("zwitterion",  "CC(CCC(=O)[O-])[NH3+]", 0, 1)

# e-Aminocaproic acid
eaca = MolecularState("zwitterion",  "C(CCC(=O)[O-])CC[NH3+]",  0, 1)

digly = MolecularState("zwitterion", "[NH3+]CC(=O)NCC(=O)[O-]",  0, 1)

trigly = MolecularState("zwitterion", "[NH3+]CC(=O)NCC(=O)NCC(=O)[O-]",  0, 1)

# Glycyl-L-alanine 
gly_ala = MolecularState("zwitterion", "C[C@@H](C(=O)[O-])NC(=O)C[NH3+]",  0, 1)

# Alanylglycine
ala_gly = MolecularState("zwitterion", "CC(C(=O)NCC(=O)[O-])[NH3+]",  0, 1)

#ALANYLALANINE
ala_ala = MolecularState("zwitterion", "CC(C(=O)NC(C)C(=O)[O-])[NH3+]",  0, 1)

ser = MolecularState("zwitterion",   "[NH3+]C(CO)C(=O)[O-]", 0, 1)

thr = MolecularState("zwitterion",   "CC(O)C([NH3+])C(=O)[O-]", 0, 1)

# l-hydroxyproline
hpro = MolecularState("zwitterion",   "C1[C@@H](C[NH2+][C@@H]1C(=O)[O-])O", 0, 1)

# sarcosine
sar = MolecularState("zwitterion",   "C[NH2+]CC(=O)[O-]", 0, 1)

# betaine
bet = MolecularState("zwitterion",   "C[N+](C)(C)CC(=O)[O-]", 0, 1)


smiles = ["C[C@@H](C(=O)[O-])NC(=O)C[NH3+]", "CC(C(=O)NCC(=O)[O-])[NH3+]", "CC(C(=O)NC(C)C(=O)[O-])[NH3+]", "C1[C@@H](C[NH2+][C@@H]1C(=O)[O-])O", 
"C[NH2+]CC(=O)[O-]", "C[N+](C)(C)CC(=O)[O-]"]
mol =mol_from_smiles(smiles)


mols_to_grid_image(mol)

