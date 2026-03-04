# COSMO-RS.jl

COSMO-RS.jl is an in-development Julia package for generating COSMO-based sigma profiles
and performing multi-state, multi-conformer thermodynamic averaging
for amino acids in solution.

The package is designed for:

- Zwitterionic systems
- Protonation-state equilibria
- Conformer ensemble averaging
- COSMO-RS preparation workflows
- Activity coefficient prediction in mixed solvents (e.g. water/urea/TMAO)

---

## Scientific Scope

This package implements:

1. Multi-state molecular ensembles (neutral, zwitterion, protonated, etc.)
2. Conformer generation (MoleculeFlow) with RMS pruning and Kabsch alignment
3. ORCA automation for COSMO surface generation with COPT fallback
4. `.orcacosmo` file assembly and parsing (including CPCM-corrected charges)
5. Gaussian sigma averaging (standard COSMO-RS preprocessing)
6. Sigma profile computation via linear interpolation binning
7. Sigma moments (0th-6th order) with hydrogen-bond donor/acceptor decomposition
8. Boltzmann conformer weighting
9. State-level free energy weighting
10. pH-dependent protonation weighting via polyprotic acid equilibria

The thermodynamic model used is:

Conformer partition function per state:

Z_s = Σ_i exp(−E_{s,i}/RT)

State free energy:

G_s = −RT ln(Z_s)

State population (energy-only):

W_s = exp(−G_s/RT) / Σ_s' exp(−G_s'/RT)

pH-dependent state population:

W_s = α_j(pH) × w_s|j

where α_j is the macroscopic protonation fraction from polyprotic pKa equilibria
and w_s|j is the Boltzmann weight among tautomers at the same protonation level.

Total sigma profile:

σ_total = Σ_s W_s ( Σ_i w_{i|s} σ_{s,i} )

---

## Why This Exists

Amino acids require:

- Explicit zwitterionic treatment
- Solution-phase stabilization
- Multi-state thermodynamics
- Rigorous ensemble averaging

Many COSMO workflows:
- assume a single molecular structure
- ignore protonation equilibria
- ignore conformer entropy

This package addresses those limitations.

---

## Planned Features

- Parallel ORCA job management
- Multi-component mixture modeling
- Activity coefficient calculation
