# COSMO-RS.jl

COSMO-RS.jl is a in development Julia package for generating COSMO-based sigma profiles
and performing multi-state, multi-conformer thermodynamic averaging
for amino acids in solution.

The package is designed for:

- Zwitterionic systems
- Protonation-state equilibria
- Conformer ensemble averaging
- COSMO-RS preparation workflows
- Activity coefficient prediction in mixed solvents (e.g. water/urea/TMAO)

This package focuses on *physically rigorous statistical mechanics* and
clean reproducible computational workflows.

---

## Scientific Scope

This package implements:

1. Multi-state molecular ensembles (neutral, zwitterion, protonated, etc.)
2. Conformer generation (MoleculeFlow)
3. RMS pruning with Kabsch alignment
4. ORCA automation for COSMO surface generation
5. Gaussian sigma averaging (standard COSMO-RS preprocessing)
6. Boltzmann conformer weighting
7. State-level free energy weighting

The thermodynamic model used is:

Conformer partition function per state:

Z_s = Σ_i exp(−E_{s,i}/RT)

State free energy:

G_s = −RT ln(Z_s)

State population:

W_s = exp(−G_s/RT) / Σ_s' exp(−G_s'/RT)

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

## Installation

```julia
using Pkg
Pkg.add(path="path/to/COSMORS")
```

## Basic Workflow

### 1️⃣ Generate Molecular States
```julia
states = generate_states("glycine")
```

### 2️⃣ Generate Conformers
```julia
confs = generate_conformers(states[1]; nconfs=200)
```

### 3️⃣ RMS Prune
```julia
confs = prune_by_rms(confs; threshold=0.5)
```

### 4️⃣ Run ORCA COSMO Calculations
```julia
run_orca("conf.xyz", charge=0, multiplicity=1)
```

### 5️⃣ Parse Sigma Surface
```julia
surface = parse_sigma_surface("conf.orcacosmo")
```

### 6️⃣ Gaussian Sigma Averaging
```julia
gaussian_sigma_averaging_fast!(surface)
```

### 7️⃣ Build Sigma Profile
```julia
profile = compute_sigma_profile(surface)
```

### 8️⃣ Multi-State Averaging
```julia
ensemble = StateEnsemble(states[1], confs)
total_profile = average_over_states([ensemble])
```

## My Physical Assumptions:
- Energies must be computed in the same dielectric environment
used for COSMO surfaces.

- Gas-phase conformer energies should not be mixed with CPCM energies.

- Segment smoothing radius default: r_av = 0.5 Å

- RMS pruning default: 0.5 Å (heavy atoms)

## Planned Features

- PH-dependent protonation weighting
- Hydrogen-bond sigma moments

In a not so near future:
- Parallel ORCA job management 
- Multi-component mixture modeling 
- Activity coefficient calculation (γᵢ) - I wish :D 