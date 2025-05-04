
# 🌱 Vegetation Community Patterns in Drylands — pde2path Implementation

This repository contains the `pde2path` implementation of a trait-structured vegetation model under dryland conditions, including a multi-species community model and a simplified single-species case. The code enables the computation of periodic solution branches, bifurcation diagrams, and Busse-Balloon structures.

---

## 🔧 Folder Structure

```
bwhRS/
├── community_model/
│   ├── sG.m               # RHS of community model equations
│   ├── sGjac.m            # Analytical Jacobian
│   ├── oosetfemops.m      # FEM matrix setup
│   └── bwhinit.m          # Model initialization
├── single_species_model/
│   ├── sG.m               # RHS for single species
│   ├── sGjac.m            # Jacobian for single species
│   ├── oosetfemops.m
│   └── bwhinit.m
├── cmds/
│   ├── cmds1.m            # Continuation: homogeneous + Turing branches
│   └── cmds2.m            # Brute-force Busse-Balloon computation
├── utils/
│   └── bfbb.m             # Stability range extraction from branches
└── results/
    └── comm/              # Output directories (e.g. hom, T1, T2, ...)
```

---

## 🚀 Getting Started

```matlab
addpath(genpath('bwhRS'));
cmds1; % Runs the main continuation for the community model
cmds2; % Computes the Busse-Balloon by extracting stability info
```

**Requirements**:
- MATLAB (R2020 or later)
- [pde2path](https://www.staff.uni-oldenburg.de/hannes.uecker/pde2path/)
- OOPDE library (included or linked from pde2path)

---

## 📊 Parameter Table

| Index | Parameter | Description                | Symbol | Typical Value |
|-------|-----------|----------------------------|--------|----------------|
| 1     | `pp`      | Precipitation              | \(P\)  | 300            |
| 2     | `Lam0`    | Base growth rate           | \(\Lambda_0\) | 8        |
| 3     | `Ga`      | Mortality scale            | \(\Gamma\)     | 10       |
| 4     | `A`       | Max infiltration rate      | \(A\)  | 3000           |
| 5     | `R`       | Evaporation feedback       | \(R\)  | 0.7            |
| 6     | `L0`      | Base evaporation rate      | \(L_0\)| 200            |
| ...   | ...       | ...                        | ...    | ...            |
| 19    | `chimin`  | Min trait value            | \(\chi_{\min}\) | 0       |
| 20    | `chimax`  | Max trait value            | \(\chi_{\max}\) | 1       |

> See `cmds1.m` line 10 for full `par` vector construction.

---

## 📎 References

- Main paper: [FPBUM24], *submitted/in preparation*
- For `pde2path`: Uecker (2021), Rademacher & Uecker (2018)


---

## 🧱 About OOPDE

**OOPDE** (*Object-Oriented Partial Differential Equations*) is a MATLAB framework for finite element discretization and PDE handling. It provides the underlying data structures used in `pde2path`, offering:

- Object-oriented setup of FEM problems
- Efficient assembly of mass and stiffness matrices
- Handling of boundary conditions and domain structures
- Support for mesh refinement and matrix reuse

Unlike MATLAB’s PDE Toolbox, OOPDE gives users full access to the internals of FEM discretization, making it ideal for custom continuation, bifurcation tracking, and PDE analysis.

**Note**: OOPDE is bundled with `pde2path` and does not require separate installation.
