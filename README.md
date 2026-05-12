# EUCLID — Linear Viscoelasticity 2D

**Numerical Implementation and Experimental Validation of the EUCLID Framework for Linear Viscoelastic Parameter Identification**

*Master's Thesis - Kali Satya Sri Charan Undamatla*

------

> **TL;DR:** A Python implementation of the EUCLID framework for identifying linear viscoelastic Prony series parameters from full-field displacement measurements. Supports plain and complex specimen geometries, with synthetic validation achieving **< 1% parameter recovery error**.

![EUCLID Identification Pipeline](figures/EUCLID_identification_Pipeline.png)

## Overview

This repository contains the Python implementation developed in the thesis for identifying linear viscoelastic Prony series parameters from full-field displacement measurements within the **EUCLID** (Efficient Unsupervised Constitutive Law Identification) framework.

The identification problem is formulated in **plane stress** for thin flat specimens under uniaxial tensile creep, supporting plain rectangles and arbitrary geometries with holes. The viscoelastic response is modeled using a dense **Prony-series library** with prescribed relaxation times, and the unknown moduli are determined through a **Non-Negative Least Squares (NNLS)** inverse problem that enforces thermodynamic admissibility.

Two complementary data paths are supported:

- **Synthetic path** : Forward FEM simulations with known ground-truth parameters θ_true generate benchmark datasets for validation
- **Experimental path** : Raw Digital Image Correlation (DIC) displacement fields are preprocessed into FE-compatible inputs

Both paths converge to a **Standardized FE Dataset** D = {**X**, conn, **U**, **F**, t} that serves as the unified input to the inverse solver.

---

## Results at a Glance

### Forward Solver: Complex Three-Hole Specimen

The forward FEM solver generates standardized datasets for arbitrary specimen geometries. Shown below: 22.83 × 59.62 mm specimen with three circular holes, meshed via pygmsh/OpenCASCADE.

**FE Mesh:**

![FE Mesh](figures/mesh.png)

**Simulated displacement fields (Ux, Uy) at the final time step:**

![Displacement Fields](figures/displacement_field.png)

### Inverse Identification: Parameter Recovery Accuracy

From Table 5.1 of the thesis (synthetic plain-rectangle specimen, N_G = N_K = 150 Prony branches):

| Parameter | Ground Truth | Relative Error |
|-----------|-------------|----------------|
| G₀ (instantaneous shear) | 3200 MPa | 0.21% |
| G∞ (equilibrium shear)   | 1500 MPa | 1.00% |
| K₀ (instantaneous bulk)  | 3767 MPa | 0.87% |
| K∞ (equilibrium bulk)    | 2000 MPa | 0.99% |

Accuracy is consistent across plain, elliptical-hole, and three-hole geometries (Table 5.3, thesis), and stable across mesh densities from 280 to 3888 elements (Table 5.2, thesis).

![Parameter Recovery](figures/Parameter_recovery.png)

**Identified Prony spectrum (NNLS + clustering, deviatoric and volumetric branches):**

![Prony Series](figures/Prony_series.png)

### Experimental Validation: Three-Hole PA6 Specimen

Measured vs EUCLID-predicted vertical displacement field (u_y) for a real PA6 creep specimen with three holes (Fig. 5.6, thesis). DIC data collected at Fraunhofer IWM.

![Three-hole field comparison](figures/Three_hole%20field%20comparision.png)

---

## Repository Structure

```
EUCLID_Lin_viscoelasticity_2D/
│
├── Forward_solver/              # FEM forward problem: synthetic dataset generation
│   ├── core/
│   │   ├── material.py          # Viscoelastic material: known Prony parameters (G, K, tau)
│   │   ├── mesh.py              # Structured CST mesh generation for rectangular specimens
│   │   ├── geometry_builder.py  # Complex geometry (holes) via pygmsh/OpenCASCADE
│   │   ├── geometry_advanced.py # Parametric geometry: true ellipses, smooth curves
│   │   ├── mesh_converter.py    # pygmsh mesh → EUCLID coord.csv / conne.txt format
│   │   ├── assembly.py          # Global stiffness K^n and history force F_hist assembly
│   │   ├── time_integration.py  # Trapezoidal update of Maxwell internal variables beta
│   │   ├── solver.py            # Forward FEM time-stepping: K^n U^n = F_ext + F_hist
│   │   └── data_generation.py   # Export results as D = {X, conn, U, F, t}
│   │
│   ├── configs/                 # Geometry configs for complex specimens (holes)
│   ├── scripts/                 # Internal: run_full_simulation.py (full pipeline implementation)
│   ├── run_simulation.py        # MAIN: run forward FEM → writes standardized FE dataset
│   └── run_verification.py      # Verification against known solutions
│
├── inverse_problem/             # EUCLID inverse problem: identify Prony parameters
│   ├── core/
│   │   ├── data.py              # Load and validate standardized FE dataset D
│   │   ├── geometry.py          # Mesh structure: nodes, CST elements, B-matrix (B_e)
│   │   ├── material.py          # Prony structure: log-spaced tau grid, projectors Pi_D, Pi_V
│   │   ├── history.py           # Container for beta history variables (beta_{e,a}, beta_{th,a})
│   │   ├── beta_computation.py  # Compute beta coefficients via trapezoidal time integration
│   │   ├── assembly.py          # Assemble global inverse system matrix A and RHS r
│   │   ├── boundary.py          # BC strategies (TopBottomForce / BottomForceBC);
│   │   │                        #   interior (lambda_int) and boundary (lambda_bnd) weights
│   │   ├── solver.py            # NNLS solver: min ||A*theta - r||^2 s.t. theta >= 0
│   │   ├── solver_lasso.py      # Alternative LASSO-regularised solver
│   │   ├── clustering.py        # Merge closely-spaced relaxation times (weighted average)
│   │   └── visualization.py     # 9 diagnostic plots: data, mesh, beta, residuals, spectra
│   │
│   ├── preprocessing/
│   │   ├── unified_preprocessor.py  # DIC point cloud → standardized FE dataset
│   │   └── dic_data_analysis.py     # DIC data quality and mapping analysis
│   │
│   ├── Postprocessing/          # Output directory (results.npz + figures written here)
│   ├── Real_data/               # Place raw DIC experiment data here
│   ├── run_experiment.py        # MAIN: run inverse identification pipeline
│   ├── inverse_problem.py       # Orchestrator: history → assembly → NNLS → clustering
│   ├── plot_relaxation_curves.py  # G(t), K(t) identified vs ground truth
│   ├── plot_prony_comparison.py   # Bar chart: identified Prony terms vs ground truth
│   └── read_results_npz.py        # Inspect saved results.npz
│
├── standard_FE_dataset/         # Standardized FE datasets D = {X, conn, U, F, t}
│   └── 8xx/                     # Synthetic experiment folders (generated by Forward_solver)
│                                # Experimental datasets (9xx-series) not included — see Data Availability
│
├── figures/                     # Result and documentation figures (PNG exports from thesis)
├── requirements.txt
└── README.md
```

---

## Standardized FE Dataset

Both synthetic and experimental data are stored in the same 5-tuple format D = {**X**, conn, **U**, **F**, **t**} (Section 4.3.1 of the thesis):

| File | Symbol | Shape | Description |
|------|--------|-------|-------------|
| `coord.csv` | **X** | N_n × 2 | Nodal coordinates (x_i, y_i) with boundary labels |
| `conne.txt` | conn | N_e × 3 | CST element connectivity (1-based node indices) |
| `U.csv` | **U** | 2N_n × N_t | Nodal displacement history: [u_x(t_1)...u_y(t_Nt)] [mm] |
| `F.csv` | **F** | 4 × N_t | Boundary forces [F_top^y, F_bottom^y, F_left^x, F_right^x] [kN] |
| `time.csv` | **t** | N_t | Monotonically increasing time discretization [s] |

---

## Installation

**Python 3.8+ required**

```bash
pip install -r requirements.txt
```

Verify:
```bash
cd Forward_solver
python -c "from core.material import ViscoelasticMaterial; print('Forward solver ready')"

cd ../inverse_problem
python -c "from inverse_problem import InverseProblem; print('Inverse solver ready')"
```

---

## Workflow 1: Synthetic Data (Forward → Inverse)

### Step 1: Generate Synthetic Dataset

```bash
cd Forward_solver
python run_simulation.py
```

Configure in `run_simulation.py`:

```python
experiment_id = 1        # Any integer — sets the output folder name in standard_FE_dataset/
nx, ny        = 5, 13    # Structured mesh density (nodes)
n_timesteps   = 600      # Number of time steps N_t
dt            = 0.5      # Time step h_n [s]
load          = 50.0     # Distributed traction q [N/mm]
```

The forward solver integrates the viscoelastic FEM system at each time step (Eq. 3.50):

```
K^n U^n = F_ext^n + F_hist^{n-1}
```

**Output:** `standard_FE_dataset/<experiment_id>/` with coord.csv, conne.txt, U.csv, F.csv, time.csv, ground_truth_material.txt, SUMMARY.txt

### Step 2: Run Inverse Identification

```bash
cd inverse_problem
python run_experiment.py
```

Configure in `run_experiment.py`:

```python
EXPERIMENT_NUMBER = 1     # Must match the experiment_id used in run_simulation.py
N_MAXWELL_SHEAR   = 150   # N_G: deviatoric Maxwell branches
N_MAXWELL_BULK    = 150   # N_K: volumetric Maxwell branches
TAU_MIN, TAU_MAX  = 1.0, 600.0   # Relaxation time library range [s]
LAMBDA_INTERIOR   = 1.0   # Interior equilibrium weight (lambda_int)
LAMBDA_BOUNDARY   = 1.0   # Boundary traction weight (lambda_bnd)
APPLY_CLUSTERING  = True  # Merge closely-spaced relaxation times
```

The solver constructs and solves the global NNLS system (Eq. 3.56):

```
min_{theta >= 0}  (1/2) || A * theta - r ||^2

theta = [G_inf, G_1, ..., G_NG, K_inf, K_1, ..., K_NK]
```

**Output:** `inverse_problem/Postprocessing/final_outputs/<EXPERIMENT_NUMBER>_.../results.npz` + 9 diagnostic figures

### Step 3: Plot Results

```bash
cd inverse_problem
python plot_relaxation_curves.py    # G(t), K(t) vs ground truth
python plot_prony_comparison.py     # Bar chart: identified vs ground truth
```

---

## Workflow 2: Experimental DIC Data (Preprocess → Inverse)

### Step 1: Preprocess DIC Data

Place raw DIC CSV files (exported from GOM Correlate) in `inverse_problem/Real_data/experiments/`.

```bash
cd inverse_problem/preprocessing
python unified_preprocessor.py
```

The preprocessor (Section 4.3.2 of the thesis) performs:
1. **Temporal consistency:** permanent DIC node selection across all frames
2. **Nearest-neighbor mapping:** DIC point cloud to FE mesh (mean mapping distance ~0.22 mm)
3. **Rigid-body motion removal:** subtraction of reference node translation
4. **Boundary labeling:** automatic classification of Top / Bottom / Left / Right nodes
5. **Export** → `standard_FE_dataset/<experiment_id>/` in standardized D format

![DIC Preprocessing Pipeline](figures/preprocessing_904.png)

### Step 2: Run Inverse Identification

```bash
cd inverse_problem
python run_experiment.py
```

Set for experimental data:

```python
EXPERIMENT_NUMBER = 1     # Must match the folder created by the preprocessor
LAMBDA_INTERIOR   = 1.0   # Use interior equilibrium equations for real data
LAMBDA_BOUNDARY   = 1.0
TAU_MAX           = 1206.0  # Match experimental creep duration [s]
```

---

## Scientific Background

### Generalized Maxwell Model (Prony Series)

The plane-stress deviatoric and volumetric relaxation functions are (Eq. 3.32):

```
G(t) = G_inf  +  sum_{a=1}^{N_G}  G_a * exp(-t / tau_Ga)    [Shear, MPa]
K(t) = K_inf  +  sum_{a=1}^{N_K}  K_a * exp(-t / tau_Ka)    [Volumetric, MPa]
```

Instantaneous moduli: G_0 = G_inf + sum(G_a),  K_0 = K_inf + sum(K_a)

### EUCLID Inverse System

For fixed relaxation times {tau_Ga}, {tau_Ka}, the time-discrete constitutive update is **linear** in the unknown moduli. Aggregating equilibrium equations across all N_t time steps yields (Eq. 3.55):

```
A * theta ≈ r
```

Solved as a non-negative constrained optimization (Eq. 3.56):

```
min_{theta >= 0}  || A * theta - r ||^2        (NNLS)
```

Non-negativity enforces thermodynamic admissibility of the identified relaxation spectrum.

![Algorithmic Flowchart](figures/Algorithmic_Flowchart.png)

### Ground-Truth Material: Synthetic Validation (Table 4.2)

```
G_prony = [200, 500, 1000] MPa,  tau_G = [5.3, 50.1, 400.2] s,  G_inf = 1500 MPa
K_prony = [500, 700,  567] MPa,  tau_K = [5.3, 50.1, 400.2] s,  K_inf = 2000 MPa
```

Synthetic validation achieves approximately **1% relative error** across multiple decades of relaxation time (Section 5.2 of the thesis).

---

## Experimental Datasets

PA6 (Polyamide 6) tensile creep specimens tested at Fraunhofer IWM:
- Uniaxial creep, pinned-and-clamped grips, loading direction +y
- DIC measurements at ~1 Hz using GOM Correlate software
- Specimen dimensions: 90 mm total length, 24 mm width, 1 mm thickness

> The experimental datasets (plain, elliptical-hole, and three-hole PA6 specimens) are excluded from this public repository. See the [Data Availability](#data-availability-and-reproducibility) section below.

---

## Key Parameters Reference

### Forward Solver (`run_simulation.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `experiment_id` | any int | Output folder name in standard_FE_dataset/ |
| `width / height` | 20 / 60 mm | Specimen gauge dimensions |
| `nx / ny` | 5 / 13 | Structured mesh node count |
| `dt` | 0.5 s | Time step h_n |
| `n_timesteps` | 600 | Total time steps N_t |
| `load` | 50 N/mm | Distributed traction q |

### Inverse Solver (`run_experiment.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `EXPERIMENT_NUMBER` | any int | Must match `experiment_id` used in Forward_solver |
| `N_MAXWELL_SHEAR/BULK` | 150 | Prony branch counts N_G, N_K |
| `TAU_MIN / TAU_MAX` | 1 / 600 s | Relaxation time library range |
| `LAMBDA_INTERIOR` | 1.0 | Interior equilibrium weight |
| `LAMBDA_BOUNDARY` | 1.0 | Boundary traction weight |
| `APPLY_CLUSTERING` | True | Merge closely-spaced tau values |
| `CLUSTERING_RANGE` | 0.3 | Merging threshold (30% relative distance) |

**Recommended settings:**
- Synthetic data: `LAMBDA_INTERIOR = 1.0, LAMBDA_BOUNDARY = 1.0`
- Real DIC data: `LAMBDA_INTERIOR = 1.0, LAMBDA_BOUNDARY = 1.0`

---

## Output Files

### `results.npz`

```python
import numpy as np
data = np.load('results.npz', allow_pickle=True)

# Full Prony spectrum (all N_G / N_K terms including zeros)
G_params = data['G_params']       # [N_G] deviatoric moduli [MPa]
tau_G    = data['tau_G']          # [N_G] relaxation times [s]

# Active terms after clustering
G_nonzero     = data['G_nonzero']
tau_G_nonzero = data['tau_G_nonzero']

# Equilibrium moduli
G_inf = data['G_inf']   # [MPa]
K_inf = data['K_inf']   # [MPa]
```

### Diagnostic Plots (9 figures)

| # | Plot | Description |
|---|------|-------------|
| 1 | Data validation | Input displacement field u(x,t) and force history F(t) |
| 2 | Mesh quality | Node distribution and CST element connectivity |
| 3 | Beta fields | History variables β at selected time steps |
| 4 | Displacement residuals | Measured vs reconstructed u |
| 5 | Force residuals | Boundary force balance check |
| 6 | Parameter distribution | All N_G + N_K identified Prony moduli |
| 7 | Clustering comparison | Relaxation spectrum before and after clustering |
| 8 | Relaxation curves | G(t), K(t) identified vs ground truth |
| 9 | Prony spectrum | Bar chart: identified vs ground truth terms |

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'core'`:** Run scripts from their own directory:
```bash
cd Forward_solver    # for forward scripts
cd inverse_problem   # for inverse scripts
```

**`FileNotFoundError: standard_FE_dataset/<N>`:** Generate data first, then set `EXPERIMENT_NUMBER` to match:
```bash
cd Forward_solver && python run_simulation.py
```

**Too few Prony terms after clustering:** Reduce `CLUSTERING_RANGE` from 0.3 to 0.2, or increase `N_MAXWELL_SHEAR` to 200.

**Large displacement residuals (> 10%):** Switch to `BOUNDARY_CONDITION = 'BottomForceBC'`, or expand `TAU_MAX`.

---

## Data Availability and Reproducibility

| Data | Reproducible from this repo | Notes |
|------|----------------------------|-------|
| Synthetic datasets (any geometry) | **Yes** | Run `Forward_solver/run_simulation.py` with any `experiment_id` |
| Ground-truth Prony parameters | **Yes** | Defined in `Forward_solver/core/material.py` |
| Plain PA6 creep specimen (DIC preprocessed) | No | Fraunhofer IWM (confidential) |
| Elliptical-hole PA6 specimen (DIC preprocessed) | No | Fraunhofer IWM (confidential) |
| Three-hole PA6 specimen (DIC preprocessed) | No | Fraunhofer IWM (confidential) |

Raw DIC displacement fields were measured at Fraunhofer Institut für Werkstoffmechanik (IWM), Freiburg, and are excluded from this public repository under a data confidentiality agreement. All synthetic-path results reported in the thesis are fully reproducible from this repository alone.

---

## Citation

```bibtex
@mastersthesis{undamatla2026euclid,
  title  = {Numerical Implementation and Experimental Validation of the
            EUCLID Framework for Linear Viscoelastic Parameter Identification},
  author = {Undamatla, Kali Satya Sri Charan},
  school = {Technische Universit\"{a}t Braunschweig,
            Institute of Structural Analysis (ISD)},
  year   = {2026},
  month  = {March},
  note   = {Supervised by Chang Kuo-I, Fraunhofer IWM. Code available at
            https://github.com/kaliundamatla/EUCLID-Linear-Viscoelsticity-2d}
}
```

---

**Supervisor:** Chang Kuo-I Msc, Fraunhofer Institute for Mechanics of Materials (IWM), Freiburg
**Examiner:** Prof. Dr.-Ing. Ursula Kowalsky, ISD, Technische Universität Braunschweig
**Python 3.8+ | numpy | scipy | matplotlib | scikit-learn | pygmsh | gmsh**
