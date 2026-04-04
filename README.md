
# Quantitative CSF–ISF Exchange Modeling from MRI

This repository contains code and data for estimating cerebrospinal fluid (CSF) to interstitial fluid (ISF) exchange in humans using quantitative MRI following intrathecal (IT) and intravenous (IV) contrast administration.

Accompanies:
"Quantitative assessment of flow between cerebrospinal and interstitial fluid compartments in humans" , Wåhlin et al. (2026), *PNAS*

---

## Overview

We implement a compartmental model to estimate:

- CSF → ISF inflow (IT-experiments) or ISF → CSF outflow (IV-experiments)
- ISF volume fraction (both experiemnts)

from time-resolved MRI-derived concentration curves.

Two complementary paradigms:

- Intrathecal (IT): CSF → tissue
- Intravenous (IV): tissue → CSF

---

## Repository Structure

- `ITexperiments.m` — Intrathecal analysis
- `IVexperiments.m` — Intravenous analysis
- `itgflow.m` — Simulation / validation
- `dataset1.mat` — IT dataset
- `dataset2.mat` — IV dataset
- `simulatedData.mat` — Synthetic data

---

## Quick Start

### Requirements
- MATLAB
- Optimization Toolbox (`fminunc`, `fmincon`)

### Run intrathecal analysis
```matlab
run ITexperiments.m
```

### Run intravenous analysis
```matlab
run IVexperiments.m
```

### Run simulation
```matlab
run itgflow.m
```

---

## IT-Model

Exchange model:

- Parameters:
  - k1
  - k2

Derived:
- q = k1 × tissue volume
- ve = k1 / k2

---

## IV-Model

Exchange model:

- Parameters:
  - k1
  - k2

Derived:
- q = k1 × CSF volume
- ve = k2 (model-specific estimate)

Note: The IV formulation uses a different parameterization than the IT model.

---

## Data

### IT (dataset1.mat)
- CSF and tissue concentrations per region
- Tissue volumes
- Time vectors

### IV (dataset2.mat)
- CSF and tissue concentrations
- CSF volume
- Time vectors

### Simulation
- Known ground truth (q, ve)

---

## Notes

- Estimates depend on model assumptions
- Optimization may be sensitive to initialization/version/optimization parameters
- Low SNR may produce unstable estimates (mostly the case for IV)

---

## Citation

if you use the code, please cite the paper "Quantitative assessment of flow between cerebrospinal and interstitial fluid compartments in humans" by Wåhlin et al. *PNAS*.

---

## Contact

Anders Wåhlin  
Umeå University  
anders.wahlin@umu.se
