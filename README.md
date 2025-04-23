# LK-99: Investigating Superconductivity with DFT, EPW, and Topological Data Analysis (TDA)

This repository contains a complete computational pipeline to evaluate the superconductivity claims of LK-99 using:

- Density Functional Theory (DFT)
- Density Functional Perturbation Theory (DFPT)
- Electron-Phonon Coupling (EPW)
- Topological Data Analysis (TDA)

The project integrates Quantum ESPRESSO, Z2Pack, and Giotto-TDA libraries to assess electronic, phononic, and topological properties of LK-99.

---

## Overview

In 2023, LK-99 was proposed as a potential room-temperature superconductor. This project evaluates that claim by computing:

- **Electronic structure** with DFT (PBEsol, SCAN, and SCAN+U)
- **Phonon dispersion & EPC** using DFPT and EPW
- **Eliashberg spectral function** \( \alpha^2 F(\omega) \) and EPC strength \( \lambda \)
- **Topological invariant** \( \mathbb{Z}_2 \) using Z2Pack
- **Atomic-scale topology** via Persistent Homology

## Dependencies

- Python 3.10+

- Quantum ESPRESSO (pw.x, ph.x, epw.x)

- Z2Pack

- Giotto-TDA

- pymatgen, numpy, matplotlib, ASE, IPython

## Project Goals
1. Verify LK-99’s potential metallicity or superconductivity

2. Compare its phonon-mediated pairing with known superconductors

3. Analyze its structural/topological features using TDA

## References

**Original claim:** arXiv:2307.12008

**Persistent homology in materials:**
Rucco et al. (2015), Chazal & Michel (2017)

**EPC theory:** Allen & Dynes (1975), Giustino et al. (2017)

## License
MIT License




