# LK-99: Investigating Superconductivity with DFT, EPW, and Topological Data Analysis (TDA)

This repository contains a complete computational pipeline to evaluate the superconductivity claims of LK-99 using Density Functional Theory (DFT), Density Functional Perturbation Theory (DFPT), Electron-Phonon coupling calculations (EPW), and Topological Data Analysis (TDA). The workflow integrates Quantum ESPRESSO, Z2Pack, and Giotto-TDA libraries, among others.

## Overview

In 2023, LK-99 was proposed as a potential room-temperature superconductor. This project critically evaluates that claim by computing:

- Electronic structure with DFT (PBEsol, SCAN, and SCAN+U)

- Phonon dispersion and electron-phonon coupling (DFPT, EPW)

- Eliashberg spectral function  and EPC strength 

- Topological invariants () using Z2Pack

- Structural topology from CIFs using Persistent Homology and TDA

## Dependencies

- Python 3.10+

- Quantum ESPRESSO (pw.x, ph.x, epw.x)

- Z2Pack

- Giotto-TDA

- pymatgen, matplotlib, numpy, ASE, IPython

## Purpose

This project aims to:

- Verify LK-99’s potential metallicity or superconductivity

- Compare its phonon-mediated pairing strength with known superconductors

- Identify structural/topological features that correlate with superconductivity

## License

MIT License
