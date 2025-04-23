#   phonon-mediatedSuperconductivity_EPW_SCAN.py
#
#   Using Quantum ESPRESSO, the self-consistent field calculation creates input for 
#   and runs the method, pw.x, using the SCAN functional for the sake of computing 
#   the ground-state electron density. The total energy is parsed from the Fermi 
#   energy, which will be used later in EPW. Then, ph.x is used to compute a q-point
#   grid with dynamical matrices that is vital for electron-phonon coupling, 
#   including phonon frequencies and modes. Finally, the epw.x method finds the 
#   Eliashberg spectral function (α²F(ω)), the coupling strength, and the averaged
#   frequency. These outputs reveal which phonon modes (Pb/Cu vibrations) mediate 
#   pairing, the total EPC strength, and the soft modes that may enhance the resultant
#   critical temperature, if the coupling strength is large, respectively. All will 
#   be useful for determining the superconductivity of LK-99. 
#
#____________________________________________________________________________________


import os
from pymatgen.core import Structure
from pymatgen.io.quantumespresso import PWInput, PhononInput, EpwInput
from IPython.display import FileLink, display
import matplotlib.pyplot as plt
import numpy as np

# loads the relaxed structure

try:

    relaxed_struct = Structure.from_file("scan_relax.out") 
    print("Successfully loaded relaxed structure")

except:

    print("Error: Could not load relaxed structure. Please run SCAN relaxation first.")
    raise

# generates SCAN SCF input for accurate charge density

input_scf = PWInput

(
    relaxed_struct, pseudopotentials = 
    
    {
        "Pb": "Pb.scan-rrkjus_psl.1.0.0.UPF",
        "Cu": "Cu.scan-rrkjus_psl.1.0.0.UPF",
        "P": "P.scan-rrkjus_psl.1.0.0.UPF", 
        "O": "O.scan-rrkjus_psl.1.0.0.UPF",
    },

    control = {"calculation": "scf", "tstress": True, "prefix": "pwscf"},

    system = 
    
    {
        "ecutwfc": 70,
        "ecutrho": 560,
        "input_dft": "SCAN",
        "occupations": "smearing",
        "smearing": "mv",
        "degauss": 0.01,
    },

    kpts=(6, 6, 6),  # dense k-grid for accurate Fermi surface
)

input_scf.write_file("scan_scf.in")

# runs SCAN SCF calculation

print("Running SCAN SCF calculation...")
!mpirun -np 4 pw.x -in scan_scf.in > scan_scf.out

if "JOB COMPLETED" not in open("scan_scf.out").read():

    print("SCF calculation failed!")
    display(FileLink("scan_scf.out"))
    raise RuntimeError("SCF calculation failed")

# generates phonon calculation input

ph_input = PhononInput

(
    relaxed_struct, inputph = 
    
    {
        "tr2_ph": 1e-12,
        "prefix": "pwscf",
        "fildyn": "ph_dyn",
        "ldisp": True,
        "nq1": 4, "nq2": 4, "nq3": 4,  # q-point grid
        "epsil": True,  # For polar materials
    }
)

ph_input.write_file("ph.in")

# runs phonon calculation

print("Running phonon calculation...")
!ph.x -in ph.in > ph.out

if "JOB DONE" not in open("ph.out").read():
    
    print("Phonon calculation failed!")
    display(FileLink("ph.out"))
    raise RuntimeError("Phonon calculation failed")

# prepares EPW input for electron-phonon coupling

epw_input = f"""&inputepw
  prefix = 'pwscf'
  outdir = './'
  amass(1) = 207.2    ! Pb
  amass(2) = 63.55     ! Cu
  amass(3) = 30.97     ! P
  amass(4) = 16.00     ! O
  
  elph = .true.
  ep_coupling = .true.
  epwwrite = .true.
  epbread = .false.
  wannierize = .true.
  nbndsub = 24         ! Bands near Fermi level
  bands_skipped = 'exclude_bands = 1:8' ! Exclude core states
  
  ! SCAN-specific settings
  input_dft = 'SCAN'
  scdm_entanglement = 'isolated'
  scdm_mu = 0.0
  scdm_sigma = 1.0
  
  ! Fine q/k meshes
  nkf1 = 16, nkf2 = 16, nkf3 = 16
  nqf1 = 8, nqf2 = 8, nqf3 = 8
/
"""

with open("epw.in", "w") as f:

    f.write(epw_input)

# runs EPW calculation

print("Running EPW calculation...")
!mpirun -np 4 epw.x -in epw.in > epw.out

if "JOB COMPLETED" not in open("epw.out").read():

    print("EPW calculation failed!")
    display(FileLink("epw.out"))
    raise RuntimeError("EPW calculation failed")

# post-processing: extracts & plots Eliashberg function

print("Extracting results...")

try:
    # parses EPW output for Eliashberg function

    with open("epw.out") as f:
        lines = f.readlines()
    
    # finds Eliashberg function data

    start_idx = None
    for i, line in enumerate(lines):
        if "Eliashberg Spectral Function" in line:
            start_idx = i + 3
            break
    
    if start_idx is None:
        raise ValueError("Could not find Eliashberg function in output")
    
    omega = []
    alpha2F = []
    for line in lines[start_idx:]:
        if not line.strip():
            break
        parts = line.split()
        omega.append(float(parts[0]))
        alpha2F.append(float(parts[1]))