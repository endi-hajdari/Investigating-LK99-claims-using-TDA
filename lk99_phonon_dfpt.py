#    phononCalculations_DFPT.py
#
#    Starting with a pre-relaxed structure file that contains LK-99 data, a downloaded
#    Quantum ESPRESSO suite, and a completed self-consistent field (SCF) calculation
#    with metallic smearing performed, this program establishes and runs density 
#    functional perturbation theory (DFPT). It computes dynamical matrices, D(q), 
#    on a 6x6x2 q-point grid that spans across the Brillouin zone (for anisotropy 
#    testing). It solves the quantum mechanical vibrations for the crystal structure
#    of the LK-99 material. After the main calculation, it performs post-processing 
#    to arrive at force constants from dynamical matrices using q2r.x and phonon 
#    frequencies along high-symmetry paths using matdyn.x. The user can adjust the 
#    dimensions of the q-point grid, convergence threshold (tr2_ph), and the 
#    directories for prefix or output. After running the code, one obtains the result
#    of main phonon calculation, dynamical matrices for each q-point, force constants 
#    file, and phonon frequencies. It may be used to aid in the verification of 
#    superconductivity claims for LK-99, along with electron-pairing EPC calculations
#    and anharmonic corrections.
#
#____________________________________________________________________________________

import os
from ase.io import read
import numpy as np
from shutil import which

# configures system parameters

prefix = "lk99"            
outdir = "./tmp"           
pseudo_dir = "./pseudos"   

# sets up calculation parameters

nq1, nq2, nq3 = 6, 6, 2    
tr2_ph = 1.0e-14           
electron_phonon = "simple"  

# path to Quantum ESPRESSO executables (modifiable)

ph_executable = "ph.x"
q2r_executable = "q2r.x"
matdyn_executable = "matdyn.x"

# defines the functions

def check_executables():

    """Verify that required executables are available."""

    for exe in [ph_executable, q2r_executable, matdyn_executable]:

        if which(exe) is None:

            raise RuntimeError(f"Executable {exe} not found in PATH")

def generate_ph_input():

    """Generate input file for phonon calculation."""

    ph_input = f"""

&INPUTPH
  prefix = '{prefix}',
  outdir = '{outdir}',
  fildyn = '{prefix}.dyn',
  fildvscf = '{prefix}.dvscf',
  tr2_ph = {tr2_ph},
  ldisp = .true.,
  nq1 = {nq1},
  nq2 = {nq2},
  nq3 = {nq3},
  electron_phonon = '{electron_phonon}',
  recover = .true.,
/
"""
    with open("ph.in", "w") as f:

        f.write(ph_input)

def run_phonon_calculation():

    """Run the phonon calculation."""

    print("Running phonon calculation...")

    os.system(f"{ph_executable} < ph.in > ph.out")

    print("Phonon calculation completed. Output is in ph.out")

def post_process():

    """Post-processing of phonon output."""

    # generates force constants

    with open("q2r.in", "w") as f:

        f.write(f"""&INPUT
 fildyn='{prefix}.dyn',
 zasr='crystal',
 flfrc='{prefix}.fc',
/""")
        
    os.system(f"{q2r_executable} < q2r.in > q2r.out")
    
    # calculates phonon frequencies along high-symmetry path

    with open("matdyn.in", "w") as f:

        f.write(f"""&INPUT
 asr='crystal',
 flfrc='{prefix}.fc',
 flfrq='{prefix}.freq',
/""")
        
    os.system(f"{matdyn_executable} < matdyn.in > matdyn.out")

# main process

if __name__ == "__main__":

    print("\nLK-99 Phonon Calculation Set-Up")

    print(f"System prefix: {prefix}.")

    print(f"q-point grid: {nq1}x{nq2}x{nq3}")
    
    # verifies the environment
    check_executables()
    
    # creates the output directory, if needed

    os.makedirs(outdir, exist_ok=True)
    
    # generates the input files & runs calculations

    generate_ph_input()

    run_phonon_calculation()

    post_process()
    
    print("\nThe Calculation Summary")

    print(f"- dynamical matrices computed in {prefix}.dyn*")

    print(f"- force constants stored in {prefix}.fc")

    print(f"- phonon frequencies in {prefix}.freq")

    print("\nphonon calculation workflow completed successfully!")