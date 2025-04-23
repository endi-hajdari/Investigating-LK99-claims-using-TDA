#   smearing_MP.py
#
#   First, the program loads the pre-relaxed LK-99 structure from a labeled CIF. 
#   Then, the input is prepared and a self-consistent field calculation is set up. 
#   The Methfessel-Paxton smearing is configured so that there is a broadening of the
#   electronic states near the Fermi level with degauss = 0.02, which is required by 
#   the (disputed) metallicity of LK-99 to stabilize the calculation. For the self-
#   consistent field component, the Kohn-Sham equations are solved to arrive at the 
#   total energy, Fermi energy, and charge density. The output consists of the former
#   two values, written into a new file for possible further analysis. The charge 
#   density can be used to identify if Cu doping has introduced charge inhomogeneity
#   in the LK-99 material. 
#
#____________________________________________________________________________________

from pymatgen.core import Structure
from pymatgen.io.quantumespresso import PWInput
from IPython.display import FileLink, display
import os

# loads the pre-relaxed LK-99 structure

try:
    relaxed_struct = Structure.from_file("relaxed_LK99.cif") 

 # from QE output

    print("Successfully loaded relaxed LK-99 structure")

except FileNotFoundError:

    print("Error: Could not find relaxed structure file.")

    print("Please provide the relaxed structure as 'relaxed_LK99.cif' or similar.")

    exit()

# sets up Quantum ESPRESSO input w/Methfessel-Paxton smearing

input_params = 
{
    "structure": relaxed_struct,
    "pseudopotentials": 

    {"Pb": "Pb.pbesol-n-rrkjus_psl.1.0.0.UPF",
     "Cu": "Cu.pbesol-n-rrkjus_psl.1.0.0.UPF",
     "P": "P.pbesol-n-rrkjus_psl.1.0.0.UPF",
     "O": "O.pbesol-n-rrkjus_psl.1.0.0.UPF",},

    "control": 

    {"calculation": "scf", 
     "restart_mode": "from_scratch",
     "prefix": "lk99",
     "outdir": "./tmp",
     "tstress": True,
     "tprnfor": True,},

    "system": 

    {"ecutwfc": 60,        
     "ecutrho": 480,      
     "input_dft": "PBEsol",
     "occupations": "smearing",
     "smearing": "mv",     
     "degauss": 0.02,      
     "nspin": 1,},         

    "electrons": 

    {"conv_thr": 1.0e-8,   
     "mixing_beta": 0.7,   
     "electron_maxstep": 100,},

    "kpts": (6, 6, 6)}    
# generates the input file

qe_input = PWInput(**input_params)

qe_input.write_file("lk99_scf.in")

# runs Quantum ESPRESSO calculation

print("Running electronic structure calculation with Methfessel-Paxton smearing...")

os.system("mpirun -np 4 pw.x -in lk99_scf.in > lk99_scf.out")

# checks results & analyzes output

if "JOB DONE" in open("lk99_scf.out").read():
    print("Calculation completed successfully!")
    
    # parses the output

    from pymatgen.io.quantumespresso import QEOutput
    qe_out = QEOutput("lk99_scf.out")
    
    print("\nKey results:")
    print(f"Fermi energy: {qe_out.fermi:.4f} eV")
    print(f"Total energy: {qe_out.final_energy:.4f} eV")
    
    # links to output file

    display(FileLink("lk99_scf.out")
