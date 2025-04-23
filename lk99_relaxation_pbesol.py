#   structuralRelaxation_PBEsol.py
#
#   The program loads the LK-99 crystallographic information file and reads it to
#   determine the initial atomic positions and lattice parameters for relation. 
#   Then, Quantum ESPRESSO-friendly input is generated for the relaxation process, 
#   where a 4x4x4 k-point grid is created for Brillouin zone sampling and the 
#   Methfessel-Paxton smearing method is applied to minimize artificial entropy 
#   effects. Running pw.x (from the QE main DFT engine) relaxes atomic positions 
#   and cell volume. This convergence calculation is verified and the result is 
#   displayed to the user, followed by a final relaxed structure, in an interactive 
#   3D viewer.
#
#____________________________________________________________________________________


from pymatgen.core import Structure 
from pymatgen.io.quantumespresso import PWInput 
from pymatgen.visualization import plotly_vis
from pymatgen.io.quantumespresso import QEOutput

# loads the LK-99 CIF 

structure = Structure.from_file("LK99.cif") 

# generates QE input for PBEsol relaxation 

input_pbesol = PWInput 
    
(structure, pseudopotentials = 
   
    { 
        "Pb": "Pb.pbesol-n-rrkjus_psl.1.0.0.UPF", 
        "Cu": "Cu.pbesol-n-rrkjus_psl.1.0.0.UPF", 
        "P": "P.pbesol-n-rrkjus_psl.1.0.0.UPF", 
        "O": "O.pbesol-n-rrkjus_psl.1.0.0.UPF", 
    }, 
    
control = {"calculation": "relax", "tprnfor": True, "tstress": True}, 
    
system = 
    
    { 
        "ecutwfc": 60, 
        "ecutrho": 480, 
        "input_dft": "PBEsol",  # setting PBEsol 
        "occupations": "smearing", 
        "smearing": "mv", 
        "degauss": 0.02, 
    }, 
    
kpts = (4, 4, 4),) 
    
input_pbesol.write_file("relax_pbesol.in") 

# runs Quantum ESPRESSO (where pw.x is in the PATH)
    
!mpirun -np 4 pw.x -in relax_pbesol.in > relax_pbesol.out

# verifies completion
    
if "JOB COMPLETE" in open("relax_pbesol.out").read():
    
    print("Relaxation worked!")
    
# displays link to download output
    
display(FileLink("relax_pbesol.out"))

else:
print("Relaxation failed. Check relax_pbesol.out.")

# visualizes strcture 

qe_out = QEOutput("relax_pbesol.out")

final_structure = qe_out.final_structure
    
plotly_vis.show_structure(final_structure)
