#   electronicStructure_SCAN.py
#
#   This Python program loads a pre-relaxed structure CIF for LK-99, then sets up a 
#   Quantum ESPRESSO (QE) input file for a self-consistent field calculation, using 
#   the SCAN functional. After this, it runs the denisty functional theory calculat-
#   ion (pw.x), including a dense 6x6x6 k-point grid for Brillouin zone sampling. 
#   The SCAN fnctional (Meta-GGA) incorporates kinetic energy density, so it is 
#   better for band gaps, strongly correlated systems, and non-covalent interactions.
#   The results of the propgram, the SCF output file (with total energy, forces, 
#   eigenvalues, etc.), one can determine whether LK-99 is a metal, semiconductor, 
#   or a superconductor. 
#____________________________________________________________________________________

from pymatgen.io.quantumespresso import PWInput

# loads relaxed structure from PBEsol

structure_relaxed = Structure.from_file("relax_pbesol.out")

# generates input for SCAN calculation

input_scan = PWInput

(
    structure_relaxed, pseudopotentials = 
    
    {
        "Pb": "Pb.scan-rrkjus_psl.1.0.0.UPF",
        "Cu": "Cu.scan-rrkjus_psl.1.0.0.UPF",
        "P": "P.scan-rrkjus_psl.1.0.0.UPF",
        "O": "O.scan-rrkjus_psl.1.0.0.UPF",
    },

    control = {"calculation": "scf", "tstress": True},
    system = 
    
    {
        "ecutwfc": 70,  # a higher cutoff for SCAN
        "ecutrho": 560,
        "input_dft": "SCAN",
        "occupations": "smearing",
        "smearing": "mv",
        "degauss": 0.01,
    },

    kpts = (6, 6, 6),  # denser k-grid (for accuracy)
)

input_scan.write_file("scan_scf.in")

# runs Quantum ESPRESSO with SCAN

!mpirun -np 4 pw.x -in scan_scf.in > scan_scf.out

# verifies completion & downloads output

if "JOB COMPLETED" in open("scan_scf.out").read():

    print("SCAN calculation worked!")
    display(FileLink("scan_scf.out")) 

else:

    print("SCAN calculation failed. Check scan_scf.out.")