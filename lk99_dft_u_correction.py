#   Hubbard+UCorrection.py
#
#   The program below applies density functional theory (DFT) calculations with 
#   Hubbard U corrections (DFT+U) on a pre-relaxed LK-99 structure, so that the 
#   electronic properties of the material can be studied thoroughly. This code 
#   is particularly useful in evaluating the superconductivity of the LK-99 
#   material. The program loads the relaxed LK-99 structure, reading the atomic 
#   position and lattice parameters from a particular output file that had been 
#   previously generated in a prior structural relaxation with the PBEsol 
#   functional. 
#____________________________________________________________________________________

    import os
    from pymatgen.core import Structure
    from pymatgen.io.quantumespresso import PWInput
    from IPython.display import FileLink, display

    # loads the relaxed structure

    relaxed_structure = Structure.from_file("relax_pbesol.out")

    # defines the Hubbard U parameters (Cu-3d)

    hubbard_u = {"Cu": {"L": 2, "U": 5.0, "J": 0.0, "Hubbard_U": 5.0}}

    # generates QE input with +U correction

    input_hubbard = PWInput

    (
        
        relaxed_struct, pseudopotentials = 
        {
            "Pb": "Pb.pbesol-n-rrkjus_psl.1.0.0.UPF",
            "Cu": "Cu.pbesol-n-rrkjus_psl.1.0.0.UPF",
            "P": "P.pbesol-n-rrkjus_psl.1.0.0.UPF",
            "O": "O.pbesol-n-rrkjus_psl.1.0.0.UPF",
        },

        control = {"calculation": "scf", "tstress": True, "prefix": "pwscf"},

        system = 
        
        {
            "ecutwfc": 60,
            "ecutrho": 480,
            "occupations": "smearing",
            "smearing": "mv",
            "degauss": 0.02,
            "lda_plus_u": True,  # Enable +U
            "hubbard_u": hubbard_u["Cu"]["Hubbard_U"],
            "hubbard_j0": hubbard_u["Cu"]["J"],
        },

        kpts=(4, 4, 4),
    
    )

    input_hubbard.write_file("hubbard_u.in")

    # runs QE with +U

    !mpirun -np 4 pw.x -in hubbard_u.in > hubbard_u.out

    # verifies completion & analyzeS results

    if "JOB COMPLETED" in open("hubbard_u.out").read():

        print("+U calculation succeeded!")
        display(FileLink("hubbard_u.out"))  
    
    # parses total energy

    from pymatgen.io.quantumespresso import QEOutput

    qe_out = QEOutput("hubbard_u.out")

    print(f"Total energy with +U: {qe_out.final_energy} eV")

    else:

        print("+U calculation failed. Check hubbard_u.out.")