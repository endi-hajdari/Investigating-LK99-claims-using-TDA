#   topologicalInvariantComputation_Z2Pack.py
#
#   Written in Python, this program first converts CIF LK-99 data into a format 
#   suitable for Quantum ESPRESSO (QE), which can be used with SCAN+U for a ground-
#   state calculation. Then, it runs Z2Pack to compute the topological invariant,
#   ℤ₂. The calculation carried out relies on time-reversal symmetry, where spin-
#   orbit coupling is important. The Wilson loop method of Z2Pack computes the Berry 
#   phase around 2D planes in the Brillouin zone. This method tracks the evolution of
#   Wannier centres. After the setting have been calibrated, the surface calculation
#   for ℤ₂ is carried out. The result, a file which contains wavefunction data, is 
#   finally saved. If LK-99 is a topological insulator, with ℤ₂ = 1, then doping or 
#   certain effects could cause superconductivity.
#____________________________________________________________________________________

import os
import z2pack
from pymatgen.core import Structure
from shutil import which

# converts CIF to Quantum ESPRESSO input

def cif_to_pw(cif_path, output_dir="qe_input"):

    os.makedirs(output_dir, exist_ok=True)

    struct = Structure.from_file(cif_path)
    
    # writes an input file

    with open(f"{output_dir}/lk99.in", "w") as f:

        f.write
        
    (f"""&CONTROL
    calculation = 'scf'
    prefix = 'lk99'
    outdir = './out'
    pseudo_dir = './pseudos'
    tprnfor = .true.
    tstress = .true.
/
&SYSTEM
    ibrav = 0
    nat = {len(struct)}
    ntyp = 4  # Pb, Cu, P, O
    ecutwfc = 60
    ecutrho = 480
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.01
    nspin = 2
    lspinorb = .true.  # Include SOC for ℤ₂
/
&ELECTRONS
    conv_thr = 1e-8
    mixing_beta = 0.7
/
ATOMIC_SPECIES
Pb 207.2 Pb.upf
Cu 63.55 Cu.upf
P  30.97  P.upf
O  16.00  O.upf
K_POINTS automatic
6 6 4 0 0 0
""")
    
        f.write("CELL_PARAMETERS angstrom\n")
    
        for vec in struct.lattice.matrix:

            f.write(f"{vec[0]:.8f} {vec[1]:.8f} {vec[2]:.8f}\n")

        f.write("ATOMIC_POSITIONS angstrom\n")
    
        for site in struct:

            f.write(f"{site.species_string} {site.x:.8f} {site.y:.8f} {site.z:.8f}\n")

# runs Z2Pack to compute ℤ₂

def compute_z2_invariant():

    # defines QE command (paths are adjustable)

    qe_command = "mpirun -np 4 pw.x -in lk99.in > lk99.out"
    
    # calibrates Z2Pack settings

    settings = 
    
    {
        "qe_executable": which("pw.x"),

        "num_lines": 8,  

        "pos_tol": 1e-2,

        "gap_tol": 0.3,

        "move_tol": 0.3,
    }
    
    # runs surface calculation for ℤ₂

    result = z2pack.surface.run
    
    (

        system=z2pack.fp.System
        (
            input_files=["lk99.in"],

            kpt_fct=z2pack.fp.kpts.qe_explicit,

            kpt_path="lk99.in",

            command=qe_command,

            executable=settings["qe_executable"],
        ),

        surface=lambda s, t: [s / 2, t],  

        **settings
    )
    
    # saves the result

    z2pack.save(result, "z2_results.json")

    print(f"ℤ₂ invariant is {z2pack.invariant.z2(result)}.")

if __name__ == "__main__":

    cif_path = "lk99_relaxed.cif"  # can be edited

    cif_to_pw(cif_path)

    compute_z2_invariant()
