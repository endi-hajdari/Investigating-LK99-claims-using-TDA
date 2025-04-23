#   SCF_highPrecision.py
#
#   This program performs a high-precision self-consistent field (SCF) calculation 
#   on LK-99 data using Quantum ESPRESSO. Its input is a pre-relaxed CIF file of the
#   material and it produces a converged charge density for upcoming DFPT/phonon
#   calculations. The settings feature incredibly tight convergence, Marzari-Vander-
#   bilt smearing for metallic systems, and hybrid functional support for improved 
#   gaps. The script is optimized for the LK-99 material, with a dense k-point grid 
#   that makes use of Fermi surface features and stabilizing advanced mixing. It 
#   handles errors by checking SCF convergence in a separate output file and if this 
#   fails, an exception is raised. The results will, ultimately, be used to estimate 
#   the superconducting transition temperature. 
#____________________________________________________________________________________


import os
from ase.io import read, write
from ase.calculators.espresso import Espresso

# configures paths & file names

cif_file = "lk99_relaxed.cif"  

output_dir = "./high_precision_scf/"

pseudo_dir = "./pseudos/"  # a directory w/pseudopotentials

# sets Quantum ESPRESSO parameters (with high-precision settings)

input_data = {
    "control": {
        "calculation": "scf",
        "prefix": "lk99",
        "outdir": output_dir,
        "pseudo_dir": pseudo_dir,
        "tprnfor": True,
        "tstress": True,
        "verbosity": "high",
        "disk_io": "low",  
    },
    "system": {
        "ecutwfc": 80, 
        "ecutrho": 640,
        "occupations": "smearing",
        "smearing": "marzari-vanderbilt",  
        "degauss": 0.02,  
        "nspin": 1, 
        "input_dft": "HSE", 
        "exx_fraction": 0.25,  
    },
    "electrons": {
        "electron_maxstep": 200, 
        "scf_must_converge": True,
        "conv_thr": 1e-10,  
        "mixing_mode": "local-TF", 
        "mixing_beta": 0.3, 
        "diagonalization": "david",  
    },
}

# pseudopotentials that can be adjusted for individual systems

pseudopotentials = 

{
    "Pb": "Pb.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "Cu": "Cu.pbe-n-kjpaw_psl.1.0.0.UPF",
    "P": "P.pbe-n-kjpaw_psl.1.0.0.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
}

# k-point grid that is dense for metals

kpts = (8, 8, 8)  

def run_high_precision_scf():

    # creates output directory
    os.makedirs(output_dir, exist_ok=True)

    # loads pre-relaxed CIF

    atoms = read(cif_file)

    # setups up Quantum ESPRESSO calculator
    calculator = Espresso
    (
        input_data=input_data,
        pseudopotentials=pseudopotentials,
        kpts=kpts,
        directory=output_dir,
    )

    atoms.calc = calc

    # runs SCF calculation & finds potential energy

    print("Starting high-precision SCF.")

    energy = atoms.get_potential_energy() 

    print(f'SCF calculation completed. The total energy is {energy} eV.')

    # checks convergence

    with open(os.path.join(output_dir, "lk99.pwo"), "r") as f:

        log = f.read()

        if "convergence HAS BEEN ACHIEVED" in log:

            print("SCF converged successfully!")

        else:

            raise RuntimeError("SCF failed to converge.")

    # saves results

    write(os.path.join(output_dir, "lk99_scf.pwi"), atoms) 

    print(f"Results have been saved to {output_dir}.")

if __name__ == "__main__":

    run_high_precision_scf()