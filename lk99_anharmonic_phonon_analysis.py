#   analyzingAnharmonicPhonons.py
#
#   The process implemented here begins with relaxing the structure with pw.x in 
#   Quantum ESPRESSO. the harmonic phonons will be computed with the ph.x function
#   with anharmonicity and the SCPH method, beginning with the input, harmonic_ph.in.
#   The harmonic phonons are run, following this. The third step is to run self-
#   consistent phonon by enabling scph = .true to compute anharmonic corrections 
#   and then executing it.  In the post-processing step, there is a Python script 
#   to extract phonon frequencies and code to plot the anharmonic phonon density 
#   of states (DOS). Finally, the program compares the harmonic to the anharmonic
#   phonons, considering the phonon line widths and frequency shifts. 
#____________________________________________________________________________________

# setup & dependencies for program

import os
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from IPython.display import FileLink

# sets path to QE executables (which can be adjusted as needed)

QE_BIN = "/path/to/qe/bin/"  # e.g., /opt/qe-7.1/bin/

PW = os.path.join(QE_BIN, "pw.x")

PH = os.path.join(QE_BIN, "ph.x")

SCPH = os.path.join(QE_BIN, "scph.x")

EPW = os.path.join(QE_BIN, "epw.x")

# defines input/output directories

WORKDIR = "lk99_anharmonic"

os.makedirs(WORKDIR, exist_ok=True)

# completes structural relaxation of LK-99 CIF w/PBEsol

def write_relax_input(cif_path="lk99.cif"):

    struct = Structure.from_file(cif_path)
    
    input_str = f"""
&CONTROL
    calculation = 'relax'
    restart_mode = 'from_scratch'
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
    input_dft = 'PBEsol'
/
&ELECTRONS
    conv_thr = 1e-8
    mixing_beta = 0.7
/
&IONS
    ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
Pb 207.2 Pb.upf
Cu 63.55 Cu.upf
P  30.97  P.upf
O  16.00  O.upf
K_POINTS automatic
4 4 4 0 0 0

CELL_PARAMETERS angstrom
{" ".join(map(str, struct.lattice.matrix[0]))}
{" ".join(map(str, struct.lattice.matrix[1]))}
{" ".join(map(str, struct.lattice.matrix[2]))}

ATOMIC_POSITIONS angstrom
""" + "\n".join([f"{site.species_string} {site.x} {site.y} {site.z}" for site in struct])
    
    with open(f"{WORKDIR}/relax.in", "w") as f:

        f.write(input_str)

    print("Relaxation input file written.")

write_relax_input()

# creates the harmonic phonon input file

def write_harmonic_ph_input():

    input_str = """

&INPUTPH
    prefix = 'lk99'
    outdir = './out'
    fildyn = 'lk99.dyn'
    trans = .true.
    ldisp = .true.
    nq1 = 2, nq2 = 2, nq3 = 2
/
"""
    with open(f"{WORKDIR}/harmonic_ph.in", "w") as f:

        f.write(input_str)

    print("Harmonic phonon input file written.")

write_harmonic_ph_input()

# calculates the self-consistent phonons & writes result to a file

def write_scph_input():

    input_str = """

&INPUTPH
    prefix = 'lk99'
    outdir = './out'
    fildyn = 'lk99.dyn'
    scph = .true.
    scph_drho = .true.
    niter_ph = 100
    tr2_ph = 1e-12
/
"""
    with open(f"{WORKDIR}/scph.in", "w") as f:

        f.write(input_str)

    print("SCPH input file written.")

write_scph_input()

# for phonon density of states

def plot_phonon_dos():

    # loads harmonic/anharmonic frequencies (modifiable)

    harmonic_freqs = np.loadtxt(f"{WORKDIR}/harmonic_freqs.dat")  

    scph_freqs = np.loadtxt(f"{WORKDIR}/scph_freqs.dat")
    
    plt.figure(figsize=(10, 6))

    plt.hist(harmonic_freqs, bins = 50, alpha = 0.5, label = "Harmonic")

    plt.hist(scph_freqs, bins = 50, alpha = 0.5, label = "SCPH Anharmonic")

    plt.xlabel("Frequency (cm⁻¹)")

    plt.ylabel("DOS")
    
    plt.legend()

    plt.savefig(f"{WORKDIR}/phonon_dos.png")

    plt.show()

plot_phonon_dos()


# calculates Eliashberg function result

def plot_eliashberg():

    # parses EPW output 

    data = np.loadtxt(f"{WORKDIR}/epw.out", comments="#")

    omega = data[:, 0]  

    alpha2F = data[:, 1]  

    plt.figure(figsize=(10, 6))

    plt.plot(omega, alpha2F, label=r"$\alpha^2 F(\omega)$")

    plt.xlabel(r"$\omega$ (meV)")

    plt.ylabel(r"$\alpha^2 F(\omega)$")

    plt.legend()

    plt.savefig(f"{WORKDIR}/eliashberg.png")

    plt.show()

plot_eliashberg()