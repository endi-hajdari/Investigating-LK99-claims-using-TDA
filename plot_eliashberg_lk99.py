#   plottingEliashberg.py
#
#   The program, written in Python, loads Eliashberg function data from a text file
#   that contains the frequency, ω, in meV, and α²F(ω) values. It plots the 
#   relationship between α²F(ω) and ω, then saves the result as .png file. In the 
#   plot, phonon frequency lies on the x-axis and electron-phonon coupling strength 
#   is on the y-axis. The EPC strength indicates whether the material, LK-99 is a 
#   strong-coupling superconductor or if it is weak-coupling. The result, λ, is used 
#   in the McMillan-Allen-Dynes formula to calculate the critical computer. The peaks
#   Analzying the plot, peaks in α²F(ω) hint at strong electron-phonon coupling at 
#   certain frequencies. The area under the curve is approximately the EPC strength,
#   which, if large enough, means that LK-99 might be a superconductor. 
#____________________________________________________________________________________

import numpy as np
import matplotlib.pyplot as plt

# loads Eliashberg function data (from EPW)

omega, alpha2F = np.loadtxt("eliashberg.dat", unpack=True)

# α²F(ω) vs. ω

plt.plot(omega, alpha2F, 'b-', lw=2, label="SCAN + U")
plt.xlabel("Frequency (meV)")
plt.ylabel(r"α² F(ω)")
plt.title("Eliashberg Function of LK-99 (SCAN Meta-GGA)")
plt.legend()

# saves the plot as an image 

plt.savefig("eliashberg_scan.png", dpi=300)