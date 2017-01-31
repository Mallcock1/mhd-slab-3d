# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 15:15:13 2016

@author: Matt
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:19:48 2016

@author: Matt
"""


############################
#asymmetric slab EXACT dispersion plots
#Slow sausage and kink modes
#K fixed
##########################

import slab_functions as sf
import solver
import numpy as np
import scipy as sc
from scipy.optimize import newton
import matplotlib.pyplot as plt

##########################


#fix R1
R1 = 1.5 ##rho_1 / rho_0
c1 = sf.c2 * sc.sqrt(sf.R2 / R1)

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)

#number of iterations
N = 1000
#Nb = 100

# Starting W values
W_init_vals = [sf.vA - 0.00000001, 0.3783] #[0.01, sf.cT - 0.01]
W_fin_vals = [sf.cT, sf.cT] #[sf.cT, sf.cT]
W_ntries = [500] * len(W_init_vals)

#Define zero vectors, which will be filled in the for loop
Kvals = np.zeros(N)
Wvals_list = []

#Kvals_backwards = np.zeros(Nb)
#Wvals_backwards_list = []

#initial K and W values
Kmin = 2.
Kmax = 5.

for i in range(len(W_init_vals)):
    W = W_init_vals[i]
    K = Kmin
    Wvals_list.append(np.zeros(N, dtype=complex))
    #for each K value, find the root of the dispersion function
    for j in range(0,N):
        if j == 0:
            Wvals_list[i][j] = solver.solver_backwards(disp_rel_asym_1var, W_init_vals[i], W_fin_vals[i], W_ntries[i])
        else:
            Wvals_list[i][j] = newton(disp_rel_asym_1var, W + 0.0000000001, tol=1e-5, maxiter=50)
        Kvals[j] = K
        if j == 0:
            W = Wvals_list[i][j]
        else:
            W = 2 * Wvals_list[i][j] - Wvals_list[i][j-1]
        K = K + (Kmax - Kmin) / N

##and backwards
#for i in range(len(W_init_vals)):
#    W = W_init_vals[i]
#    K = Kmin
#    Wvals_backwards_list.append(np.zeros(Nb, dtype=complex))
#    #for each K value, find the root of the dispersion function
#    for j in range(0,Nb):
#        if j == 0:
#            Wvals_backwards_list[i][j] = solver.solver_backwards(disp_rel_asym_1var, W_init_vals[i], W_fin_vals[i], W_ntries[i])
#        else:
#            Wvals_backwards_list[i][j] = newton(disp_rel_asym_1var, W, tol=1e-5, maxiter=50)
#        Kvals_backwards[j] = K
#        if j == 0:
#            W = Wvals_backwards_list[i][j]
#        else:
#            W = 2 * Wvals_backwards_list[i][j] - Wvals_backwards_list[i][j-1]
#        K = K - (Kmin) / Nb

###############
##############
############
##########
########

Kvals_full = np.linspace(0, Kmax)

fig = plt.figure(figsize=[5, 8])
#plot them both and label axes etc
plt.plot(Kvals_full, sf.cT * np.ones_like(Kvals_full), 'k-.')
plt.plot(Kvals_full, sf.vA * np.ones_like(Kvals_full), 'k-.')
plt.plot(Kvals_full, sf.c0 * np.ones_like(Kvals_full), 'k-.')
plt.plot(Kvals_full, sf.c2 * np.ones_like(Kvals_full), 'k-.')
plt.plot(Kvals_full, c1 * np.ones_like(Kvals_full), 'k-.')
K = 2.
plt.plot(K, solver.solver_backwards(disp_rel_asym_1var, sf.vA-0.000001, sf.cT, 100), 'o')
plt.plot(K, newton(disp_rel_asym_1var, 0.3783, tol=1e-5, maxiter=50), 'ro')

for i in  range(len(W_init_vals)):
    if solver.what_mode_s(Wvals_list[i][N/2], Kvals[N/2], R1)[0] == 'SAUSAGE':
        plt.plot(Kvals, Wvals_list[i], 'k-')
#        plt.plot(Kvals_backwards, Wvals_backwards_list[i], 'k-')
    elif solver.what_mode_k(Wvals_list[i][N/2], Kvals[N/2], R1)[0] == 'KINK':
        plt.plot(Kvals, Wvals_list[i], 'k--')
#        plt.plot(Kvals_backwards, Wvals_backwards_list[i], 'k--')
    else:
        print('ERROR: saus or kink?')
        
plt.xlabel(r'$kx_0$', fontsize='xx-large')
plt.ylabel(r'$\omega/kc_0$', fontsize='xx-large')
plt.ylim(0., sf.c0 + 0.05)
plt.annotate(r'$c_T$', xy=(0.97 * Kmax, sf.cT+0.005),
             xytext=(0.92 * Kmax, sf.cT + 0.005), fontsize='x-large')
plt.legend(loc=4)

plt.tight_layout()
plt.show()
#plt.savefig(u'D:/my_work/visualisations/xix_of_x_disp/disp.png')
