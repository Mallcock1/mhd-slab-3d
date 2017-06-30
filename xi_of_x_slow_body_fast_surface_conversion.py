# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:35:03 2017

@author: Matt
"""

import slab_functions as sf
import numpy as np
import scipy as sc
from scipy.optimize import newton
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib
import dispersion_diagram
import toolbox as tool

##########################

#*****************REQUIRES SBS SPEEDS IN slab_functions

show_xi_of_x = False
show_disp = False

show_xi_of_x = True
show_disp = True


mode_options = ['slow-kink-body-1']

alfven_mode_options = []
kink_mode_options = []
saus_mode_options = []
slow_surf_mode_options = []
fast_surf_mode_options = []
fast_kink_mode_options = []
fast_saus_mode_options = []
slow_body_1_mode_options = []
slow_body_2_mode_options = []
slow_body_3_mode_options = []
slow_body_mode_options = []
fast_body_1_mode_options = []
fast_body_2_mode_options = []
fast_body_3_mode_options = []

for mode in mode_options:
    if 'alfven' in mode:
        alfven_mode_options.append(mode)
    if 'kink' in mode:
        kink_mode_options.append(mode)
    if 'saus' in mode:
        saus_mode_options.append(mode)
    if 'slow' in mode and 'surf' in mode:
        slow_surf_mode_options.append(mode)
    if 'fast' in mode and 'surf' in mode:
        fast_surf_mode_options.append(mode)
    if 'fast' in mode and 'kink' in mode:
        fast_kink_mode_options.append(mode)
    if 'fast' in mode and 'saus' in mode:
        fast_saus_mode_options.append(mode)
    if 'fast' in mode and 'body-1' in mode:
        fast_body_1_mode_options.append(mode)
    if 'fast' in mode and 'body-2' in mode:
        fast_body_2_mode_options.append(mode)
    if 'fast' in mode and 'body-3' in mode:
        fast_body_3_mode_options.append(mode)
    if 'slow' in mode and 'body-1' in mode:
        slow_body_1_mode_options.append(mode)
    if 'slow' in mode and 'body-2' in mode:
        slow_body_2_mode_options.append(mode)
    if 'slow' in mode and 'body-3' in mode:
        slow_body_3_mode_options.append(mode)
    if 'slow' in mode and 'body' in mode:
        slow_body_mode_options.append(mode)

#fix R1 and R2
R1=1.5
R2=2. ##rho_2 / rho_0

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1) / (sf.c0**2 - W**2)

def disp_rel_asym_2var(W, K):
    return sf.disp_rel_asym(W, K, R1) / (sf.c0**2 - W**2)

Kmin = 0.01
Kmax = 10.

#if show_disp == True:
#    dispersion_diagram.density_diagram(disp_rel_asym_2var, K, W, R1, Kmin, Kmax,
#                                       just_dots=True)
#        plt.tight_layout() # seems to make it chop the sides off with this
#    plt.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_dispersion_diagrams\\'
#                + 'R1_' + str(R1) + '_' + mode + '.png')   
#    plt.close()


if show_xi_of_x == True:
    #number of iterations
    N=500
    
#    #Define zero vectors, which will be filled in the for loop
#    R1vals=np.zeros(N)
#    Wvals_saus=np.zeros(N)
#    Wvals_kink=np.zeros(N)
    
    #Choose which values of R1 you want to produce a 'xix of x' plot for.
    K_xix_vals = np.linspace(2., 6., 5)
    
    #initial K and W values
    K = Kmin
    
    
    
    
    
    font = {'size': 15}
    matplotlib.rc('font', **font)
    
    if show_disp == True:
        #density plot initiate
        plt.figure(num=None, figsize=(10, 11), facecolor='w', edgecolor='k')
        ax = plt.subplot()
    
    n = 1
    #xi of x plot initiate
    fig2 = plt.figure(figsize=(10, 3))
   
    gs = gridspec.GridSpec(n, len(K_xix_vals))#, height_ratios=[2,1]) 
    plt.subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel(r'$kx$', fontsize='x-large')
    plt.ylabel(r'$\widehat{\xi}_x(x)$', fontsize='x-large')#, labelpad=10)
    yaxmax = 1.
    nx = 1000
    z = 0.
    t = 0.
    plt.tight_layout()
    axes = []
    
    
    
    
    
    K_guess = [1.84] * n
    W_guess = [0.920]
    step = [0.001] * n
    K_start = [0.2] * n
    K_end = [Kmax] * n
    K_fast_body_circles = []
    W_fast_body_circles = []
    K_transition_circles = []
    W_transition_circles = []

        
    lw = 1.5
    num_modes = range(n)
    
    W_vals_list = []
    K_xix_N_float = (K_xix_vals-K_start[0]) / step[0]
    K_xix_N = K_xix_N_float.astype(int)
    
    plot_number = 0
    for i in num_modes:#
        mode = mode_options[i]
        K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, K_guess[i], 
                                                      W_guess[i], step[i], K_start[i],
                                                      K_end[i], (None))
        if show_disp == True:
            ax.plot(K_values, root_array, linewidth=lw)
        
        W_vals_list.append(np.zeros(len(K_xix_vals)))
        for j in range(len(K_xix_vals)):
            K = K_xix_vals[j]
            xmin = -2.*K
            xmax = 2.*K
            xvals = np.linspace(xmin, xmax, nx)
            K = K_xix_vals[j]
            W_vals_list[i][j] = newton(disp_rel_asym_1var, root_array[K_xix_N[j]], 
                                       tol=1e-5, maxiter=50)
            W = W_vals_list[i][j]
            print('W = '+str(W))
            print(root_array[K_xix_N[j]])
#            import pdb; pdb.set_trace()
            # xixvals are the \xi_x values, i.e. the displacement in the x direction.
            xixvals = sf.xix_hat(mode, xvals, W, K, R1)            
            if mode in saus_mode_options: #not sure why we need to do this but for some reason the saus xix vals are imag.
                xixvals = 1j*xixvals
            axes.append(fig2.add_subplot(gs[j+len(K_xix_vals)*i]))
            plt.tick_params(top='off', bottom='off')
            if i == 0: #top line of plots
                plt.title(r'$kx_0 = $' + str(K))
#            if j != 0: #all but left most plots
#                plt.setp(axes[plot_number].get_yticklabels(), visible=False)
#            if i != len(R1_xix_vals): #all but bottom plots
#                plt.setp(axes[plot_number].get_xticklabels(), visible=False)
            plt.setp(axes[plot_number].get_xticklabels(), visible=False)
            plt.setp(axes[plot_number].get_yticklabels(), visible=False)
            plt.plot(xvals,xixvals, 'k-')
            plt.plot((-K,-K), (-yaxmax, yaxmax), 'k--')
            plt.plot((K,K), (-yaxmax, yaxmax), 'k--')
            plt.ylim(0, yaxmax)
            plt.xlim(xmin, xmax)
            plt.xticks((-K, K))
            plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                    where=xixvals <= np.zeros_like(xvals),
                                    facecolor='blue', alpha=0.1)
            plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                    where=xixvals >= np.zeros_like(xvals),
                                    facecolor='red', alpha=0.1)
            plot_number += 1
            