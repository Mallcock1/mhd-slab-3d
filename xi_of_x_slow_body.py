# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 08:54:21 2017

@author: Matt
"""



import slab_functions as sf
import numpy as np
import scipy as sc
from scipy.optimize import newton
import matplotlib.pyplot as plt
from matplotlib import gridspec
import dispersion_diagram
import toolbox as tool
import matplotlib

##########################

#*****************REQUIRES SBB SPEEDS IN slab_functions

show_xi_of_x = False
show_disp = False

show_xi_of_x = True
#show_disp = True


##Define the sound speeds and alfven speeds.
#c0=1.
#c2=0.7
#vA=0.4
#cT=sc.sqrt(c0**2*vA**2*(c0**2+vA**2)**(-1))

mode_options = ['slow-kink-body-1', 'slow-saus-body-1', 'slow-kink-body-2', 
                'slow-saus-body-2', 'slow-kink-body-3', 'slow-saus-body-3']

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

#fix K and R2
K=6 #1.5 ## kx_0
R2=2. ##rho_2 / rho_0

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)

def disp_rel_asym_2var(W, R1):
    return sf.disp_rel_asym(W, K, R1)

R1min = 0.001
R1max = 4.

if show_disp == True:
    dispersion_diagram.density_diagram(disp_rel_asym_2var, K, W, R1, R1min, R1max,
                                       just_dots=True)
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
    R1_xix_vals = np.linspace(1., 3., 5)
    
    #initial K and W values
    R1 = R1min
    
    
    
    
    
    font = {'size': 15}
    matplotlib.rc('font', **font)
    
    if show_disp == True:
        #density plot initiate
        plt.figure(num=None, figsize=(10, 11), facecolor='w', edgecolor='k')
        ax = plt.subplot()
    
    n = 6
    #xi of x plot initiate
    fig2 = plt.figure(figsize=(10, 11))
    gs = gridspec.GridSpec(n, len(R1_xix_vals))#, height_ratios=[2,1]) 
    plt.subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.xlabel(r'$kx$', fontsize='x-large')
    plt.ylabel(r'$\widehat{\xi}_x(x)$', fontsize='x-large')#, labelpad=10)
    yaxmax = 1.
    xmin = -2.*K
    xmax = 2.*K
    nx = 1000
    z = 0.
    t = 0.
    xvals = np.linspace(xmin, xmax, nx)
    axes = []
    
    
    
    
    
    R1_guess = [0.2096] * n
    W_guess = [0.8694, 0.8040, 0.7577, 0.7289, 0.7111, 0.6998]
    step = [0.01] * n
    R1_start = [R1min] * n
    R1_end = [R1max] * n
    R1_fast_body_circles = []
    W_fast_body_circles = []
    R1_transition_circles = []
    W_transition_circles = []

        
    lw = 1.5
    num_modes = range(n)
    
    W_vals_list = []
    R1_xix_N_float = (R1_xix_vals-R1min)*N / (R1max - R1min)
    R1_xix_N = R1_xix_N_float.astype(int)
    
    plot_number = 0
    for i in num_modes:#
        mode = mode_options[i]
        R1_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, R1_guess[i], 
                                                      W_guess[i], step[i], R1_start[i],
                                                      R1_end[i], (None))
        if show_disp == True:
            ax.plot(R1_values, root_array, linewidth=lw)
        
        W_vals_list.append(np.zeros(len(R1_xix_vals)))
        for j in range(len(R1_xix_vals)):
            
            R1 = R1_xix_vals[j]
            W_vals_list[i][j] = newton(disp_rel_asym_1var, root_array[R1_xix_N[j]], 
                                       tol=1e-5, maxiter=50)
            W = W_vals_list[i][j]
            # xixvals are the \xi_x values, i.e. the displacement in the x direction.
            xixvals = sf.xix_hat(mode, xvals, W, K, R1)            
            if mode in saus_mode_options: #not sure why we need to do this but for some reason the saus xix vals are imag.
                xixvals = 1j*xixvals
            axes.append(fig2.add_subplot(gs[j+len(R1_xix_vals)*i]))
            plt.tick_params(top='off', bottom='off')
            if i == 0: #top line of plots
                plt.title(r'$\rho_1/\rho_0 = $' + str(R1))
#            if j != 0: #all but left most plots
#                plt.setp(axes[plot_number].get_yticklabels(), visible=False)
#            if i != len(R1_xix_vals): #all but bottom plots
#                plt.setp(axes[plot_number].get_xticklabels(), visible=False)
            plt.setp(axes[plot_number].get_xticklabels(), visible=False)
            plt.setp(axes[plot_number].get_yticklabels(), visible=False)
            plt.plot(xvals,xixvals, 'k-')
            plt.plot((-K,-K), (-yaxmax, yaxmax), 'k--')
            plt.plot((K,K), (-yaxmax, yaxmax), 'k--')
            plt.ylim(-yaxmax, yaxmax)
            plt.xticks((-K, K))
            plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                    where=xixvals <= np.zeros_like(xvals),
                                    facecolor='blue', alpha=0.1)
            plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                    where=xixvals >= np.zeros_like(xvals),
                                    facecolor='red', alpha=0.1)
            plot_number += 1