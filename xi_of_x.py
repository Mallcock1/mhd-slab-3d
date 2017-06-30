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

##########################

show_xi_of_x = False
show_disp = False

#show_xi_of_x = True
show_disp = True


##Define the sound speeds and alfven speeds.
#c0=1.
#c2=0.7
#vA=0.4
#cT=sc.sqrt(c0**2*vA**2*(c0**2+vA**2)**(-1))



#fix K and R2
K=6. #1.5 ## kx_0
R2=2. ##rho_2 / rho_0

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)

def disp_rel_asym_2var(W, R1):
    return sf.disp_rel_asym(W, K, R1)

if show_disp == True:
    dispersion_diagram.density_diagram(disp_rel_asym_2var, K, W, R1, 1., 3.,
                                       just_dots=True)
#        plt.tight_layout() # seems to make it chop the sides off with this
#    plt.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_dispersion_diagrams\\'
#                + 'R1_' + str(R1) + '_' + mode + '.png')   
#    plt.close()


if show_xi_of_x == True:

    #number of iterations
    N=500
    
    #Define zero vectors, which will be filled in the for loop
    R1vals=np.zeros(N)
    Wvals_saus=np.zeros(N)
    Wvals_kink=np.zeros(N)
    
    #Choose which values of R1 you want to produce a 'xix of x' plot for.
    R1_xix_vals = np.linspace(1., 3., 5)
    
    #initial K and W values
    R1min = 0.01
    R1max = 5
    R1 = R1min
    
    W_init_saus = cT - 0.0001
    W_init_kink = 0.1
    
    W = W_init_saus
    
    R1_xix_N_float = (R1_xix_vals-R1min)*N / (R1max - R1min)
    R1_xix_N = R1_xix_N_float.astype(int)
    
    #for each K value, find the root of the dispersion function, and put it in the vector Wvals_saus
    for i in range(0,N):
        Wvals_saus[i]=newton(disp_rel_asym_1var, W, tol=1e-5, maxiter=50)
        R1vals[i]=R1
        W=Wvals_saus[i]
        R1=R1+(R1max-R1min)*N**(-1)
        
    # Find sausage W values for each of the chosen R1 values
    W_xix_vals_saus = np.zeros(len(R1_xix_vals))
    for i in range(len(R1_xix_vals)):
        R1 = R1_xix_vals[i]
        W_xix_vals_saus[i] = newton(disp_rel_asym_1var, Wvals_saus[R1_xix_N[i]], 
                                   tol=1e-5, maxiter=50)
    
    #initial K and W values
    R1=R1min
    W = W_init_kink
    
    #for each K value, find the root of the dispersion function, and put it in the vector Wvals_kink
    for i in range(0,N):
        Wvals_kink[i]=newton(disp_rel_asym_1var, W, tol=1e-5, maxiter=50)
        R1vals[i]=R1
        W=Wvals_kink[i]
        R1=R1+(R1max-R1min)*N**(-1)
    
    # Find kink W values for each of the chosen R1 values
    W_xix_vals_kink = np.zeros(len(R1_xix_vals))
    for i in range(len(R1_xix_vals)):
        R1 = R1_xix_vals[i]
        W_xix_vals_kink[i] = newton(disp_rel_asym_1var, Wvals_kink[R1_xix_N[i]],
                                    tol=1e-5, maxiter=50)
    
    ###############
    ##############
    ############
    ##########
    ########
    
    fig1 = plt.figure(figsize=[5, 4])
    
    #plot them both and label axes etc
    Plot_cT=plt.plot(R1vals,cT*np.ones(N),'k-.')
    plt.plot(R1vals,Wvals_saus, 'k-', label='Quasi-sausage')
    plt.plot(R1vals,Wvals_kink, 'k--', label='Quasi-kink')
    plt.xlabel(r'$\rho_1/\rho_0$', fontsize='xx-large')
    plt.ylabel(r'$\omega/kc_0$', fontsize='xx-large')
    plt.ylim(0., vA)
    plt.annotate(r'$c_T$', xy=(0.97*R1max, cT+0.005), xytext=(0.92*R1max, cT+0.005), fontsize='x-large')
    plt.legend(loc=4)
    
    #annotate each chosen R1 value
    for i in range(len(R1_xix_vals)):
        plt.plot(R1_xix_vals, W_xix_vals_saus, 'bo')
        plt.plot(R1_xix_vals, W_xix_vals_kink, 'ro')
    
    plt.tight_layout()
    plt.show()
    #plt.savefig(u'D:/my_work/visualisations/xix_of_x_disp/disp.png')
    
    #############
    
    
    
    fig2 = plt.figure()
    gs = gridspec.GridSpec(2, len(R1_xix_vals), height_ratios=[2,1]) 
    
    plt.subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    
    plt.xlabel(r'$kx$', fontsize='x-large')
    plt.ylabel(r'$\widehat{\xi}_x(x)$', fontsize='x-large', labelpad=10)
    
    
    yaxmax = 1.
    
    xmin = -2.*K
    xmax = 2.*K
    nx = 1000
    
    z = 0.
    t = 0.
    
    xvals = np.linspace(xmin, xmax, nx)
    
    ax = []
    
    # plot the sausage xix plots
    mode = 'slow-saus-surf'
    for i in range(len(R1_xix_vals)):  
        W = W_xix_vals_saus[i]
        R1 = R1_xix_vals[i]
        
    # xixvals are the \xi_x values, i.e. the displacement in the x direction.
        xixvals = sf.xix_hat(mode, xvals, W, K, R1)
        
    
        if i == 0:
            ax.append(fig2.add_subplot(gs[i]))
        else:
            ax.append(fig2.add_subplot(gs[i]))
            plt.setp(ax[i].get_yticklabels(), visible=False)
        plt.tick_params(top='off', bottom='off')
        plt.setp(ax[i].get_xticklabels(), visible=False)
        plt.plot(xvals,xixvals, 'k-')
        plt.plot((-K,-K), (-yaxmax, yaxmax), 'k--')
        plt.plot((K,K), (-yaxmax, yaxmax), 'k--')
    
        plt.ylim(-yaxmax, yaxmax)
        plt.title(r'$\rho_1/\rho_0 = $' + str(R1))
        plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                where=xixvals <= np.zeros_like(xvals),
                                facecolor='blue', alpha=0.1)
        plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                where=xixvals >= np.zeros_like(xvals),
                                facecolor='red', alpha=0.1)
        
    # plot the kink xix plots
    mode = 'slow-kink-surf'
    for i in range(len(R1_xix_vals)):
        W = W_xix_vals_kink[i]
        R1 = R1_xix_vals[i]
    
    # xixvals are the \xi_x values, i.e. the displacement in the x direction.
        xixvals = sf.xix_hat(mode, xvals, W, K, R1)
        
        if i == 0:
            ax.append(fig2.add_subplot(gs[len(R1_xix_vals)+i]))
        else:
            ax.append(fig2.add_subplot(gs[len(R1_xix_vals)+i]))
            plt.setp(ax[len(R1_xix_vals)+i].get_yticklabels(), visible=False)
        plt.tick_params(top='off', bottom='off')
        plt.yticks(np.arange(0, 1.5, 0.5))
        plt.xticks((-1.5, 1.5))
        plt.plot(xvals,xixvals, 'k-')
        plt.plot((-K,-K), (0, yaxmax), 'k--')
        plt.plot((K,K), (0, yaxmax), 'k--')
        plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                where=xixvals <= np.zeros_like(xvals),
                                facecolor='blue', alpha=0.1)
        plt.fill_between(xvals, xixvals, np.zeros_like(xvals),
                                where=xixvals >= np.zeros_like(xvals),
                                facecolor='red', alpha=0.1)
