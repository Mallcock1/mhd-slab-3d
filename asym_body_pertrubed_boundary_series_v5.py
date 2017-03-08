
#import sys
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\ffmpeg')
##sys.path.append(u'W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/')

#import pdb # pause code for debugging at pdb.set_trace()

import numpy as np

import Toolbox as tool
import slab_functions_perturbed_boundary_v3 as sf

from pysac.plot.mayavi_seed_streamlines import SeedStreamline

import matplotlib.pyplot as plt
import matplotlib

from mayavi import mlab
#mlab.options.offscreen = True

import img2vid as i2v

###############################################################################

# What mode do you want? OPTIONS:
mode_options = ['slow-kink-surf', 'slow-saus-surf', 'slow-saus-body-3',
                'slow-kink-body-3', 'slow-saus-body-2', 'slow-kink-body-2', 
                'slow-saus-body-1', 'slow-kink-body-1', 'fast-saus-body-1',
                'fast-kink-body-1', 'fast-saus-body-2', 'fast-kink-body-2',
                'fast-saus-body-3', 'fast-kink-body-3', 'fast-kink-surf',
                'fast-saus-surf', 'shear-alfven', 'shear-alfven-broadband']

alfven_mode_options = ['shear-alfven', 'shear-alfven-broadband']
                
kink_mode_options = ['slow-kink-surf', 'slow-kink-body-1', 'slow-kink-body-2',
                     'slow-kink-body-3', 'fast-kink-body-1', 'fast-kink-body-2',
                     'fast-kink-body-3', 'fast-kink-surf']
saus_mode_options = ['slow-saus-surf', 'slow-saus-body-1', 'slow-saus-body-2',
                     'slow-saus-body-3', 'fast-saus-body-1', 'fast-saus-body-2',
                     'fast-saus-body-3', 'fast-saus-surf']
slow_surf_mode_options = ['slow-kink-surf', 'slow-saus-surf']
fast_surf_mode_options = ['fast-kink-surf', 'fast-saus-surf']
fast_kink_mode_options = ['fast-kink-surf', 'fast-kink-body-3', 'fast-kink-body-2', 
                          'fast-kink-body 1']
fast_saus_mode_options = ['fast-saus-surf', 'fast-saus-body-3', 'fast-saus-body-2', 
                          'fast-saus-body-1']
slow_body_1_mode_options = ['slow-kink-body-1', 'slow-saus-body-1']
slow_body_2_mode_options = ['slow-kink-body-2', 'slow-saus-body-2']
slow_body_3_mode_options = ['slow-kink-body-3', 'slow-saus-body-3']
slow_body_mode_options = slow_body_1_mode_options + slow_body_2_mode_options + slow_body_3_mode_options
fast_body_1_mode_options = ['fast-kink-body-1', 'fast-saus-body-1']
fast_body_2_mode_options = ['fast-kink-body-2', 'fast-saus-body-2']
fast_body_3_mode_options = ['fast-kink-body-3', 'fast-saus-body-3']



# Which angle shall we view from?
view_options = ['front', 'front-parallel', 'top', 'top-parallel' 'front-top',
                'front-side']
view = 'front'
#view = 'front-parallel'
#view = 'top'
#view = 'top-parallel'
#view = 'front-top'
#view = 'front-side'

# Uniform lighting?
#uniform_light = True
uniform_light = False

show_density = False
show_density_pert = False
show_mag = False
show_mag_scale = False
show_mag_fade = False
show_mag_vec = False
show_vel_front = False
show_vel_front_pert = False
show_vel_top = False
show_vel_top_pert = False
show_disp_top = False
show_disp_front = False
show_axes = False
show_axis_labels = False
show_mini_axis = False
show_boundary = False


# Set to True if you would like the dispersion diagram with this mode highlighted.
show_dispersion = False
#show_dispersion = True

# Wanna see the animation? Of course you do
#show_animation = False
show_animation = True

# Uncomment the parametrer you would like to see
# No density perturbations or vel/disp pert for alfven modes.
#show_density = True
#show_density_pert = True
show_mag = True
#show_mag_scale = True #must also have show_mag = True
show_mag_fade = True
#show_mag_vec = True
#show_vel_front = True
#show_vel_front_pert = True
#show_vel_top = True
#show_vel_top_pert = True
#show_disp_top = True
#show_disp_front = True
show_axes = True
#show_axis_labels = True
show_mini_axis = True
show_boundary = True

# Video resolution
#res = (1920,1080)
res = tuple(101 * np.array((16,9)))
#res = tuple(51 * np.array((16,9)))
#res = tuple(21 * np.array((16,9)))

number_of_frames = 1#50

make_video = False
#make_video = True

#
##
###
####
#####
######
#######
########
#########

#if np.array([show_density, show_vel_front_pert, show_vel_top_pert]).any()  == True:
#    raise NameError('Cannot show density or vel/disp pert for this mode')

#for mode_ind in range(14): # for all others. REMEMBER SBB pparameters
#for mode_ind in [14,15]: #for fast body surf. REMEMBER SBS parameters
for mode_ind in [14]: #for an individual mode
    if mode_ind not in range(len(mode_options)):
        raise NameError('Mode not in mode_options')
        
    # choose your mode: (note that fast surface modes, i.e. 14 and 15, can only be 
    # found with SBS parameters in slab_functions...)
    mode = mode_options[mode_ind] #All working, 0-15
#    mode = alfven_mode_options[mode_ind]
    
    print('Starting ' + mode)
    
    # Specify oscillation parameters
    if mode in slow_surf_mode_options + alfven_mode_options:
        K = 2.
    elif mode in slow_body_mode_options:
        K = 8.
    elif mode in fast_body_1_mode_options:
        K = 8.
    elif mode in fast_body_2_mode_options:
        K = 15.
    elif mode in fast_body_3_mode_options:
        K = 22.
    elif mode in fast_surf_mode_options:
        K = 8. #6.
    else:
        raise NameError('Mode not found')
            
#    R1 = 1.5 # Higher denisty on left than right
    R1 = 1.8
#    R1 = 2. # Symmetric slab
        
    def disp_rel_asym_2var(W, K):
        return sf.disp_rel_asym(W, K, R1)
    
    # find eigenfrequencies W (= omega/k) within the range Wrange for the given parameters.
    Wrange1 = np.linspace(0., sf.cT, 11)
    Wrange2 = np.linspace(sf.cT, sf.c0, 401)
    Wrange3 = np.linspace(sf.c0, sf.c2, 11)
    
    Woptions_slow_surf = np.real(tool.point_finder_scipy(disp_rel_asym_2var, np.array(K), Wrange1, args=None).transpose())
    Woptions_slow_body = np.real(tool.point_finder_scipy(disp_rel_asym_2var, np.array(K), Wrange2, args=None).transpose())
    Woptions_fast = np.real(tool.point_finder_scipy(disp_rel_asym_2var, np.array(K), Wrange3, args=None).transpose())
    
    # remove W values that are very close to characteristic speeds
    tol = 1e-2
    indices_to_rm = []
    
    for i in range(len(Woptions_slow_surf)):
        w = Woptions_slow_surf[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < 0 or w > sf.cT:
            indices_to_rm.append(i)
    Woptions_slow_surf = np.delete(Woptions_slow_surf, indices_to_rm)
    Woptions_slow_surf.sort()
    
    indices_to_rm = []
    for i in range(len(Woptions_slow_body)):
        w = Woptions_slow_body[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < sf.cT or w > sf.c0:
            indices_to_rm.append(i)
    Woptions_slow_body = np.delete(Woptions_slow_body, indices_to_rm)
    Woptions_slow_body.sort()
    
    indices_to_rm = []
    for i in range(len(Woptions_fast)):
        w = Woptions_fast[i]
        if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < sf.c0 or w > min(sf.c1, sf.c2):
            indices_to_rm.append(i)
    Woptions_fast = np.delete(Woptions_fast, indices_to_rm)
    Woptions_fast.sort()
    
    # remove any higher order slow body modes
    if len(Woptions_slow_body) > 6:
        Woptions_slow_body = np.delete(Woptions_slow_body, range(len(Woptions_slow_body) - 6))
    
    Woptions = np.concatenate((Woptions_slow_surf, Woptions_slow_body, Woptions_fast))
    
    
    # set W to be the eigenfrequency for the requested mode
    if mode in ['fast-saus-body-1', 'fast-saus-body-2', 'fast-saus-body-3', 'fast-kink-surf']:
        W = Woptions_fast[-2]
    elif mode in ['fast-kink-body-1', 'fast-kink-body-2', 'fast-kink-body-3', 'fast-saus-surf']:
        W = Woptions_fast[-1]
    elif mode in slow_surf_mode_options:
        W = Woptions_slow_surf[mode_ind]
    elif mode in slow_body_mode_options:
        W = Woptions_slow_body[mode_ind-2]
    
    if mode in alfven_mode_options:
        W = sf.vA
    else:
        W = np.real(W)
        
# Quick plot to see if we are hitting correct mode    
#    plt.plot([K] * len(Woptions), Woptions, '.')
#    plt.plot(K+0.5, W, 'go')
#    plt.xlim([0,23])
    
    # Dependent variables:
    # x = k*x
    # y = k*y
    # z = k*z
    # W = omega/k
    # K = k*x_0
    # t = omega*t
    
    #################################################################################
    
    if show_dispersion == True:
        if mode in alfven_mode_options:
            raise NameError('Disperion plot requested for an alfven mode. Cant do it')
        # Plot the dispersion diagram with the chosen mode highlighted
        
        Kmin = 0.
        Kmax = 23.#10.
        
        Wmin = 0.
        Wmax = sf.c2
        
        K_range = np.linspace(Kmin, Kmax, 51)
        W_range = np.linspace(Wmin, Wmax, 51)
        
        
        font = {'size': 15}
        matplotlib.rc('font', **font)
        
        plt.figure(num=None, figsize=(10, 11), facecolor='w', edgecolor='k')
        ax = plt.subplot()
        
        # colours 
        def colour(mode):
            colours = []
            if mode in fast_surf_mode_options:
                for i in range(len(mode_options)):
                    if i <= 7:
                        if mode_options[i] == mode:
                            colours.append('g')
                        elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                            colours.append('r')
                        else:
                            colours.append('b')
                    elif i >= 14:
                        if mode_options[i] == mode:
                            colours.append('g')
                        elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                            colours.append('r')
                        else:
                            colours.append('b')
            else:
                for i in range(len(mode_options)):
                    if mode_options[i] == mode:
                        colours.append('g')
                    elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                        colours.append('r')
                    else:
                        colours.append('b')
            return colours
         
        def ls(mode):
            linestyles = []
            if mode in fast_surf_mode_options:
                for i in range(8):
                    if mode_options[i] in kink_mode_options:
                        linestyles.append('--')
                    else:
                        linestyles.append('-')
                for i in range(8, len(mode_options)):
                    if mode_options[i-1] in kink_mode_options:
                        linestyles.append('--')
                    else:
                        linestyles.append('-')
            else:
                for i in range(len(mode_options)):
                    if mode_options[i] in kink_mode_options:
                        linestyles.append('--')
                    else:
                        linestyles.append('-')
            return linestyles
            
        
#        #Plot the dots
#        W_array = tool.point_finder_scipy(disp_rel_asym_2var, K_range, W_range,
#                                          args=(None))
#        ax.plot(K_range, W_array, '.', color = 'b')
        
        
        if R1 == 1.5:
            if mode in fast_surf_mode_options:
                K_guess = [1., 1., 10.12, 8.74, 5.98, 5.98, 3.68, 1.84, 6., 6.]
                W_guess = [0.64, 0.72, 0.863, 0.867, 0.855, 0.885, 0.885,
                           0.920, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 4.22, 
                           0.67]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 4.2, Kmax, 
                         Kmax]
                K_fast_body_circles = [-10, 0.67]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = [4.2]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 1.84,
                           6.44, 9.2, 15.18, 19.78, 22.01]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.085, 1.131, 1.181,
                           1.154, 1.156, 1.181]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.47, 
                           4.51, 8.55, 12.55, 16.5, 20.53]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [0.47, 4.51, 8.55, 12.55, 16.5, 20.53]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []
 
                
        elif R1 == 1.8:
            if mode in fast_surf_mode_options:
                K_guess = [1., 1., 10.12, 8.74, 5.98, 5.98, 3.68, 1.84, 6., 6.]
                W_guess = [0.64, 0.72, 0.863, 0.867, 0.855, 0.885, 0.885,
                           0.920, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 4.752, 
                           0.365]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 4.752, Kmax, 
                         Kmax]
                K_fast_body_circles = [-10, 0.365]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = [4.752]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 1.84,
                           6.44, 9.2, 15.18, 19.78, 22.01]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.085, 1.131, 1.181,
                           1.154, 1.156, 1.181]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.336,
                           4.306, 8.346, 12.318, 16.309, 20.342]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [0.336, 4.306, 8.346, 12.318, 16.309,
                                       20.342]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []

        elif R1 == 2.:
            if mode in fast_surf_mode_options:
                K_guess = [1., 1., 11.96, 6., 4., 4., 2., 2., 6., 6.]
                W_guess = [0.64, 0.72, 0.877, 0.836, 0.826, 0.847, 0.83, 
                           0.915, 1.01567, 1.05754]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01]
                K_start = [0.0001, 0.0001, 2., 1.6, 1.2, 0.9, 0.6, 0.2, 5.24, 
                           0.001]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, 5.2, Kmax, 
                         Kmax]
                K_fast_body_circles = []
                W_fast_body_circles = []
                K_transition_circles = [5.2]
                W_transition_circles = [sf.c0]
            else:
                K_guess = [1., 1., 8.74, 8.74, 5.98, 5.98, 3.68, 1.84, 2., 6., 
                           9., 14.73, 19.8, 23.4]
                W_guess = [0.473, 0.547, 0.719463, 0.733764, 0.722183, 
                           0.746659, 0.741942, 0.749326, 1.16675, 1.156, 
                           1.17, 1.1547, 1.16, 1.16]
                step = [0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                        0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
                K_start = [0.0001, 0.0001, 2.3, 1.8, 1.3, 1., 0.8, 0.3, 0.0001, 
                           4.1, 8.2, 12.2, 16.2, 20.2]
                K_end = [Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax, Kmax,
                         Kmax, Kmax, Kmax, Kmax, Kmax]
                K_fast_body_circles = [-10, 4.1, 8.2, 12.2, 16.2, 20.2]
                W_fast_body_circles = [sf.c2] * len(K_fast_body_circles)
                K_transition_circles = []
                W_transition_circles = []

        
        lw = 1.5
        disp_modes = range(len(K_guess))
#        disp_modes = [-1]
        for i in disp_modes:
            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, K_guess[i], 
                                                         W_guess[i], step[i], K_start[i],
                                                         K_end[i], (None))
            ax.plot(K_values, root_array, color=colour(mode)[i], linewidth=lw, 
                    linestyle=ls(mode)[i])
        for i in range(len(K_fast_body_circles)):
            ax.plot(K_fast_body_circles[i], W_fast_body_circles[i], marker='o', 
                    markerfacecolor='None', markeredgecolor=colour(mode)[i+8], 
                    markeredgewidth=lw, markersize=8)
        for i in range(len(K_transition_circles)):
            ax.plot(K_transition_circles[i], W_transition_circles[i], marker='o',
                    markerfacecolor='None', markeredgecolor=colour(mode)[i+8],
                    markeredgewidth=lw, markersize=8)
        
        ax.set_ylabel(r'$v_{ph}/c_0$', fontsize = 20)
        ax.set_xlabel(r'$k x_0$', fontsize = 20)
        
        ax.set_xlim(K_range[0], K_range[-1])
        ax.set_ylim(0., 1.41)
        
        if mode in fast_surf_mode_options:
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.42, 1.42),
                edgecolor='gray', linestyle='-.', color='None', hatch='/',
                linewidth=2)

            #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='-.', linewidth=2)
            ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
#            ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
            if R1 == 2.:
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            else:
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
        else:
            K_range_for_fill = np.append(K_range, K_range[-1]+0.1)
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c0,sf.c0), (sf.vA,sf.vA),
                            edgecolor='gray', linestyle='-.', color='None', hatch='/',
                            linewidth=2)
            ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.42, 1.42),
                            edgecolor='gray', linestyle='-.', color='None', hatch='/',
                            linewidth=2)
            
            #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='-.', linewidth=2)
            ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
            ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
    #            ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
            if R1 == 2.:
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
            else:
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)

                    
        ax.plot(K, W, 'go', markersize=10)
#        plt.tight_layout() # seems to make it chop the sides off with this
        plt.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_dispersion_diagrams\\'
                    + 'R1_' + str(R1) + '_' + mode + '.png')   
#        plt.close()
                    

    
    ##############################################################################
    
    if show_animation == True:
        xmin = -2.*K
        xmax = 2.*K
        ymin = 0.
        ymax = 4.
        zmin = 0.
        zmax = 2*np.pi
        
        # You can change ny but be careful changing nx, nz.
        nx = 100 #100
        ny = 100 #100#20 #100
        nz = 100 #100
        nt = number_of_frames
        
        if nz % nt != 0:
            print("nt doesnt divide nz so there may be a problem with chopping in z direction for each time step")
        
        t_start = 0.
        t_end = zmax
        
        t = t_start
        
        xvals = np.linspace(xmin, xmax, nx)
        yvals = np.linspace(ymin, ymax, ny)
        zvals = np.linspace(zmin, zmax, nz, endpoint=False)
        
        x_spacing = max(nx, ny, nz) / nx
        y_spacing = max(nx, ny, nz) / ny
        z_spacing = max(nx, ny, nz) / nz
        
        
        # For masking points
        mod = int(3 * nx / 100)
        mod_top = int(np.ceil(mod / y_spacing))
        
        if show_disp_top == True or show_disp_front == True:
            xixvals = np.real(np.repeat(sf.xix(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xizvals = np.real(np.repeat(sf.xiz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xiyvals = np.real(np.repeat(sf.xiz(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
        
                                
        if show_vel_front == True or show_vel_top == True:
            vxvals = np.real(np.repeat(sf.vx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.real(np.repeat(sf.vy(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
        
        
        if show_vel_front_pert == True or show_vel_top_pert == True:
            vxvals = np.real(np.repeat(sf.vx_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.zeros_like(vxvals)
        
        bxvals = np.real(np.repeat(sf.bx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        byvals = np.real(np.repeat(sf.by(mode, xvals, zvals, t, K)[:, :, np.newaxis], ny, axis=2))
        bz_eq3d = np.repeat(sf.bz_eq(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2)
        bzvals = np.real(np.repeat(-sf.bz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2) +
                         bz_eq3d)
                         
        
        if show_boundary == True:
            xix_boundary_r_vals = np.real(np.repeat(K + sf.xix_boundary(mode, zvals, t, W, K, R1, boundary='r')[:, np.newaxis], ny, axis=1))
            xix_boundary_l_vals = np.real(np.repeat(-K + sf.xix_boundary(mode, zvals, t, W, K, R1, boundary='l')[:, np.newaxis], ny, axis=1))
        
        if show_density == True:
            rho_vals = np.real(np.repeat(sf.rho(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        
        if show_density_pert == True:
            rho_vals = np.real(np.repeat(sf.rho_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        
        
        
        
        for t_ind in range(nt):
            if t_ind == 0:
                bxvals_t = bxvals
                byvals_t = byvals
                bzvals_t = bzvals
                
                if show_disp_top == True or show_disp_front == True:
                    xixvals_t = xixvals
                    xiyvals_t = xiyvals
                    xizvals_t = xizvals
                    
                if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:
                    vxvals_t = vxvals
                    vyvals_t = vyvals
                    vzvals_t = vzvals
                
                if show_boundary == True:
                    xix_boundary_r_vals_t = xix_boundary_r_vals
                    xix_boundary_l_vals_t = xix_boundary_l_vals
                    
                if show_density == True or show_density_pert == True:
                    rho_vals_t = rho_vals
        
                
            else:
                bxvals_split = np.split(bxvals, [nz - (nz / nt) * t_ind], axis=1)
                byvals_split = np.split(byvals, [nz - (nz / nt) * t_ind], axis=1)
                bzvals_split = np.split(bzvals, [nz - (nz / nt) * t_ind], axis=1)
                
                bxvals_t = np.concatenate((bxvals_split[1], bxvals_split[0]), axis=1)
                byvals_t = np.concatenate((byvals_split[1], byvals_split[0]), axis=1)
                bzvals_t = np.concatenate((bzvals_split[1], bzvals_split[0]), axis=1)
                
                
                if show_disp_top == True or show_disp_front == True:            
                    xixvals_split = np.split(xixvals, [nz - (nz / nt) * t_ind], axis=1)
                    xiyvals_split = np.split(xiyvals, [nz - (nz / nt) * t_ind], axis=1)
                    xizvals_split = np.split(xizvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    xixvals_t = np.concatenate((xixvals_split[1], xixvals_split[0]), axis=1)
                    xiyvals_t = np.concatenate((xiyvals_split[1], xiyvals_split[0]), axis=1)
                    xizvals_t = np.concatenate((xizvals_split[1], xizvals_split[0]), axis=1)
                                        
                if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:            
                    vxvals_split = np.split(vxvals, [nz - (nz / nt) * t_ind], axis=1)
                    vyvals_split = np.split(vyvals, [nz - (nz / nt) * t_ind], axis=1)
                    vzvals_split = np.split(vzvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    vxvals_t = np.concatenate((vxvals_split[1], vxvals_split[0]), axis=1)
                    vyvals_t = np.concatenate((vyvals_split[1], vyvals_split[0]), axis=1)
                    vzvals_t = np.concatenate((vzvals_split[1], vzvals_split[0]), axis=1)
                
                if show_boundary == True:
                    xix_boundary_r_vals_split = np.split(xix_boundary_r_vals, [nz - (nz / nt) * t_ind], axis=0)
                    xix_boundary_l_vals_split = np.split(xix_boundary_l_vals, [nz - (nz / nt) * t_ind], axis=0)
        
                    xix_boundary_r_vals_t = np.concatenate((xix_boundary_r_vals_split[1], xix_boundary_r_vals_split[0]), axis=0)
                    xix_boundary_l_vals_t = np.concatenate((xix_boundary_l_vals_split[1], xix_boundary_l_vals_split[0]), axis=0)
                
                if show_density == True or show_density_pert == True:            
                    rho_vals_split = np.split(rho_vals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    rho_vals_t = np.concatenate((rho_vals_split[1], rho_vals_split[0]), axis=1)                         
            
            #Masking points
            if show_mag_vec == True:
                bxvals_mask_front_t = np.copy(bxvals_t)
                byvals_mask_front_t = np.copy(byvals_t)
                bzvals_mask_front_t = np.copy(bzvals_t)
                
                for i in range(bxvals_t.shape[0]):
                    for j in range(bxvals_t.shape[1]):
                        for k in range(bxvals_t.shape[2]):
                            if (i%mod) != 1 or (j%mod) != 1:
                                bxvals_mask_front_t[i,j,k] = 0.
                                bzvals_mask_front_t[i,j,k] = 0.
        
        
            if show_disp_top == True:    
                xixvals_mask_top_t = np.copy(xixvals_t)
                xiyvals_mask_top_t = np.copy(xiyvals_t)
                xizvals_mask_top_t = np.copy(xizvals_t)
                
                for i in range(xixvals_t.shape[0]):
                    for j in range(xixvals_t.shape[1]):
                        for k in range(xixvals_t.shape[2]):
                            if (i%mod) != 1 or (k%mod_top) != 1:
                                xixvals_mask_top_t[i,j,k] = 0.
                                xizvals_mask_top_t[i,j,k] = 0.
                                
            if show_disp_front == True:
                xixvals_mask_front_t = np.copy(xixvals_t)
                xiyvals_mask_front_t = np.copy(xiyvals_t)
                xizvals_mask_front_t = np.copy(xizvals_t)
                
                for i in range(xixvals_t.shape[0]):
                    for j in range(xixvals_t.shape[1]):
                        for k in range(xixvals_t.shape[2]):
                            if (i%mod)!=1 or (j%mod)!=1:
                                xixvals_mask_front_t[i,j,k] = 0.
                                xizvals_mask_front_t[i,j,k] = 0.    
            
            
            if show_vel_top == True or show_vel_top_pert == True:    
                vxvals_mask_top_t = np.copy(vxvals_t)
                vyvals_mask_top_t = np.copy(vyvals_t)
                vzvals_mask_top_t = np.copy(vzvals_t)
                
                for i in range(vxvals_t.shape[0]):
                    for j in range(vxvals_t.shape[1]):
                        for k in range(vxvals_t.shape[2]):
                            if (i%mod) != 1 or (k%mod_top) != 1:
                                vxvals_mask_top_t[i,j,k] = 0
                                vyvals_mask_top_t[i,j,k] = 0
                                vzvals_mask_top_t[i,j,k] = 0
                                                  
                                
            if show_vel_front == True or show_vel_front_pert == True:
                vxvals_mask_front_t = np.copy(vxvals_t)
                vyvals_mask_front_t = np.copy(vyvals_t)
                vzvals_mask_front_t = np.copy(vzvals_t)
                
                for i in range(vxvals_t.shape[0]):
                    for j in range(vxvals_t.shape[1]):
                        for k in range(vxvals_t.shape[2]):
                            if (i%mod) != 1 or (j%mod) != 1:
                                vxvals_mask_front_t[i,j,k] = 0
                                vzvals_mask_front_t[i,j,k] = 0   
            
        
            zvals, yvals = np.mgrid[0:nz:(nz)*1j,
                                    0:ny:(ny)*1j]
        
            #
            ##
            ###
            ####
            #####
            ####
            ###
            ##
            #
            
            fig = mlab.figure(size=res) # (1920, 1080) for 1080p , tuple(101 * np.array((16,9))) #16:9 aspect ratio for video upload
            
            
            
            spacing =  np.array([x_spacing, z_spacing, y_spacing])
            
            if show_boundary == True: # maybe +1 after nx???? did have this, now removed. but still unsure.
                ext_min_r = ((nx) * (xix_boundary_r_vals_t.min() - xmin) / (xmax - xmin)) * x_spacing
                ext_max_r = ((nx) * (xix_boundary_r_vals_t.max() - xmin) / (xmax - xmin)) * x_spacing
                
                ext_min_l = ((nx) * (xix_boundary_l_vals_t.min() - xmin) / (xmax - xmin)) * x_spacing #plus 2 after (xmax-xmin)?
                ext_max_l = ((nx) * (xix_boundary_l_vals_t.max() - xmin) / (xmax - xmin)) * x_spacing #plus 2 after (xmax-xmin)?
                                                   
                if view == 'front-parallel':
                    lut = np.reshape(np.array([150, 150, 150, 255]*256), (256,4))
                    fade_value = 125
                    lut[:fade_value,-1] = np.linspace(0, 255, fade_value)
                    lut[-fade_value:,-1] = np.linspace(255, 0, fade_value)
                    
                    boundary_r_thick = mlab.mesh(xix_boundary_r_vals_t, zvals, yvals,
                                                 extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                                 opacity=1., representation='wireframe',
                                                 line_width=12., scalars=zvals)
#                    boundary_r_thick.enable_contours = True
                    boundary_l_thick = mlab.mesh(xix_boundary_l_vals_t, zvals, yvals,
                                                 extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                                 opacity=1., representation='wireframe',
                                                 line_width=12., scalars=zvals)
#                    boundary_l_thick.enable_contours = True
                    
                    boundary_r_thick.module_manager.scalar_lut_manager.lut.table = lut
                    boundary_l_thick.module_manager.scalar_lut_manager.lut.table = lut
                    boundary_r_thick.actor.property.lighting = False
                    boundary_r_thick.actor.property.shading = False
                    boundary_l_thick.actor.property.lighting = False
                    boundary_l_thick.actor.property.shading = False
                                                 
                else:
                    lut = np.reshape(np.array([150, 150, 150, 255]*256), (256,4))
                    fade_value = 20
                    lut[:fade_value,-1] = np.linspace(0, 255, fade_value)
                    lut[-fade_value:,-1] = np.linspace(255, 0, fade_value)
                    
                    boundary_r = mlab.mesh(xix_boundary_r_vals_t, zvals, yvals,
                                           extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                           opacity=0.7, scalars=zvals)

                    boundary_l = mlab.mesh(xix_boundary_l_vals_t, zvals, yvals,
                                           extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                           opacity=0.7, scalars=zvals)

                    boundary_r.module_manager.scalar_lut_manager.lut.table = lut
                    boundary_l.module_manager.scalar_lut_manager.lut.table = lut
                    boundary_r.actor.property.lighting = False
                    boundary_r.actor.property.shading = False
                    boundary_l.actor.property.lighting = False
                    boundary_l.actor.property.shading = False                                           
                                           
                                                 
            if show_density == True or show_density_pert == True:
                # Scalar field density   
                rho = mlab.pipeline.scalar_field(rho_vals_t, name="density", figure=fig)
                #scalar_cut_plane = ScalarCutPlane()
                #fig.parent.add_filter(scalar_cut_plane, sca)
                rho.spacing = spacing
                minr = rho_vals_t.min()
                maxr = rho_vals_t.max()
                
                #Volume for high pressure
                rvmin1 = minr + 0.5 * (maxr - minr)
                rvmax1 = minr + 1. * (maxr - minr)
                rvol1 = mlab.pipeline.volume(rho, vmin=rvmin1, vmax=rvmax1)
                
                # Changing the ctf:
                from tvtk.util.ctf import ColorTransferFunction
                ctf1 = ColorTransferFunction()
                ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
                ctf1.add_rgb_point(rvmin1 + 0.4 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
                ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
                # ...
                rvol1._volume_property.set_color(ctf1)
                rvol1._ctf = ctf1
                rvol1.update_ctf = True
                
                #Changing the opacity of the volume vol1
                ## Changing the otf:
                from tvtk.util.ctf import PiecewiseFunction
                otf = PiecewiseFunction()
                otf.add_point(rvmin1, 0)
                otf.add_point(rvmin1 + (rvmax1-rvmin1)*0.2, 0.012)
                otf.add_point(rvmin1 + (rvmax1-rvmin1)*0.5, 0.05)
                otf.add_point(rvmax1, 0.15)
                ##vol1._otf = otf
                rvol1._volume_property.set_scalar_opacity(otf)
                
                # exempt volume from shading and improve overall look by increasing opacity
                rvol1.volume_property.shade = False
                rvol1.volume_property.scalar_opacity_unit_distance = 2.0
                
                
                #Volume for low pressure
                rvmin2 = minr + 0. * (maxr - minr)
                rvmax2 = minr + 0.5 * (maxr - minr)
                rvol2 = mlab.pipeline.volume(rho, vmin=rvmin2, vmax=rvmax2)
                
                # Changing the ctf:
                ctf2 = ColorTransferFunction()
                ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
                ctf2.add_rgb_point(rvmin2 + 0.6 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
                ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
                # ...
                rvol2._volume_property.set_color(ctf2)
                rvol2._ctf = ctf2
                rvol2.update_ctf = True
            
                #Changing the opacity of the volume vol2
                ## Changing the otf:
                otf = PiecewiseFunction()
                otf.add_point(rvmax2, 0)
                otf.add_point(rvmax2 - (rvmax2-rvmin2)*0.2, 0.012)
                otf.add_point(rvmax2 - (rvmax2-rvmin2)*0.5, 0.05)
                otf.add_point(rvmin2, 0.15)
                ##vol1._otf = otf
                rvol2._volume_property.set_scalar_opacity(otf)
                
                # exempt volume from shading and improve overall look by increasing opacity
                rvol2.volume_property.shade = False
                rvol2.volume_property.scalar_opacity_unit_distance = 2.0
            
                
            
            # Vector field bxvals, bzvals, byvals
#            field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
#                                                   figure=fig)
#            field.spacing = spacing
#                
#            if show_mag == True:
#                #contours = mlab.pipeline.iso_surface(magnitude,
#                #                                        contours=range(2, 14, 3),
#                #                                        transparent=True,
#                #                                        opacity=0.4,
#                #                                        colormap='YlGnBu',
#                #                                        vmin=0, vmax=14)
#                
#                # Create an array of points for which we want mag field seeds
#                nx_seed = 9 #7
#                ny_seed = 13 #10
#                start_x = 30. #38
#                end_x = nx+1 - start_x
#                start_y = 1.
#                if ny == 20:
#                    end_y = ny - 1 #ny-2 for ny = 100
#                elif ny == 100:
#                    end_y = ny - 2
#                else:
#                    end_y = ny - 1
#                seeds=[]
#                dx_res = (end_x - start_x) / (nx_seed-1)
#                dy_res = (end_y - start_y) / (ny_seed-1)
#                for j in range(0,ny_seed):
#                    for i in range(0,nx_seed):
#                        x = start_x + (i * dx_res) * x_spacing
#                        y = start_y + (j * dy_res) * y_spacing
#                        z = 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
#                        seeds.append((x,z,y))
#                
#                if mode in alfven_mode_options:
#                    for i in range(nx_seed):
#                        del seeds[0]
#                        del seeds[-1]
#                                                             
#                field_lines = SeedStreamline(seed_points=seeds)
#                field_lines.stream_tracer.integration_direction='both'
#                field_lines.streamline_type = 'tube'
#                
#                magnitude = mlab.pipeline.extract_vector_norm(field)
#                magnitude.add_child(field_lines)
#                module_manager = field_lines.parent
#                module_manager.scalar_lut_manager.lut_mode = 'Reds'
#                module_manager.scalar_lut_manager.data_range=[-30,25]
            
#                field_lines.stream_tracer.maximum_propagation = 500.
#                field_lines.tube_filter.number_of_sides = 20
#                field_lines.tube_filter.radius = 0.7
#                field_lines.tube_filter.capping = True
#                
#                if show_mag_scale == True:
#                    module_manager.scalar_lut_manager.lut_mode = 'jet'
#                    module_manager.scalar_lut_manager.data_range=[7,18]
            
            xvals, zvals, yvals = np.mgrid[0:nx:(nx)*1j,
                                           0:nz:(nz)*1j,
                                           0:ny:(ny)*1j]
                
            field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
                                                   figure=fig, scalars=zvals)
            field.spacing = spacing
                
                #contours = mlab.pipeline.iso_surface(magnitude,
                #                                        contours=range(2, 14, 3),
                #                                        transparent=True,
                #                                        opacity=0.4,
                #                                        colormap='YlGnBu',
                #                                        vmin=0, vmax=14)
            if show_mag == True:
                # Create an array of points for which we want mag field seeds
                nx_seed = 9 #7
                ny_seed = 13 #10
                start_x = 30. #38
                end_x = nx+1 - start_x
                start_y = 1.
                if ny == 20:
                    end_y = ny - 1 #ny-2 for ny = 100
                elif ny == 100:
                    end_y = ny - 2
                else:
                    end_y = ny - 1
                seeds=[]
                dx_res = (end_x - start_x) / (nx_seed-1)
                dy_res = (end_y - start_y) / (ny_seed-1)
                for j in range(0,ny_seed):
                    for i in range(0,nx_seed):
                        x = start_x + (i * dx_res) * x_spacing
                        y = start_y + (j * dy_res) * y_spacing
                        z = 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                        seeds.append((x,z,y))
                
                if mode in alfven_mode_options:
                    for i in range(nx_seed):
                        del seeds[0]
                        del seeds[-1]
                                                             
                field_lines = SeedStreamline(seed_points=seeds)
                field_lines.stream_tracer.integration_direction='both'
                field_lines.streamline_type = 'tube'
                
#                module_manager = field_lines.parent
#                module_manager.scalar_lut_manager.lut_mode = 'Reds'
#                module_manager.scalar_lut_manager.data_range=[-30,25]    
                
#                magnitude = mlab.pipeline.add_dataset(field)
                field.add_child(field_lines)
                module_manager = field_lines.parent
#                module_manager.scalar_lut_manager.lut_mode = 'Reds'
#                module_manager.scalar_lut_manager.data_range=[-30,25]
                
                field_lines.stream_tracer.maximum_propagation = 500.
                field_lines.tube_filter.number_of_sides = 20
                field_lines.tube_filter.radius = 0.7
                field_lines.tube_filter.capping = True
                field_lines.actor.property.opacity = 1.0
                
                if show_mag_scale == True:
                    module_manager.scalar_lut_manager.lut_mode = 'jet'
                    module_manager.scalar_lut_manager.data_range=[7,18]
                else:
                    mag_lut = module_manager.scalar_lut_manager.lut.table.to_array()
                    mag_lut[:,0] = [220]*256
                    mag_lut[:,1] = [20]*256
                    mag_lut[:,2] = [20]*256
                    module_manager.scalar_lut_manager.lut.table = mag_lut
#                    module_manager.scalar_lut_manager.data_range=[-1500,500]
                if show_mag_fade == True:
                    mag_lut = module_manager.scalar_lut_manager.lut.table.to_array()
#                    mag_fade_value = fade_value
#                    mag_lut[:mag_fade_value,-1] = np.linspace(0, 255, mag_fade_value)
#                    mag_lut[-mag_fade_value:,-1] = np.linspace(255, 0, mag_fade_value)
#                    module_manager.scalar_lut_manager.lut.table = mag_lut
            
            scalefactor = 4. # scale factor for direction field vectors
            
            if show_mag_vec == True:
                bdirfield_front = mlab.pipeline.vector_field(bxvals_mask_front_t, bzvals_mask_front_t,
                                                             byvals_mask_front_t, name="B field front",
                                                             figure=fig)
                bdirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(bdirfield_front, 
                                                                  scale_factor=scalefactor)
                vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
                vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
                vector_cut_plane_front.glyph.color_mode = 'no_coloring'
                vector_cut_plane_front.implicit_plane.widget.enabled = False
                vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
            
            
            if show_vel_top == True or show_vel_top_pert == True:
                vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top_t, np.zeros_like(vxvals_mask_top_t),
                                                            vyvals_mask_top_t, name="V field top",
                                                            figure=fig)
                vdirfield_top.spacing = spacing
                vector_cut_plane_top = mlab.pipeline.vector_cut_plane(vdirfield_top, 
                                                                  scale_factor=scalefactor)
                vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
                vector_cut_plane_top.glyph.color_mode = 'no_coloring'
                vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-0.1, 50.5])
                vector_cut_plane_top.implicit_plane.widget.enabled = False
                vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
                
            if show_vel_front == True or show_vel_front_pert == True:
                vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front_t, vzvals_mask_front_t,
                                                             vyvals_mask_front_t, name="V field front",
                                                             figure=fig)
                vdirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                                  scale_factor=scalefactor)
                vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
                vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
                vector_cut_plane_front.glyph.color_mode = 'no_coloring'
                vector_cut_plane_front.implicit_plane.widget.enabled = False
                vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
            
            if show_disp_top == True:
                xidirfield_top = mlab.pipeline.vector_field(xixvals_mask_top_t, np.zeros_like(xixvals_mask_top_t),
                                                            xiyvals_mask_top_t, name="Xi field top",
                                                            figure=fig)
                xidirfield_top.spacing = spacing
                vector_cut_plane_top = mlab.pipeline.vector_cut_plane(xidirfield_top, 
                                                                  scale_factor=scalefactor)
                vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
                vector_cut_plane_top.glyph.color_mode = 'no_coloring'
                vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-0.1, 50.5])
                vector_cut_plane_top.implicit_plane.widget.enabled = False
                vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
                
            if show_disp_front == True:
                xidirfield_front = mlab.pipeline.vector_field(xixvals_mask_front_t, xizvals_mask_front_t,
                                                             xiyvals_mask_front_t, name="Xi field front",
                                                             figure=fig)
                xidirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(xidirfield_front, 
                                                                  scale_factor=scalefactor)
                vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
                vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
                vector_cut_plane_front.glyph.color_mode = 'no_coloring'
                vector_cut_plane_front.implicit_plane.widget.enabled = False
                vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
            
            #Set viewing angle
            if view == 'front-parallel':
                field.scene.z_plus_view()
                field.scene.parallel_projection = True
                field.scene.camera.zoom(1.65) # Parallel projection zoom is done in this way, different to perspective projection
            if view == 'front':
                field.scene.z_plus_view()
                field.scene.camera.view_angle = 21.
            if view == 'top':
                field.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
                field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
                field.scene.camera.view_angle = 14.
                field.scene.camera.view_up = [-0, 0, -1]
                field.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
            if view == 'top-parallel':
                field.scene.parallel_projection = True
                field.scene.camera.zoom(2.)
                field.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
                field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
    #            field.scene.camera.view_angle = 14.
                field.scene.camera.view_up = [-0, 0, -1]
                field.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
            if view == 'front-top':
                field.scene.camera.position = [48.764852970361503, 223.64895482756552, 498.62216293273576]
                field.scene.camera.focal_point = [50.821544647216797, 46., 50.159849926829338]
                field.scene.camera.view_angle = 16.0
                field.scene.camera.view_up = [-0.002418791139063777, 0.93281530024654913, -0.36034672896443193]
                field.scene.camera.clipping_range = [345.97885880654962, 650.71850659694883]
             
            if view == 'front-side':
                field.scene.camera.position = [126.6 * nx / 100., 60.5 * nz / 100., 524.8 * ny / 100.]
                field.scene.camera.focal_point = [50.8 * nx / 100., 50.4 * nz / 100., 50.2 * ny / 100.]
                field.scene.camera.view_angle = 14.
                field.scene.camera.view_up = [-0.01695 * nx / 100., 0.999686 * nz / 100., -0.0184509 * ny / 100.]
                field.scene.camera.clipping_range = [366.21083458278804, 631.07664372567524]
            
            if show_axes == True:
                axes = mlab.axes(field, nb_labels=1, line_width=3)
                axes.axes.label_format = ''
                if show_axis_labels == True:
                    axes.axes.x_label = 'x'
                    if view == 'front-parallel':
                        axes.axes.y_label = 'z'
                        axes.axes.z_label = ''
                    elif view == 'top-parallel':
                        axes.axes.y_label = ''
                        axes.axes.z_label = 'y'
                    else:
                        axes.axes.y_label = 'z'
                        axes.axes.z_label = 'y'
                else:
                    axes.axes.x_label = ''
                    axes.axes.y_label = ''
                    axes.axes.z_label = ''
                    
            if show_mini_axis == True:
                
                oa = mlab.orientation_axes(xlabel='x', ylabel='z', zlabel='y')
                oa.marker.set_viewport(0,0,0.25,0.25) # minx, miny, maxx, maxy
    
            if uniform_light == True:
                #uniform lighting, but if we turn shading of volumes off, we are ok without
                field.scene.light_manager.number_of_lights = 6
                
                camera_light1 = field.scene.light_manager.lights[0]
                camera_light1.activate = True
                camera_light1.intensity = 0.7
                camera_light1.elevation = 90.
                camera_light1.azimuth = 0.
            
                camera_light2 = field.scene.light_manager.lights[1]
                camera_light2.activate = True
                camera_light2.intensity = 0.7
                camera_light2.elevation = -90.
                camera_light2.azimuth = 0.
            
                camera_light3 = field.scene.light_manager.lights[2]
                camera_light3.activate = True
                camera_light3.intensity = 0.7
                camera_light3.elevation = 0.
                camera_light3.azimuth = -90
            
                camera_light4 = field.scene.light_manager.lights[3]
                camera_light4.activate = True
                camera_light4.intensity = 0.7
                camera_light4.elevation = 0.
                camera_light4.azimuth = 0.
            
                camera_light5 = field.scene.light_manager.lights[4]
                camera_light5.activate = True
                camera_light5.intensity = 0.7
                camera_light5.elevation = 0.
                camera_light5.azimuth = 90.
            
                camera_light6 = field.scene.light_manager.lights[5]
                camera_light6.activate = True
                camera_light6.intensity = 0.7
                camera_light6.elevation = 0.
                camera_light6.azimuth = 180.
            
            #Black background
            field.scene.background = (0., 0., 0.)
    
            
            
    # Trying and failing to sort out memory issues.
    #        mlab.gcf()
    #        mlab.clf()
    #        mlab.close()
    #        gc.collect()
    #        del fig
    #        engine_manager.current_engine = None
    #        registry.engines = {}
            
            if make_video == True:
                prefix = 'R1_'+str(R1)+'_'+view + '_' + mode
                mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
                             + prefix + str(t_ind+1) + '.png')
                mlab.close(fig)
            
        t = t + (t_end - t_start) / nt
        
    if make_video == True:
        i2v.image2video(prefix=prefix, output_name=prefix+'_video', out_extension='mp4', fps=20, n_loops=4, delete_images=True, delete_old_videos=True, res=res[1])
        
    print('Finished ' + mode)
    