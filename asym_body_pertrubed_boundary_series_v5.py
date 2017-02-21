
#import sys
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\ffmpeg')
##sys.path.append(u'W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/')

import pdb # pause code for debugging at pdb.set_trace()
import gc

import numpy as np

import Toolbox as tool
import slab_functions_perturbed_boundary_v3 as sf

from pysac.plot.mayavi_seed_streamlines import SeedStreamline

import matplotlib.pyplot as plt
import matplotlib

from mayavi import mlab
#mlab.options.offscreen = True

from mayavi.tools.engine_manager import engine_manager
from mayavi.core.registry import registry



import img2vid as i2v

###############################################################################


# What mode do you want? OPTIONS:
mode_options = ['slow-kink-surf', 'slow-saus-surf', 'slow-saus-body-3',
                'slow-kink-body-3', 'slow-saus-body-2', 'slow-kink-body-2', 
                'slow-saus-body-1', 'slow-kink-body-1', 'fast-saus-body-1',
                'fast-kink-body-1', 'fast-saus-body-2', 'fast-kink-body-2',
                'fast-saus-body-3', 'fast-kink-body-3', 'fast-kink-surf',
                'fast-saus-surf']
                
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
fast_body_1_mode_options = ['fast-kink-body-1', 'fast-saus-body-1']
fast_body_2_mode_options = ['fast-kink-body-2', 'fast-saus-body-2']
fast_body_3_mode_options = ['fast-kink-body-3', 'fast-saus-body-3']




# Which angle shall we view from?
view_options = ['front', 'front-parallel', 'top', 'top-parallel' 'front-top',
                'front-side']
#view = 'front'
#view = 'front-parallel'
#view = 'top'
#view = 'top-parallel'
view = 'front-top'
#view = 'front-side'

# Uniform lighting?
#uniform_light = True
uniform_light = False

show_density = False
show_density_lagrang = False
show_density_pert = False
show_density_pert2 = False
show_pressure = False
show_mag = False
show_mag_scale = False
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
show_animation = False
#show_animation = True

# Uncomment the parametrer you would like to see
show_density = True
#show_density_pert = True
show_mag = True
#show_mag_scale = True #must also have show_mag = True
#show_mag_vec = True
show_vel_front = True
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


#for mode_ind in range(14): # for all others. REMEMBER SBB pparameters
#for mode_ind in [14,15]: #for fast body surf. REMEMBER SBS parameters
for mode_ind in [15]: #for an individual mode

    # choose your mode: (note that fast surface modes, i.e. 14 and 15, can only be 
    # found with SBS parameters in slab_functions...)
    mode = mode_options[mode_ind] #All working, 0-15    
    
    # Specify oscillation parameters
    if mode in slow_surf_mode_options:
        K = 2.
    elif mode in slow_body_1_mode_options + slow_body_2_mode_options + slow_body_3_mode_options:
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
        print('mode not found')
            
#    R1 = 1.5 #1.8 # Higher denisty on left than right
    R1 = 2. # Symmetric slab
        
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
    else:
        W = Woptions_slow_body[mode_ind-2]
            
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
        # Plot the dispersion diagram with the chosen mode hghlighted
        
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
            for i in range(len(mode_options)):
                if mode_options[i] == mode:
                    colours.append('g')
                elif mode_options[i] in slow_surf_mode_options + fast_surf_mode_options:
                    colours.append('r')
                else:
                    colours.append('b')
            return colours
    
        
        ##Plot the dots
    #    W_array = tool.point_finder_scipy(disp_rel_asym_2var, K_range, W_range,
    #                                      args=(None))
    #    ax.plot(K_range, W_array, '.', color = 'b')
        
        if R1 == 1.8:
            if mode in fast_surf_mode_options:
                ## Line trace each solutions
                lw = 1.5
                ##Slow surface modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.473,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colous(mode)[0], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.547,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colous(mode)[1], linewidth=lw)
                
                ##Slow body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.758137, 
                                                             0.001, 0.25, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 3., 0.725553, 
                                                             0.001, 0.6, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.747275, 
                                                             0.001, 0.9, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.722607, 
                                                             0.001, 1.2, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 8.2, 0.728825, 
                                                             0.001, 1.6, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                #Fast body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 1.07888, 
                                                             0.01, 0.25, Kmax, (None))
                ax.plot(K_values, root_array, color='b', linewidth=lw)
                
        #        K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.143, 
        #                                                     0.01, 4.3, Kmax, (None))
        #        ax.plot(K_values, root_array, color='b', linewidth=lw)
                
        #        K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 9., 1.18647, 
        #                                                     0.01, 8.3, Kmax, (None))
        #        ax.plot(K_values, root_array, color='b', linewidth=lw)
                
                ##plot some circles at start of fast body modes
                ax.plot(0.25, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b',
                        markeredgewidth=lw)
                ax.plot(4.3, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b',
                        markeredgewidth=lw)
                ax.plot(8.3, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b', 
                        markeredgewidth=lw)        
                    
                ##annotate each line
                #s_labels_x = [1., 3.5, 6.]
                #s_labels_y = [0.76, 0.77, 0.76]
                #for i in range(len(s_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                #f_labels_x = [1., 6., 9.]
                #f_labels_y = [1.05, 1.09, 1.14]
                #for i in range(len(f_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                
                ax.set_ylabel(r'$v_{ph}/c_0$', rotation=0, fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
                #
                #K_range_for_fill = np.append(K_range, K_range[-1]+0.1)
                #ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c0,sf.c0), (sf.vA,sf.vA),
                #                edgecolor='gray', linestyle='--', color='None', hatch='/',
                #                linewidth=2)
                #ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.41, 1.41),
                #                edgecolor='gray', linestyle='--', color='None', hatch='/',
                #                linewidth=2)
                
                #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                #ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
                ##Annotate some things
                #ax.annotate('Leaky modes', xy=(xmax/2. - 1, (sf.vA + sf.c0)/2. - 0.01),
                #            xycoords='data', annotation_clip=False, fontsize=15)
                #ax.annotate('Leaky modes', xy=(xmax/2. - 1, (sf.c2 + 1.41) / 2. - 0.01),
                #            xycoords='data', annotation_clip=False, fontsize=15)
            
            else:
                ## Line trace each solutions
                #lw = 1.5
                ##Slow surface modes
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 1., 0.473,
                #                                             0.01, 0.0001, 10., (None))
                #ax.plot(K_values, root_array, color='r', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 1., 0.547,
                #                                             0.01, 0.0001, 10., (None))
                #ax.plot(K_values, root_array, color='r', linewidth=lw)
                #
                ##Slow body modes
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 2., 0.758137, 
                #                                             0.001, 0.25, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 3., 0.725553, 
                #                                             0.001, 0.6, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 6., 0.747275, 
                #                                             0.001, 0.9, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 6., 0.722607, 
                #                                             0.001, 1.2, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 8.2, 0.728825, 
                #                                             0.001, 1.6, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                ##Fast body modes
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 2., 1.07888, 
                #                                             0.01, 0.25, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 6., 1.143, 
                #                                             0.01, 4.3, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                #K_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 9., 1.18647, 
                #                                             0.01, 8.3, 10., (None))
                #ax.plot(K_values, root_array, color='b', linewidth=lw)
                #
                ###plot some circles at start of fast body modes
                #ax.plot(0.25, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b',
                #        markeredgewidth=lw)
                #ax.plot(4.3, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b',
                #        markeredgewidth=lw)
                #ax.plot(8.3, sf.c2, marker='o', markerfacecolor='None', markeredgecolor='b', 
                #        markeredgewidth=lw)        
                    
                #annotate each line
                s_labels_x = [1., 3.5, 6.]
                s_labels_y = [0.76, 0.77, 0.76]
                for i in range(len(s_labels_x)):
                    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                                xycoords='data', annotation_clip=False, fontsize=15)
                f_labels_x = [1., 6., 9.]
                f_labels_y = [1.05, 1.09, 1.14]
                for i in range(len(f_labels_x)):
                    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                                xycoords='data', annotation_clip=False, fontsize=15)
                
                ax.set_ylabel(r'$v_{ph}/c_0$', rotation=0, fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
                
                K_range_for_fill = np.append(K_range, K_range[-1]+0.1)
                ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c0,sf.c0), (sf.vA,sf.vA),
                                edgecolor='gray', linestyle='--', color='None', hatch='/',
                                linewidth=2)
                ax.fill_between((Kmin-0.1, Kmax+0.1), (sf.c2, sf.c2), (1.41, 1.41),
                                edgecolor='gray', linestyle='--', color='None', hatch='/',
                                linewidth=2)
            
                #ax.plot([K_range[0], K_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$v_A$', xy=(K_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_T$', xy=(K_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                #ax.plot([K_range[0], K_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_0$', xy=(K_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                #ax.plot([K_range[0], K_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='--', linewidth=2)
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
                #Annotate some things
                ax.annotate('Leaky modes', xy=(Kmax/2. - 1, (sf.vA + sf.c0)/2. - 0.01),
                            xycoords='data', annotation_clip=False, fontsize=15)
                ax.annotate('Leaky modes', xy=(Kmax/2. - 1, (sf.c2 + 1.41) / 2. - 0.01),
                            xycoords='data', annotation_clip=False, fontsize=15)
        
        elif R1 == 2.:
            if mode in fast_surf_mode_options:
                # Line trace each solutions
                lw = 1.5
    #            ##Slow surface modes
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.64,
    #                                                         0.01, 0.0001, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[0], linewidth=lw, linestyle='--')
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.72,
    #                                                         0.01, 0.0001, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[1], linewidth=lw)
    #            
    #            ##Slow body modes
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.9301, 
    #                                                         0.001, 0.2, 5.24, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[7], linewidth=lw, linestyle='--')
    #            #
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.831416, 
    #                                                         0.001, 0.6, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[6], linewidth=lw)
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 4., 0.848421, 
    #                                                         0.001, 0.9, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[5], linewidth=lw, linestyle='--')
    #            #
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 4., 0.826689, 
    #                                                         0.001, 1.2, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[4], linewidth=lw)
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.837169, 
    #                                                         0.001, 1.6, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[3], linewidth=lw, linestyle='--')
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 10.56, 0.866878, 
    #                                                         0.001, 2., Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[2], linewidth=lw)
    #            
    #            
    #            
    #            ##Fast surface modes
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.01567, 
    #                                                         0.01, 5.24, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[14], linewidth=lw, linestyle='--')
    #            #
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.05754, 
    #                                                         0.01, 0.001, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[15], linewidth=lw)
    #            
                
                ##plot some circles at start of fast surf modes
                
    #            ax.plot(0.6, sf.c2, marker='o', markerfacecolor='None', markeredgecolor=colour(mode)[15],
    #                    markeredgewidth=lw, markersize=10)
                ax.plot(5.24, sf.c0, marker='o', markerfacecolor='None', markeredgecolor=colour(mode)[14],
                        markeredgewidth=lw, markersize=10)    
                    
                ##annotate each line
                #s_labels_x = [1., 3.5, 6.]
                #s_labels_y = [0.76, 0.77, 0.76]
                #for i in range(len(s_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                #f_labels_x = [1., 6., 9.]
                #f_labels_y = [1.05, 1.09, 1.14]
                #for i in range(len(f_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                
                ax.set_ylabel(r'$v_{ph}/c_0$', rotation=0, fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
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
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
            
            else:
                # Line trace each solutions
                lw = 1.5
                #Slow surface modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.473,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[0], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.547,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[1], linewidth=lw)
                
                #Slow body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1.84, 0.749326, 
                                                             0.001, 0.3, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[7], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 3.68, 0.741942, 
                                                             0.001, 0.8, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[6], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 5.98, 0.746659, 
                                                             0.001, 1., Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[5], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 5.98, 0.722183, 
                                                             0.001, 1.3, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[4], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 8.74, 0.733764, 
                                                             0.001, 1.8, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[3], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 8.74, 0.719463, 
                                                             0.001, 2.3, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[2], linewidth=lw)
            
                
                #Fast body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 1.07888, 
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[8], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.143, 
                                                             0.01, 4.1, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[9], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 9., 1.18647, 
                                                             0.01, 8.2, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[10], linewidth=lw)
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 15., 1.15692, 
                                                             0.01, 12.2, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[11], linewidth=lw, linestyle='--')
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 19.8, 1.156, 
                                                             0.01, 16.2, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[12], linewidth=lw)
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 23.4, 1.16675, 
                                                             0.01, 20.2, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[13], linewidth=lw, linestyle='--')
                
                ##plot some circles at start of fast body modes
                K_marker = [4.1, 8.2, 12.2, 16.2, 20.2]
                y_marker = [sf.c2] * len(K_marker)
                for i in range(len(K_marker)):
                    ax.plot(K_marker[i], y_marker[i], marker='o', markersize=10, 
                            markerfacecolor='None', markeredgecolor=colour(mode)[9 + i], markeredgewidth=lw)
                        
                ##annotate each line
                #s_labels_x = [1., 3.5, 6.]
                #s_labels_y = [0.76, 0.77, 0.76]
                #for i in range(len(s_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                #f_labels_x = [1., 6., 9.]
                #f_labels_y = [1.05, 1.09, 1.14]
                #for i in range(len(f_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                
                #ax.set_ylabel(r'$v_{ph}$', rotation=0, fontsize = 20)
                ax.set_ylabel(r'$v_{ph}/c_0$', fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                #ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
                
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
                ax.annotate(r'$c_1=c_2$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
    
        
        elif R1 == 1.5:
            if mode in fast_surf_mode_options:
                # Line trace each solutions
                lw = 1.5
                ##Slow surface modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.64,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[0], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.72,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[1], linewidth=lw)
                
                ##Slow body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.9301, 
                                                             0.001, 0.2, 4.199, (None))
                ax.plot(K_values, root_array, color=colour(mode)[7], linewidth=lw, linestyle='--')
                #
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.831416, 
                                                             0.001, 0.6, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[6], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 4., 0.848421, 
                                                             0.001, 0.9, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[5], linewidth=lw, linestyle='--')
                #
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 4., 0.826689, 
                                                             0.001, 1.2, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[4], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.837169, 
                                                             0.001, 1.6, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[3], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 10.56, 0.866878, 
                                                             0.001, 2., Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[2], linewidth=lw)
                
                #
                #
                ##Fast surface modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.01567, 
                                                             0.01, 4.19, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[14], linewidth=lw, linestyle='--')
                #
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.05754, 
                                                             0.01, 0.6, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[15], linewidth=lw)
                
                
                ##plot some circles at start of fast surf modes
                
                ax.plot(0.6, sf.c2, marker='o', markerfacecolor='None', markeredgecolor=colour(mode)[15],
                        markeredgewidth=lw, markersize=10)
                ax.plot(4.19, sf.c0, marker='o', markerfacecolor='None', markeredgecolor=colour(mode)[14],
                        markeredgewidth=lw, markersize=10)    
                    
                ##annotate each line
                #s_labels_x = [1., 3.5, 6.]
                #s_labels_y = [0.76, 0.77, 0.76]
                #for i in range(len(s_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                #f_labels_x = [1., 6., 9.]
                #f_labels_y = [1.05, 1.09, 1.14]
                #for i in range(len(f_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                
                ax.set_ylabel(r'$v_{ph}/c_0$', rotation=0, fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
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
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                
            
            else:
                # Line trace each solutions
                lw = 1.5
                #Slow surface modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.473,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[0], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 1., 0.547,
                                                             0.01, 0.0001, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[1], linewidth=lw)
    #            
    #            #Slow body modes
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 0.760531, 
    #                                                         0.001, 0.2, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[7], linewidth=lw, linestyle='--')
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 3., 0.726328, 
    #                                                         0.001, 0.5, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[6], linewidth=lw)
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.747901, 
    #                                                         0.001, 0.8, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[5], linewidth=lw, linestyle='--')
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 0.722958, 
    #                                                         0.001, 1.1, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[4], linewidth=lw)
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 8.2, 0.729129, 
    #                                                         0.001, 1.5, Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[3], linewidth=lw, linestyle='--')
    #            
    #            K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 9., 0.721913, 
    #                                                         0.001, 2., Kmax, (None))
    #            ax.plot(K_values, root_array, color=colour(mode)[2], linewidth=lw)
            
                
                #Fast body modes
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 2., 1.07888, 
                                                             0.01, 0.45, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[8], linewidth=lw)
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 6., 1.143, 
                                                             0.01, 4.5, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[9], linewidth=lw, linestyle='--')
                
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 9., 1.18647, 
                                                             0.01, 8.45, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[10], linewidth=lw)
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 15., 1.15692, 
                                                             0.01, 12.49, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[11], linewidth=lw, linestyle='--')
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 19.8, 1.156, 
                                                             0.01, 16.5673, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[12], linewidth=lw)
                K_values, root_array = tool.line_trace_scipy(disp_rel_asym_2var, 23.4, 1.16675, 
                                                             0.01, 20.6113, Kmax, (None))
                ax.plot(K_values, root_array, color=colour(mode)[13], linewidth=lw, linestyle='--')
                
                ##plot some circles at start of fast body modes
                K_marker = [0.45, 4.5, 8.45, 12.49, 16.5673, 20.6113]
                y_marker = [sf.c2] * len(K_marker)
                for i in range(len(K_marker)):
                    ax.plot(K_marker[i], y_marker[i], marker='o', markersize=10, 
                            markerfacecolor='None', markeredgecolor=colour(mode)[8 + i], markeredgewidth=lw)
                        
                ##annotate each line
                #s_labels_x = [1., 3.5, 6.]
                #s_labels_y = [0.76, 0.77, 0.76]
                #for i in range(len(s_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(s_labels_x[i], s_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                #f_labels_x = [1., 6., 9.]
                #f_labels_y = [1.05, 1.09, 1.14]
                #for i in range(len(f_labels_x)):
                #    ax.annotate(r'$j={}$'.format(i+1), xy=(f_labels_x[i], f_labels_y[i]),
                #                xycoords='data', annotation_clip=False, fontsize=15)
                
                #ax.set_ylabel(r'$v_{ph}$', rotation=0, fontsize = 20)
                ax.set_ylabel(r'$v_{ph}/c_0$', fontsize = 20)
                ax.set_xlabel(r'$k x_0$', fontsize = 20)
                #ax.yaxis.labelpad = 17
                
                ax.set_xlim(K_range[0], K_range[-1])
                ax.set_ylim(0., 1.41)
                
                
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
                ax.annotate(r'$c_2$', xy=(K_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
                ax.plot([K_range[0], K_range[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='-.', linewidth=2)
                ax.annotate(r'$c_1$', xy=(K_range[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
    
        ax.plot(K, W, 'go', markersize=10)
#        plt.tight_layout() # seems to make it chop the sides off with this
        plt.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_dispersion_diagrams\\'
                    + mode + '.png')
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
        nx = 100
        ny = 100#20 #100
        nz = 100
        nt = number_of_frames
        
        if nz % nt != 0:
            print("nt doesnt divide nz so there may be a problem")
        
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
        mod = 4.
        mod_top = np.ceil(4. / y_spacing)
        
        if show_disp_top == True or show_disp_front == True:
            xixvals = np.real(np.repeat(sf.xix(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xizvals = np.real(np.repeat(sf.xiz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            xiyvals = np.zeros_like(xixvals)
        
                                
        if show_vel_front == True or show_vel_top == True:
            vxvals = np.real(np.repeat(sf.vx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.zeros_like(vxvals)
        
            
        if show_vel_front_pert == True or show_vel_top_pert == True:
            vxvals = np.real(np.repeat(sf.vx_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vzvals = np.real(np.repeat(sf.vz_pert(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
            vyvals = np.zeros_like(vxvals)
        
        bxvals = np.real(np.repeat(sf.bx(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
        bz_eq3d = np.repeat(sf.bz_eq(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2)
        bzvals = np.real(np.repeat(-sf.bz(mode, xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2) +
                         bz_eq3d)
                         
                         
        byvals = np.zeros_like(bxvals)
        
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
                bzvals_split = np.split(bzvals, [nz - (nz / nt) * t_ind], axis=1)
                
                bxvals_t = np.concatenate((bxvals_split[1], bxvals_split[0]), axis=1)
                byvals_t = byvals
                bzvals_t = np.concatenate((bzvals_split[1], bzvals_split[0]), axis=1)
                
                
                if show_disp_top == True or show_disp_front == True:            
                    xixvals_split = np.split(xixvals, [nz - (nz / nt) * t_ind], axis=1)
                    xizvals_split = np.split(xizvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    xixvals_t = np.concatenate((xixvals_split[1], xixvals_split[0]), axis=1)
                    xiyvals_t = xiyvals
                    xizvals_t = np.concatenate((xizvals_split[1], xizvals_split[0]), axis=1)
                                        
                if np.array([show_vel_top, show_vel_top_pert, show_vel_front, show_vel_front_pert]).any() == True:            
                    vxvals_split = np.split(vxvals, [nz - (nz / nt) * t_ind], axis=1)
                    vzvals_split = np.split(vzvals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    vxvals_t = np.concatenate((vxvals_split[1], vxvals_split[0]), axis=1)
                    vyvals_t = vyvals
                    vzvals_t = np.concatenate((vzvals_split[1], vzvals_split[0]), axis=1)
                
                if show_boundary == True:
                    xix_boundary_r_vals_split = np.split(xix_boundary_r_vals, [nz - (nz / nt) * t_ind], axis=0)
                    xix_boundary_l_vals_split = np.split(xix_boundary_l_vals, [nz - (nz / nt) * t_ind], axis=0)
        
                    xix_boundary_r_vals_t = np.concatenate((xix_boundary_r_vals_split[1], xix_boundary_r_vals_split[0]), axis=0)
                    xix_boundary_l_vals_t = np.concatenate((xix_boundary_l_vals_split[1], xix_boundary_l_vals_split[0]), axis=0)
                
                if show_density == True or show_density_pert == True:            
                    rho_vals_split = np.split(rho_vals, [nz - (nz / nt) * t_ind], axis=1)
                    
                    rho_vals_t = np.concatenate((rho_vals_split[1], rho_vals_split[0]), axis=1)                         
                
                
            bxvals_mask_front_t = np.copy(bxvals_t)
            byvals_mask_front_t = np.copy(byvals_t)
            bzvals_mask_front_t = np.copy(bzvals_t)
            
            for i in range(bxvals_t.shape[0]):
                for j in range(bxvals_t.shape[1]):
                    for k in range(bxvals_t.shape[2]):
                        if (i%mod) != 0 or (j%mod) != 0:
                            bxvals_mask_front_t[i,j,k] = 0.
                            bzvals_mask_front_t[i,j,k] = 0.
        
        
            if show_disp_top == True:    
                xixvals_mask_top_t = np.copy(xixvals_t)
                xiyvals_mask_top_t = np.copy(xiyvals_t)
                xizvals_mask_top_t = np.copy(xizvals_t)
                
                for i in range(xixvals_t.shape[0]):
                    for j in range(xixvals_t.shape[1]):
                        for k in range(xixvals_t.shape[2]):
                            if (i%mod) != 0 or (k%mod_top) != 0:
                                xixvals_mask_top_t[i,j,k] = 0.
                                xizvals_mask_top_t[i,j,k] = 0.
            if show_disp_front == True:
                xixvals_mask_front_t = np.copy(xixvals_t)
                xiyvals_mask_front_t = np.copy(xiyvals_t)
                xizvals_mask_front_t = np.copy(xizvals_t)
                
                for i in range(xixvals_t.shape[0]):
                    for j in range(xixvals_t.shape[1]):
                        for k in range(xixvals_t.shape[2]):
                            if (i%mod)!=0 or (j%mod)!=0:
                                xixvals_mask_front_t[i,j,k] = 0.
                                xizvals_mask_front_t[i,j,k] = 0.    
            
            
            if show_vel_top == True or show_vel_top_pert == True:    
                vxvals_mask_top_t = np.copy(vxvals_t)
                vyvals_mask_top_t = np.copy(vyvals_t)
                vzvals_mask_top_t = np.copy(vzvals_t)
                
                for i in range(vxvals_t.shape[0]):
                    for j in range(vxvals_t.shape[1]):
                        for k in range(vxvals_t.shape[2]):
                            if (i%mod) != 0 or (k%mod_top) != 0:
                                vxvals_mask_top_t[i,j,k] = 0
                                vzvals_mask_top_t[i,j,k] = 0
                                                  
                                
            if show_vel_front == True or show_vel_front_pert == True:
                vxvals_mask_front_t = np.copy(vxvals_t)
                vyvals_mask_front_t = np.copy(vyvals_t)
                vzvals_mask_front_t = np.copy(vzvals_t)
                
                for i in range(vxvals_t.shape[0]):
                    for j in range(vxvals_t.shape[1]):
                        for k in range(vxvals_t.shape[2]):
                            if (i%mod) != 0 or (j%mod) != 0:
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
                    boundary_r_thick = mlab.mesh(xix_boundary_r_vals_t, zvals, yvals,
                                                 extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                                 color=(1.,1.,1.), opacity=1., representation='wireframe',
                                                 line_width=12.)
                
                    boundary_l_thick = mlab.mesh(xix_boundary_l_vals_t, zvals, yvals,
                                                 extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                                 color=(1.,1.,1.), opacity=1., representation='wireframe',
                                                 line_width=12.)
                                                 
                else:
                    boundary_r = mlab.mesh(xix_boundary_r_vals_t, zvals, yvals,
                                           extent=[ext_min_r, ext_max_r, 1, nz, 0, (ny-1) * y_spacing],
                                           color=(0.5,0.5,0.5), opacity=0.7)
                    
                    boundary_l = mlab.mesh(xix_boundary_l_vals_t, zvals, yvals,
                                           extent=[ext_min_l, ext_max_l, 1, nz, 0, (ny-1) * y_spacing],
                                           color=(0.5,0.5,0.5), opacity=0.7)
    #        if show_density == True or show_density_pert == True:
    #            # Scalar field density   
    #            rho = mlab.pipeline.scalar_field(rho_vals_t, name="density", figure=fig)
    #            #scalar_cut_plane = ScalarCutPlane()
    #            #fig.parent.add_filter(scalar_cut_plane, sca)
    #            rho.spacing = spacing
    #            minr = rho_vals_t.min()
    #            maxr = rho_vals_t.max()
    #            
    #            #Volume for high pressure
    #            rvmin1 = minr + 0.55 * (maxr - minr)
    #            rvmax1 = minr + 1. * (maxr - minr)
    #            rvol1 = mlab.pipeline.volume(rho, vmin=rvmin1, vmax=rvmax1)
    #            
    #            # Changing the ctf:
    #            from tvtk.util.ctf import ColorTransferFunction
    #            ctf1 = ColorTransferFunction()
    #            ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    #            ctf1.add_rgb_point(rvmin1 + 0.5 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    #            ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    #            # ...
    #            rvol1._volume_property.set_color(ctf1)
    #            rvol1._ctf = ctf1
    #            rvol1.update_ctf = True
    #            
    #            #Changing the opacity of the volume vol1
    #            ## Changing the otf:
    #            from tvtk.util.ctf import PiecewiseFunction
    #            otf = PiecewiseFunction()
    #            otf.add_point(rvmin1, 0)
    #            otf.add_point(rvmax1, 0.10)
    #            ##vol1._otf = otf
    #            rvol1._volume_property.set_scalar_opacity(otf)
    #            
    #            # exempt volume from shading and improve overall look by increasing opacity
    #            rvol1.volume_property.shade = False
    #            rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    #            
    #            
    #            #Volume for low pressure
    #            rvmin2 = minr + 0. * (maxr - minr)
    #            rvmax2 = minr + 0.45 * (maxr - minr)
    #            rvol2 = mlab.pipeline.volume(rho, vmin=rvmin2, vmax=rvmax2)
    #            
    #            # Changing the ctf:
    #            ctf2 = ColorTransferFunction()
    #            ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    #            ctf2.add_rgb_point(rvmin2 + 0.5 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    #            ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    #            # ...
    #            rvol2._volume_property.set_color(ctf2)
    #            rvol2._ctf = ctf2
    #            rvol2.update_ctf = True
    #        
    #            #Changing the opacity of the volume vol1
    #            ## Changing the otf:
    #            otf = PiecewiseFunction()
    #            otf.add_point(rvmax2, 0)
    #            otf.add_point(rvmin2, 0.10)
    #            ##vol1._otf = otf
    #            rvol2._volume_property.set_scalar_opacity(otf)
    #            
    #            # exempt volume from shading and improve overall look by increasing opacity
    #            rvol2.volume_property.shade = False
    #            rvol2.volume_property.scalar_opacity_unit_distance = 2.0
                                                 
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
            field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
                                                   figure=fig)
            field.spacing = spacing
                
            if show_mag == True:
                magnitude = mlab.pipeline.extract_vector_norm(field)
                #contours = mlab.pipeline.iso_surface(magnitude,
                #                                        contours=range(2, 14, 3),
                #                                        transparent=True,
                #                                        opacity=0.4,
                #                                        colormap='YlGnBu',
                #                                        vmin=0, vmax=14)
                
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
                for i in range(0,nx_seed):
                    for j in range(0,ny_seed):
                        x = start_x + (i * dx_res) * x_spacing
                        y = start_y + (j * dy_res) * y_spacing
                        z = 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                        seeds.append((x,z,y))
                                                             
                field_lines = SeedStreamline(seed_points=seeds)
                field_lines.stream_tracer.integration_direction='both'
                field_lines.streamline_type = 'tube'
                
                magnitude = mlab.pipeline.extract_vector_norm(field)
                magnitude.add_child(field_lines)
                module_manager = field_lines.parent
                module_manager.scalar_lut_manager.lut_mode = 'Reds'
                module_manager.scalar_lut_manager.data_range=[-30,25]
                field_lines.stream_tracer.maximum_propagation = 500.
                field_lines.tube_filter.number_of_sides = 20
                field_lines.tube_filter.radius = 0.7
                field_lines.tube_filter.capping = True
                
                if show_mag_scale == True:
                    module_manager.scalar_lut_manager.lut_mode = 'jet'
                    module_manager.scalar_lut_manager.data_range=[7,18]
                
            
            if show_mag_vec == True:
                bdirfield_front = mlab.pipeline.vector_field(bxvals_mask_front_t, bzvals_mask_front_t,
                                                             byvals_mask_front_t, name="B field front",
                                                             figure=fig)
                bdirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(bdirfield_front, 
                                                                  scale_factor=8.)
                vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
                vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
                vector_cut_plane_front.glyph.color_mode = 'no_coloring'
                vector_cut_plane_front.implicit_plane.widget.enabled = False
                vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
            
            
            if show_vel_top == True or show_vel_top_pert == True:
                vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top_t, vyvals_mask_top_t,
                                                            vyvals_mask_top_t, name="V field top",
                                                            figure=fig)
                vdirfield_top.spacing = spacing
                vector_cut_plane_top = mlab.pipeline.vector_cut_plane(vdirfield_top, 
                                                                  scale_factor=8.)
                vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
                vector_cut_plane_top.glyph.color_mode = 'no_coloring'
                vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-1, 50.5])
                vector_cut_plane_top.implicit_plane.widget.enabled = False
                vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
                
            if show_vel_front == True or show_vel_front_pert == True:
                vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front_t, vzvals_mask_front_t,
                                                             vyvals_mask_front_t, name="V field front",
                                                             figure=fig)
                vdirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                                  scale_factor=8.)
                vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
                vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
                vector_cut_plane_front.glyph.color_mode = 'no_coloring'
                vector_cut_plane_front.implicit_plane.widget.enabled = False
                vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
            
            if show_disp_top == True:
                xidirfield_top = mlab.pipeline.vector_field(xixvals_mask_top_t, xiyvals_mask_top_t,
                                                            xiyvals_mask_top_t, name="Xi field top",
                                                            figure=fig)
                xidirfield_top.spacing = spacing
                vector_cut_plane_top = mlab.pipeline.vector_cut_plane(xidirfield_top, 
                                                                  scale_factor=8.)
                vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
                vector_cut_plane_top.glyph.color_mode = 'no_coloring'
                vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-1, 50.5])
                vector_cut_plane_top.implicit_plane.widget.enabled = False
                vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
                vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
                
            if show_disp_front == True:
                xidirfield_front = mlab.pipeline.vector_field(xixvals_mask_front_t, xizvals_mask_front_t,
                                                             xiyvals_mask_front_t, name="Xi field front",
                                                             figure=fig)
                xidirfield_front.spacing = spacing
                vector_cut_plane_front = mlab.pipeline.vector_cut_plane(xidirfield_front, 
                                                                  scale_factor=8.)
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
                field.scene.camera.position = [126.62380038836577, 60.458994158678557, 524.80329455823346]
                field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
                field.scene.camera.view_angle = 14.
                field.scene.camera.view_up = [-0.016952270974394879, 0.9996860427028168, -0.018450922307365961]
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
    
            prefix = 'R1_'+str(R1)+'_'+view + '_' + mode
            
            mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
                         + prefix + str(t_ind+1) + '.png')
            
    # Trying and failing to sort out memory issues.
    #        mlab.gcf()
    #        mlab.clf()
    #        mlab.close()
    #        gc.collect()
    #        del fig
    #        engine_manager.current_engine = None
    #        registry.engines = {}
            
            if make_video == True:
                mlab.close(fig)
    
            t = t + (t_end - t_start) / nt
        if make_video == True:
            i2v.image2video(prefix=prefix, output_name=prefix+'_video', out_extension='mp4', fps=20, n_loops=4, delete_images=True, delete_old_videos=True, res=res[1])