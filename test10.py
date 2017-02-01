# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 12:18:23 2017

@author: Matt
"""
import sys
sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3d/ vis')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')

import Toolbox as tool
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from scipy.optimize import newton
import solver

import slab_functions_perturbed_boundary as sf


###############################################################################

# Which angle shall we view from?
#view = 'front'
#view = 'front parallel'
#view = 'top'
#view = 'top parallel'
#view = 'front top'
view = 'front side'

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
show_boundary = False

show_density = True
#show_density_lagrang = True
#show_density_pert = True
#show_density_pert2 = True #THIS ONE WORKS
#show_pressure = True
show_mag = True
#show_mag_scale = True
#show_mag_vec = True
show_vel_front = True
#show_vel_front_pert = True
#show_vel_top = True
#show_vel_top_pert = True
#show_disp_top = True
#show_disp_front = True
show_axes = True
show_boundary = True

# Specify oscillation parameters
K = 2.
R1 = 1.8 #1.8 #2.

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)
    
def disp_rel_asym_2var(W, K):
    return sf.disp_rel_asym(W, K, R1)


Wrange = np.linspace(0., sf.c2, 51)
Wvals = tool.point_finder_scipy(disp_rel_asym_2var, np.array(K), Wrange, args=None).transpose()
tol = 1e-2
indices_to_rm = []
for i in range(len(Wvals)):
    w = Wvals[i]
    if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < 0 or w > max(sf.c0, sf.vA, sf.c1, sf.c2):
        indices_to_rm.append(i)
Wvals = np.delete(Wvals, indices_to_rm)


Kvals = np.linspace(0., 3, 51)

plt.plot(K * np.ones_like(Wvals), Wvals, '.', color = 'b')
plt.ylabel(r'$\omega/k c_0$', fontsize = 30)
plt.xlabel(r'$k x_0$', fontsize = 30)

plt.xlim(Kvals[0], Kvals[-1])
plt.ylim(Wrange[0], Wrange[-1])

plt.plot([Kvals[0], Kvals[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
plt.annotate(r'$v_A$', xy=(Kvals[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
plt.plot([Kvals[0], Kvals[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='--', linewidth=2)
plt.annotate(r'$c_T$', xy=(Kvals[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
plt.plot([Kvals[0], Kvals[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
plt.annotate(r'$c_0$', xy=(Kvals[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
plt.plot([Kvals[0], Kvals[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
plt.annotate(r'$c_2$', xy=(Kvals[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
plt.plot([Kvals[0], Kvals[-1]], [sf.c1(R1), sf.c1(R1)], color = '0.5', linestyle='--', linewidth=2)
plt.annotate(r'$c_1$', xy=(Kvals[-1] + 0.03, sf.c1(R1) - 0.01), xycoords='data', annotation_clip=False, fontsize=20)

