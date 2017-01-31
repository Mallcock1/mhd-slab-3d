from __future__ import division
import Toolbox as tool
import numpy as np
import scipy as sc
#import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib
import slab_functions as sf
#import pickle

R1 = 1.5

def disp_rel_asym_1var(W, K):
    return sf.disp_rel_asym(W, K, R1)
# We define W = c_ph / v_A.
# c0 = c_0 / v_A
# cT = c_T / v_A
    
c1=sf.c2*np.sqrt(sf.R2*R1**(-1))

x_range = np.linspace(0, 5, 51)
y_range = np.linspace(sf.cT, sf.vA, 51)

y_array_saus = tool.point_finder_scipy(disp_rel_asym_1var, x_range, y_range, args=(None))
#y_array_kink = tool.point_finder_scipy(disp_rel_kink, x_range, y_range, args=(None))

#y_array_saus = tool.point_finder_mpmath(disp_rel_saus_mp, x_range, y_range, args=(None))
#y_array_kink = tool.point_finder_mpmath(disp_rel_kink_mp, x_range, y_range, args=(None))

#pickle.dump(x_range, open("x_range.p", "wb"))
#pickle.dump(y_array, open("y_array.p", "wb"))
#
#x_range = pickle.load(open("x_range.p", "rb"))
#y_array = pickle.load(open("y_array.p", "rb")) 

font = {'size': 15}
matplotlib.rc('font', **font)

plt.figure(num=None, figsize=(8, 11), dpi=80, facecolor='w', edgecolor='k')
ax = plt.subplot()

#x_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 1., 0.378456, 0.0001, 0.35, 5., (None))
#ax.plot(x_values, root_array, color='g')

#x_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 2., 0.378275, 0.0001, 0.9, 5., (None))
#ax.plot(x_values, root_array, color='g')

#x_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 3., 0.378242, 0.0001, 1.4, 5., (None))
#ax.plot(x_values, root_array, color='g')

x_values, root_array = tool.line_trace_scipy(disp_rel_asym_1var, 3.4, 0.376691, 0.0001, 2., 5., (None))
ax.plot(x_values, root_array, color='g')

ax.plot(x_range, y_array_saus, '.', color = 'b')
#ax.plot(x_range, y_array_kink, '.', color = 'g')

ax.set_ylabel(r'$\omega/k c_0$', fontsize = 30)
ax.set_xlabel(r'$k x_0$', fontsize = 30)

ax.set_xlim(x_range[0], x_range[-1])
ax.set_ylim(y_range[0], y_range[-1])

ax.plot([x_range[0], x_range[-1]], [sf.vA, sf.vA], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$v_A$', xy=(x_range[-1] + 0.03, sf.vA - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [sf.cT, sf.cT], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_T$', xy=(x_range[-1] + 0.03, sf.cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [sf.c0, sf.c0], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_0$', xy=(x_range[-1] + 0.03, sf.c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [sf.c2, sf.c2], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_2$', xy=(x_range[-1] + 0.03, sf.c2 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [c1, c1], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_1$', xy=(x_range[-1] + 0.03, c1 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
