from __future__ import division
import Toolbox as tool
import numpy as np
import scipy as sc
import mpmath as mp
import matplotlib.pyplot as plt
import matplotlib
#import pickle

def m0(W):
    return sc.sqrt((1 - W**2) * (c0**2 - W**2) / ((1 + c0**2) *  (cT**2 - W**2)))

def me(W):
    return sc.sqrt((vAe**2 - W**2) * (ce**2 - W**2) / ((vAe**2 + c0**2) * (cTe**2 - W**2)))

def disp_rel_saus(W, K):
    return R * (vAe**2 - W**2) * m0(W) * sc.tanh(m0(W) * K) + (1 - W**2) * me(W)

def disp_rel_kink(W, K):
    return R * (vAe**2 - W**2) * m0(W) * sc.tanh(m0(W) * K)**(-1) + (1 - W**2) * me(W)

def m0_mp(W):
    return mp.sqrt((1 - W**2) * (c0**2 - W**2) / ((1 + c0**2) * (cT**2 - W**2)))

def me_mp(W):
    return mp.sqrt(1 - W**2 / ce**2)

def disp_rel_saus_mp(W, K):
    return (1 - W**2) * me_mp(W) - R * W**2 * m0_mp(W) * mp.tanh(m0_mp(W) * K)
    
def disp_rel_kink_mp(W, K):
    return (1 - W**2) * me_mp(W) - R * W**2 * m0_mp(W) * mp.coth(m0_mp(W) * K)

# We define W = c_ph / v_A.
# c0 = c_0 / v_A
# cT = c_T / v_A

c0 = 3
ce = 5
vAe = 2.5
R = c0**2 / ce**2 + 5/6 * ce**(-2)
cT = sc.sqrt(c0**2 / (c0**2 + 1))
cTe = sc.sqrt((ce**2 * vAe**2) / (ce**2 + vAe**2))

x_range = np.linspace(0, 5, 11)
y_range = np.linspace(0, 5, 11)

y_array_saus = tool.point_finder_scipy(disp_rel_saus, x_range, y_range, args=(None))
y_array_kink = tool.point_finder_scipy(disp_rel_kink, x_range, y_range, args=(None))

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

ax.plot(x_range, y_array_saus, '.', color = 'b')
ax.plot(x_range, y_array_kink, '.', color = 'g')

ax.plot(x_range, np.imag(y_array_saus), '.', color = 'r')
ax.plot(x_range, np.imag(y_array_kink), '.', color = 'r')

ax.set_ylabel(r'$\omega/k v_A$', fontsize = 30)
ax.set_xlabel(r'$k x_0$', fontsize = 30)

ax.set_xlim(x_range[0], x_range[-1])
ax.set_ylim(y_range[0], y_range[-1])

ax.plot([x_range[0], x_range[-1]], [cT, cT], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_T$', xy=(x_range[-1] + 0.03, cT - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [c0, c0], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_0$', xy=(x_range[-1] + 0.03, c0 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [ce, ce], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$c_e$', xy=(x_range[-1] + 0.03, ce - 0.01), xycoords='data', annotation_clip=False, fontsize=20)
ax.plot([x_range[0], x_range[-1]], [1, 1], color = '0.5', linestyle='--', linewidth=2)
ax.annotate(r'$v_A$', xy=(x_range[-1] + 0.03, 1 - 0.01), xycoords='data', annotation_clip=False, fontsize=20)