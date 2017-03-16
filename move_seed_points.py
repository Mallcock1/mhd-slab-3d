# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 16:48:53 2017

@author: Matt

Move seed points by displacement field
"""
from scipy.optimize import fsolve


def move_seeds(seeds, disp_x, disp_y, disp_z, mins, maxes):
    """
    Move seed points by displacement field
    
    Parameters
    ----------
    seeds: list of tuples
        original seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: 3d array
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
        
    """
    
    moved_seeds = []
    for seed in seeds:
        seed[0] = seed[0] + disp_x[seed[0],seed[1],seed[2]] * (disp_x.shape[0] / (maxes[0]-mins[0]))
        seed[1] = seed[1] + disp_z[seed[0],seed[1],seed[2]] * (disp_z.shape[2] / (maxes[2]-mins[2]))
        seed[2] = seed[2] + disp_y[seed[0],seed[1],seed[2]] * (disp_y.shape[1] / (maxes[1]-mins[1]))
        moved_seeds.append(seed)
        
    return moved_seeds
    
def original_seeds(moved_seeds, disp_x, disp_y, disp_z, mins, maxes):
    """
    Find original seed points
    
    Parameters
    ----------
    moved_seeds: list of tuples
        moved seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: 3d array
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
        
    """
    
#    seed[0] = seed[0] + disp_x[seed[0],seed[1],seed[2]] * (disp_x.shape[0] / (maxes[0]-mins[0]))
#    seed[1] = seed[1] + disp_z[seed[0],seed[1],seed[2]] * (disp_z.shape[0] / (maxes[2]-mins[2]))
#    seed[2] = seed[2] + disp_y[seed[0],seed[1],seed[2]] * (disp_y.shape[0] / (maxes[1]-mins[1]))    
#    
#    def func(r):
#        return [r[0] - x[i] + xix(mode, r[0], r[1], t, W, K, R1), 
#                r[1] - z[j] + xiz(mode, r[0], r[1], t, W, K, R1)]
#    sol = np.real(fsolve(func, [x[i],z[j]], xtol=1e-03))
    
    seeds = []
    for seed in moved_seeds:
        def func(orig_seed):
            return (orig_seed[0] - seed[0] + disp_x[orig_seed] * (disp_x.shape[0] / (maxes[0]-mins[0])),
                    orig_seed[1] - seed[1] + disp_x[orig_seed] * (disp_x.shape[2] / (maxes[2]-mins[2])),
                    orig_seed[2] - seed[2] + disp_x[orig_seed] * (disp_x.shape[1] / (maxes[1]-mins[1])))
        seeds.append(fsolve(func, seed, xtol=1e-03))
    return seeds