# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 16:48:53 2017

@author: Matt

Move seed points by displacement field
"""
from scipy.optimize import fsolve
import slab_functions as sf
import numpy as np

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
        seed = list(seed)
        seed[0] = seed[0] + disp_x[seed[0],seed[1],seed[2]] * (disp_x.shape[0] / (maxes[0]-mins[0]))
        seed[1] = seed[1] + disp_z[seed[0],seed[1],seed[2]] * (disp_z.shape[2] / (maxes[2]-mins[2]))
        seed[2] = seed[2] + disp_y[seed[0],seed[1],seed[2]] * (disp_y.shape[1] / (maxes[1]-mins[1]))
        seed = tuple(seed)        
        moved_seeds.append(seed)
        
    return moved_seeds
    
def move_seeds_non_int(seeds, mins, maxes, n, 
                      mode, x, z, t, W, K, R1):
    """
    Move seed points by displacement field
    
    Parameters
    ----------
    seeds: list of tuples
        original seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: functions
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. [xmin, ymin, zmin]
    """
    
    moved_seeds = []
    for seed in seeds:
        new_seed = list(seed)
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)))
#        print(np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)) * (n[0] / (maxes[0]-mins[0])))
        new_seed[0] = seed[0] + np.real(sf.xix(mode, x[seed[0]], z[seed[1]], t, W, K, R1)) * (n[0] / (maxes[0]-mins[0]))
        new_seed[1] = seed[1] + np.real(sf.xiz(mode, x[seed[0]], z[seed[1]], t, W, K, R1)) * (n[2] / (maxes[2]-mins[2]))
        new_seed[2] = seed[2] + np.real(sf.xiy(mode, x[seed[0]], z[seed[1]], t, W, K))     * (n[1] / (maxes[1]-mins[1]))       
        moved_seeds.append(tuple(new_seed))
        
    return moved_seeds
    
def original_seeds_non_int(moved_seeds, mins, maxes, n, 
                           mode, x, z, t, W, K, R1):
    """
    Find original seed points
    
    Parameters
    ----------
    moved_seeds: list of tuples
        moved seed points, e.g. [(x,y,z),(x,y,z)]
        
    disp_x, disp_y, disp_z: 3d array
        x, y, z displacement
    
    mins, maxes: list
        mins and maxes of x, y, z coords. e.g. [xmin, ymin, zmin]
        
    """
    
    seeds = []
    for seed in moved_seeds:
        seed = list(seed)
        def function(orig_seed):
            return [orig_seed[0] - seed[0] + np.real(sf.xix(mode, x[orig_seed[0]], z[orig_seed[1]], t, W, K, R1)) * (n[0] / (maxes[0]-mins[0])),
                    orig_seed[1] - seed[1] + np.real(sf.xiz(mode, x[orig_seed[0]], z[orig_seed[1]], t, W, K, R1)) * (n[2] / (maxes[2]-mins[2])),
                    orig_seed[2] - seed[2] + np.real(sf.xiy(mode, x[orig_seed[0]], z[orig_seed[1]], t, W, K)) *     (n[1] / (maxes[1]-mins[1]))]
        original_seed = list(np.real(fsolve(function, seed, xtol=1e-03)))
        seeds.append(tuple(original_seed))
    return seeds
    