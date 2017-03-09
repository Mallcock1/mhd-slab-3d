# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 16:48:53 2017

@author: Matt

Move seed points by displacement field
"""

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
        seed[1] = seed[1] + disp_z[seed[0],seed[1],seed[2]] * (disp_z.shape[0] / (maxes[2]-mins[2]))
        seed[2] = seed[2] + disp_y[seed[0],seed[1],seed[2]] * (disp_y.shape[0] / (maxes[1]-mins[1]))
        moved_seeds.append(seed)
        
    return moved_seeds