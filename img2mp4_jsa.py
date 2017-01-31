# Given a filepath, function sequentially reads all numpy arrays contained within and spits out a video file saved in
# the same location

# --- Add functionality for different video outputs?

# jsallcock 11/16

# Preamble

import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import matplotlib.colors
import glob

# jsallcock 11/16

def img2mp4(filepath=None, prefix='I0', leading_zeros=4, extension='png', output_name=None, fps=5):
    if filepath == None:
        filepath = '/Users/jsallcock/Documents/physics/phd/code/python_port/data/HL2A_HFS/29657/demod/I0/media/'
    if output_name == None:
        output_name = prefix
    in_fps = fps
    out_fps = fps

    start_frame = 1

    command_string = 'ffmpeg -framerate '+str(in_fps)+' -start_number '+str(start_frame)+' -i \
'+filepath+prefix+'%0'+str(leading_zeros)+'d.'+extension+' \
-c:v libx264 -r '+str(out_fps)+' -pix_fmt yuv420p \
'+filepath+output_name+'.mp4'

    os.system(command_string)

    return