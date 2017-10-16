# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 08:45:35 2017

@author: Matt
"""

import img2vid as i2v

save_directory = 'D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations'

fps = 20
res = (1616, 909)

mode_options = ['slow-kink-surf', 'slow-saus-surf', 'slow-saus-body-3',
                'slow-kink-body-3', 'slow-saus-body-2', 'slow-kink-body-2', 
                'slow-saus-body-1', 'slow-kink-body-1', 'fast-saus-body-1',
                'fast-kink-body-1', 'fast-saus-body-2', 'fast-kink-body-2',
                'fast-saus-body-3', 'fast-kink-body-3', 'fast-kink-surf',
                'fast-saus-surf', 'shear-alfven', 'shear-alfven-broadband']

modes = mode_options[:16]

vis_mod_string = 'mag_mag_fade_vel_front_'
R1 = 2.0
view = 'front-parallel'
#for mode in modes:
#    prefix = 'R1_'+str(R1) + '_' + mode + '_' + vis_mod_string + view + '_'
#    i2v.image2video(filepath=save_directory, prefix=prefix, 
#                    output_name=prefix+'video', out_extension='mp4', 
#                    fps=fps, n_loops=4, delete_images=True, 
#                    delete_old_videos=True, res=res[1])



prefix = 'R1_2.0_slow-kink-surf_mag_mag_fade_front-parallel_'
i2v.image2video(filepath=save_directory, prefix=prefix, 
                output_name=prefix+'video', out_extension='mp4', 
                fps=fps, n_loops=4, delete_images=True, 
                delete_old_videos=True, res=res[1])                 
