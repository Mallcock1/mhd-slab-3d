
#import sys
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\ffmpeg')
##sys.path.append(u'W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/')

#import pdb # pause code for debugging at pdb.set_trace()

import numpy as np

import slab_functions_perturbed_boundary_v3 as sf

import move_seed_points as msp

from pysac.plot.mayavi_seed_streamlines import SeedStreamline

from mayavi import mlab
#mlab.options.offscreen = True

import mayavi_plotting_functions as mpf

import img2vid as i2v

###############################################################################

vA2 = 1.
vA1 = 0.9

def vA_func(x):
    return (vA2 - vA1)/4 * x + vA1

max_amplitude = 1.


# Which angle shall we view from?
view_options = ['front', 'front-parallel', 'top', 'top-parallel' 'front-top',
                'front-side', 'front-top-side']
#view = 'front'
#view = 'front-parallel'
#view = 'top'
#view = 'top-parallel'
#view = 'front-top'
#view = 'front-side'
view = 'front-top-side'

# Uniform lighting?
uniform_light = True
#uniform_light = False

show_mag = False
show_mag_scale = False
show_mag_fade = False
show_mag_vec = False
show_vel_top = False
show_disp_top = False
show_axes = False
show_axis_labels = False
show_mini_axis = False

# Uncomment the parametrer you would like to see
# No density perturbations or vel/disp pert for alfven modes.
show_mag = True
#show_mag_scale = True #must also have show_mag = True
show_mag_fade = True
#show_mag_vec = True
show_vel_top = True
#show_disp_top = True
show_axes = True
#show_axis_labels = True
show_mini_axis = True

# Video resolution
#res = (1920,1080)
res = tuple(101 * np.array((16,9)))
#res = tuple(51 * np.array((16,9)))
#res = tuple(21 * np.array((16,9)))

number_of_frames = 1

make_video = False
#make_video = True

mode = 'alfven-mixed-driver'

#
##
###
####
#####
######
#######
########
#########

    
print('Starting ' + mode)
    
# Specify oscillation parameters
K = 2.

W = vA1
    
# Dependent variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t
    
#################################################################################
    

xmin = 0.
xmax = 4.
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
t_end = 2*zmax
        
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




for t_ind in range(nt):
    xvals = np.linspace(xmin, xmax, nx)
    yvals = np.linspace(ymin, ymax, ny)
    zvals = np.linspace(zmin, zmax, nz, endpoint=False)
    
    xixvals_t = np.real(np.repeat(sf.xix_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    xizvals_t = np.real(np.repeat(sf.xiz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    xiyvals_t = np.real(np.repeat(sf.xiy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    
    if show_vel_top == True:
        vxvals_t = np.real(np.repeat(sf.vx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        vzvals_t = np.real(np.repeat(sf.vz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
        vyvals_t = np.real(np.repeat(sf.vy_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
            
    bxvals_t = np.real(np.repeat(sf.bx_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    byvals_t = np.real(np.repeat(sf.by_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2))
    bz_eq3d = np.repeat(sf.bz_eq_amd(xvals, zvals)[:, :, np.newaxis], ny, axis=2)
    bzvals_t = np.real(np.repeat(-sf.bz_amd(xvals, zvals, t, vA_func)[:, :, np.newaxis], ny, axis=2) +
                       bz_eq3d)
    
    
    
#    bxvals_t = bxvals
#    byvals_t = byvals
#    bzvals_t = bzvals
    
#    if show_disp_top == True:
#        xixvals_t = xixvals
#        xiyvals_t = xiyvals
#        xizvals_t = xizvals
#        
#    if show_vel_top == True:
#        vxvals_t = vxvals
#        vyvals_t = vyvals
#        vzvals_t = vzvals

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
    
    if show_vel_top == True: 
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

#    zvals, yvals = np.mgrid[0:nz:(nz)*1j,
#                            0:ny:(ny)*1j]
        
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
    
                
#    xvals, zvals, yvals = np.mgrid[xmin:xmax:(nx)*1j,
#                                   zmax:zmax:(nz)*1j,
#                                   ymin:ymax:(ny)*1j]
    xvals, zvals, yvals = np.mgrid[0:nx:(nx)*1j,
                                   0:nz:(nz)*1j,
                                   0:ny:(ny)*1j]
        
    field = mlab.pipeline.vector_field(bxvals_t, bzvals_t, byvals_t, name="B field", 
                                           figure=fig, scalars=zvals)
    field.spacing = spacing
                
#    if show_mag == True:
#        # Create an array of points for which we want mag field seeds
#        nx_seed = 15 #7
#        ny_seed = 13 #10
#        start_x = xmin #38
#        end_x = xmax+(xmax-xmin)/nx - start_x
#        start_y = ymin+(ymax-ymin)/ny
#        if ny == 20:
#            end_y = ymax - (ymax-ymin)/ny
#        elif ny == 100:
#            end_y = ymax - 2 * (ymax-ymin)/ny
#        else:
#            end_y = ymax - (ymax-ymin)/ny
#        seeds=[]
#        dx_res = (end_x - start_x) / (nx_seed-1)
#        dy_res = (end_y - start_y) / (ny_seed-1)
#        for j in range(0,ny_seed):
#            for i in range(0,nx_seed):
#                x = start_x + (i * dx_res) * x_spacing
#                y = start_y + (j * dy_res) * y_spacing
#                z = 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
#                seeds.append((x,z,y))
                
    if show_mag == True:
        # Create an array of points for which we want mag field seeds
        nx_seed = 9 #7
        ny_seed = 13 #10
        start_x = 2. #38
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
                z = 1. #+ (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                seeds.append((x,z,y))
        
        for i in range(nx_seed):
            del seeds[0]
            del seeds[-1]
        

        seeds = msp.move_seeds(seeds, xixvals_t, xiyvals_t, xizvals_t, 
                               [xmin, ymin, zmin], [xmax, ymax, zmax])
        
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
            mag_fade_value = 20
            mag_lut[:mag_fade_value,-1] = np.linspace(0, 255, mag_fade_value)
            mag_lut[-mag_fade_value:,-1] = np.linspace(255, 0, mag_fade_value)
            module_manager.scalar_lut_manager.lut.table = mag_lut
            
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
            
            
    if show_vel_top == True:
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
    
    if view == 'front-top-side':
        field.scene.camera.position = [400.86836795744739, 181.09643412881843, 495.9729949005914]
        field.scene.camera.focal_point = [50.799999999999997, 50.399999999999999, 50.200000000000003]
        field.scene.camera.view_angle = 14.0
        field.scene.camera.view_up = [-0.14482357103326407, 0.9744012643898321, -0.17195438124301951]
        field.scene.camera.clipping_range = [418.37366265114053, 789.30998655093481]

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
        mpf.uniform_lighting(field)
    
    #Black background
    mpf.background_colour(field, (0., 0., 0.))

        
        
# Trying and failing to sort out memory issues.
#        mlab.gcf()
#        mlab.clf()
#        mlab.close()
#        gc.collect()
#        del fig
#        engine_manager.current_engine = None
#        registry.engines = {}
        
    if make_video == True:
        prefix = 'amd_' + view + '_' + mode
        mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\3D_vis_animations\\'
                 + prefix + str(t_ind+1) + '.png')
        mlab.close(fig)
        
    t = t + (t_end - t_start) / nt
    del vxvals_t
if make_video == True:
#    i2v.image2video(prefix=prefix, output_name=prefix+'_video', out_extension='mp4', fps=20, n_loops=4, delete_images=True, delete_old_videos=True, res=res[1])
    i2v.image2video(prefix='amd_front-top-side_alfven-mixed-driver', output_name='video', 
                    out_extension='mp4', fps=20, n_loops=1, delete_images=True,
                    delete_old_videos=True, cover_page=True)
print('Finished ' + mode)
    