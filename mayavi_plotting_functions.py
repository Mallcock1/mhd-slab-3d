# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 16:16:47 2017

@author: Matt
"""

from mayavi import mlab

def uniform_lighting(scene):
    scene.scene.light_manager.number_of_lights = 6
    
    camera_light1 = scene.scene.light_manager.lights[0]
    camera_light1.activate = True
    camera_light1.intensity = 0.7
    camera_light1.elevation = 90.
    camera_light1.azimuth = 0.

    camera_light2 = scene.scene.light_manager.lights[1]
    camera_light2.activate = True
    camera_light2.intensity = 0.7
    camera_light2.elevation = -90.
    camera_light2.azimuth = 0.

    camera_light3 = scene.scene.light_manager.lights[2]
    camera_light3.activate = True
    camera_light3.intensity = 0.7
    camera_light3.elevation = 0.
    camera_light3.azimuth = -90

    camera_light4 = scene.scene.light_manager.lights[3]
    camera_light4.activate = True
    camera_light4.intensity = 0.7
    camera_light4.elevation = 0.
    camera_light4.azimuth = 0.

    camera_light5 = scene.scene.light_manager.lights[4]
    camera_light5.activate = True
    camera_light5.intensity = 0.7
    camera_light5.elevation = 0.
    camera_light5.azimuth = 90.

    camera_light6 = scene.scene.light_manager.lights[5]
    camera_light6.activate = True
    camera_light6.intensity = 0.7
    camera_light6.elevation = 0.
    camera_light6.azimuth = 180.

def background_colour(scene, colour):
    """
    colour = tuple colour code
    e.g (0.,0.,0.) for black
    """
    scene.scene.background = colour

def mini_axes():
    oa = mlab.orientation_axes(xlabel='x', ylabel='z', zlabel='y')
    oa.marker.set_viewport(0,0,0.25,0.25) # minx, miny, maxx, maxy

def axes(scene, show_axis_labels, view):
    axes = mlab.axes(scene, nb_labels=1, line_width=3)
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
        
def view_position(scene, view, nx, ny, nz):
    #Set viewing angle
    if view == 'front-parallel':
        scene.scene.z_plus_view()
        scene.scene.parallel_projection = True
        scene.scene.camera.zoom(1.65) # Parallel projection zoom is done in this way, different to perspective projection
    if view == 'front':
        scene.scene.z_plus_view()
        scene.scene.camera.view_angle = 21.
    if view == 'top':
        scene.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
        scene.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
        scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [-0, 0, -1]
        scene.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
    if view == 'top-parallel':
        scene.scene.parallel_projection = True
        scene.scene.camera.zoom(2.)
        scene.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
        scene.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
#            scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [-0, 0, -1]
        scene.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
    if view == 'front-top':
        scene.scene.camera.position = [48.764852970361503, 223.64895482756552, 498.62216293273576]
        scene.scene.camera.focal_point = [50.821544647216797, 46., 50.159849926829338]
        scene.scene.camera.view_angle = 16.0
        scene.scene.camera.view_up = [-0.002418791139063777, 0.93281530024654913, -0.36034672896443193]
        scene.scene.camera.clipping_range = [345.97885880654962, 650.71850659694883]
     
    if view == 'front-side':
        scene.scene.camera.position = [126.6 * nx / 100., 60.5 * nz / 100., 524.8 * ny / 100.]
        scene.scene.camera.focal_point = [50.8 * nx / 100., 50.4 * nz / 100., 50.2 * ny / 100.]
        scene.scene.camera.view_angle = 14.
        scene.scene.camera.view_up = [-0.01695 * nx / 100., 0.999686 * nz / 100., -0.0184509 * ny / 100.]
        scene.scene.camera.clipping_range = [366.21083458278804, 631.07664372567524]
        

def mask_points(var_x, var_y, var_z, mod=4, front_or_top):
    if front_or_top == 'front':
        var_x_mask = np.copy(var_x)
        var_y_mask = np.copy(var_y)
        var_z_mask = np.copy(var_z)
        
        for i in range(var_x_mask.shape[0]):
            for j in range(var_x_mask.shape[1]):
                for k in range(var_x_mask.shape[2]):
                    if (i%mod) != 1 or (j%mod) != 1:
                        var_x_mask[i,j,k] = 0.
                        var_z_mask[i,j,k] = 0.
                        
    if front_or_top == 'top':
        var_x_mask = np.copy(var_x)
        var_y_mask = np.copy(var_y)
        var_z_mask = np.copy(var_z)
        
        for i in range(var_x_mask.shape[0]):
            for j in range(var_x_mask.shape[1]):
                for k in range(var_x_mask.shape[2]):
                    if (i%mod) != 1 or (k%mod_top) != 1:
                        var_x_mask[i,j,k] = 0.
                        var_z_mask[i,j,k] = 0.
                        
    return [var_x_mask, var_y_mask, var_z_mask]

def vector_cut_plane(vec_field, scale_factor=4, spacing):
        
    scale_factor = 4. # scale factor for direction field vectors
    
    if front_or_top == 'front':
        vecfield_front = mlab.pipeline.vector_field(var_x, var_y, var_z, figure=fig)
        vecfield_front.spacing = spacing
        vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                          scale_factor=scalefactor)
        vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
        vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
        vector_cut_plane_front.glyph.color_mode = 'no_coloring'
        vector_cut_plane_front.implicit_plane.widget.enabled = False
        vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
    
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
    
    
    if show_vel_top == True or show_vel_top_pert == True:
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
        
    if show_vel_front == True or show_vel_front_pert == True:
        vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front_t, vzvals_mask_front_t,
                                                     vyvals_mask_front_t, name="V field front",
                                                     figure=fig)
        vdirfield_front.spacing = spacing
        vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                          scale_factor=scalefactor)
        vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
        vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
        vector_cut_plane_front.glyph.color_mode = 'no_coloring'
        vector_cut_plane_front.implicit_plane.widget.enabled = False
        vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
    
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
        
    if show_disp_front == True:
        xidirfield_front = mlab.pipeline.vector_field(xixvals_mask_front_t, xizvals_mask_front_t,
                                                     xiyvals_mask_front_t, name="Xi field front",
                                                     figure=fig)
        xidirfield_front.spacing = spacing
        vector_cut_plane_front = mlab.pipeline.vector_cut_plane(xidirfield_front, 
                                                          scale_factor=scalefactor)
        vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
        vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
        vector_cut_plane_front.glyph.color_mode = 'no_coloring'
        vector_cut_plane_front.implicit_plane.widget.enabled = False
        vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'
    