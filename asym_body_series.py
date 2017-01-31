
import sys
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
#sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')
sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\ffmpeg')
#sys.path.append(u'W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/')

import numpy as np
import scipy as sc

from scipy.optimize import newton
import solver

import slab_functions as sf
from pysac.plot.mayavi_seed_streamlines import SeedStreamline

from mayavi import mlab
from mayavi.modules.image_plane_widget import ImagePlaneWidget

import img2vid as i2v

# Which angle shall we view from?
#view = 'front'
#view = 'front parallel'
#view = 'top'
#view = 'top parallel'
view = 'front top'

show_density = False
show_mag = False
show_mag_scale = False
show_vel_front = False
show_vel_top = False
show_axes = False

#show_density = True
show_mag = True
#show_mag_scale = True
show_vel_front = True
show_vel_top = True
show_axes = True


# Specify oscillation parameters
K = 2.
R1 = 1.8

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)

# Find a solution for W
ntries = 50

W_init = 0.2
W_end = sf.cT
W = W_init

W = solver.solver_forwards(disp_rel_asym_1var, W_init, W_end, ntries)


# Dependent variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

xmin = -2.*K
xmax = 2.*K
ymin = 0.
ymax = 4.
zmin = 0.
zmax = 2*np.pi

nx = 99
ny = 99
nz = 99
nt = 20

t_start = 0.
t_end = zmax

t = t_start
for t_ind in range(nt):
    xvals = np.linspace(xmin, xmax, nx)
    yvals = np.linspace(ymin, ymax, ny+1)
    zvals = np.linspace(zmin, zmax, nz+1)
    
    vxvals = np.real(np.repeat(sf.vx_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2))
    vzvals = np.real((1j * sf.c0**2 / (sf.c0**2 - W**2)) * 
                     np.repeat(sf.vx_dash_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2))
    vyvals = np.zeros_like(vxvals)
    
#    vxvals_mask = np.copy(vxvals)
#    vyvals_mask = np.copy(vyvals)
#    vzvals_mask = np.copy(vzvals)
#
#    for i in range(vxvals.shape[0]):
#        for j in range(vxvals.shape[1]):
#            for k in range(vxvals.shape[2]):
#                if (i%mod)!=0 or (j%mod)!=0 or (k%mod)!=0:
#                    vxvals_mask[i,j,k] = 0
#                    vzvals_mask[i,j,k] = 0
    
    
    mod = 7
    vxvals_mask_top = np.copy(vxvals)
    vyvals_mask_top = np.copy(vyvals)
    vzvals_mask_top = np.copy(vzvals)
    
    for i in range(vxvals.shape[0]):
        for j in range(vxvals.shape[1]):
            for k in range(vxvals.shape[2]):
                if (i%mod)!=0 or (k%mod)!=0:
                    vxvals_mask_top[i,j,k] = 0
                    vzvals_mask_top[i,j,k] = 0
    if show_vel_front == True:
        vxvals_mask_front = np.copy(vxvals)
        vyvals_mask_front = np.copy(vyvals)
        vzvals_mask_front = np.copy(vzvals)
        
        for i in range(vxvals.shape[0]):
            for j in range(vxvals.shape[1]):
                for k in range(vxvals.shape[2]):
                    if (i%mod)!=0 or (j%mod)!=0:
                        vxvals_mask_front[i,j,k] = 0
                        vzvals_mask_front[i,j,k] = 0

    bxvals = sf.B0*vxvals / W
    bz_eq2d = np.repeat(sf.bz_eq(xvals, K)[:, np.newaxis], nz+1, axis=1)
    bz_eq3d = np.repeat(bz_eq2d[:, :, np.newaxis], ny+1, axis=2)
    bzvals = np.real(((-1j*sf.B0 / W)*np.repeat(sf.vx_dash_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2) +
            bz_eq3d))
    # Maybe there should be a - in bzvals? Seems to work without though.
    byvals = np.zeros_like(bxvals)

    p_totvals = np.real(np.repeat(sf.p_tot_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2))
    
    xvals, zvals, yvals = np.mgrid[xmin:xmax:(nx)*1j,
                                   zmin:zmax:(nz+1)*1j,
                                   ymin:ymax:(ny+1)*1j]



    fig = mlab.figure(size=tuple(30 * np.array((16,9)))) # size=(1000,1000))
    
    if show_density == True:
        # Scalar field p_tot    
        sca = mlab.pipeline.scalar_field(p_totvals, name="p_tot", figure=fig)
        #scalar_cut_plane = ScalarCutPlane()
        #fig.parent.add_filter(scalar_cut_plane, sca)
        
        minp = p_totvals.min()
        maxp = p_totvals.max()
        
        #Volume for high pressure
        pvmin1 = minp+ 0.55 * (maxp - minp)
        pvmax1 = minp + 1. * (maxp - minp)
        vol1 = mlab.pipeline.volume(sca, vmin=pvmin1, vmax=pvmax1)
        
        # Changing the ctf:
        from tvtk.util.ctf import ColorTransferFunction
        ctf1 = ColorTransferFunction()
        ctf1.add_rgb_point(pvmin1, 1., 1., 0.5)
        ctf1.add_rgb_point(pvmin1 + 0.5 * (pvmax1 - pvmin1), 1, 0.3, 0.1)
        ctf1.add_rgb_point(pvmax1, 1., 0., 0.)
        # ...
        vol1._volume_property.set_color(ctf1)
        vol1._ctf = ctf1
        vol1.update_ctf = True
        
        #Changing the opacity of the volume vol1
        ## Changing the otf:
        from tvtk.util.ctf import PiecewiseFunction
        otf1 = PiecewiseFunction()
        otf1.add_point(pvmin1, 0)
        otf1.add_point(pvmax1, 0.20)
        ##vol1._otf = otf
        vol1._volume_property.set_scalar_opacity(otf1)
        
        
        #Volume for low pressure
        pvmin2 = minp+ 0. * (maxp - minp)
        pvmax2 = minp + 0.45 * (maxp - minp)
        vol2 = mlab.pipeline.volume(sca, vmin=pvmin2, vmax=pvmax2)
        
        # Changing the ctf:
        ctf2 = ColorTransferFunction()
        ctf2.add_rgb_point(pvmin2, 0., 0.5, 1.)
        ctf2.add_rgb_point(pvmin2 + 0.5 * (pvmax2 - pvmin2), 0.1, 0.7, 1.)
        ctf2.add_rgb_point(pvmax2, 0.5, 1., 1.)
        # ...
        vol2._volume_property.set_color(ctf2)
        vol2._ctf = ctf2
        vol2.update_ctf = True
        
        #Changing the opacity of the volume vol2
        ## Changing the otf:
        otf2 = PiecewiseFunction()
        otf2.add_point(pvmax2, 0)
        otf2.add_point(pvmin2, 0.20)
        ##vol1._otf = otf
        vol2._volume_property.set_scalar_opacity(otf2)
    
    
    
                                       
    #image_plane_widget = ImagePlaneWidget()
    #fig.parent.add_filter(image_plane_widget, sca)
    #image_plane_widget.ipw.plane_orientation = 'y_axes'
    #
    #image_plane_widget2 = ImagePlaneWidget()
    #fig.parent.add_filter(image_plane_widget2, sca)
    #image_plane_widget2.ipw.plane_orientation = 'z_axes'
    #
    #module_manager = image_plane_widget.parent
    #module_manager.scalar_lut_manager.lut_mode = 'RdBu'
    #module_manager.scalar_lut_manager.reverse_lut = True
    
    if show_mag == True:
        # Vector field bxvals, bzvals, byvals
        field = mlab.pipeline.vector_field(bxvals, bzvals, byvals, name="B field",
                                           figure=fig)
        
        
        magnitude = mlab.pipeline.extract_vector_norm(field)
        #contours = mlab.pipeline.iso_surface(magnitude,
        #                                        contours=range(2, 14, 3),
        #                                        transparent=True,
        #                                        opacity=0.4,
        #                                        colormap='YlGnBu',
        #                                        vmin=0, vmax=14)
        
        # Create an array of points for which we want mag field seeds
        nx_seed = 5
        ny_seed = 10
        start_x = 38 #38
        end_x = nx - start_x
        start_y = 0
        end_y = ny
        seeds=[]
        dx_res = (end_x - start_x) / nx_seed
        dy_res = (end_y - start_y) / ny_seed
        for i in range(0,nx_seed+2):
            for j in range(0,ny_seed+2):
                x= start_x + (i * dx_res)
                y= start_y + (j * dy_res)
                z= 1. + (t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
                seeds.append((x,z,y))
                                                     
        field_lines = SeedStreamline(seed_points=seeds)
        field_lines.stream_tracer.integration_direction='both'
        field_lines.streamline_type = 'tube'
        
        magnitude = mlab.pipeline.extract_vector_norm(field)
        magnitude.add_child(field_lines)
        module_manager = field_lines.parent
        module_manager.scalar_lut_manager.lut_mode = 'Reds'
        module_manager.scalar_lut_manager.data_range=[-30,25]
        field_lines.stream_tracer.maximum_propagation = 500.
        field_lines.tube_filter.number_of_sides = 5
        field_lines.tube_filter.radius = 0.7
        field_lines.tube_filter.capping = True
        if show_mag_scale == True:
            module_manager.scalar_lut_manager.lut_mode = 'jet'
            module_manager.scalar_lut_manager.data_range=[7,18]
                
    
    
    
    
    #vdirfield = mlab.pipeline.vector_field(vxvals_mask, vzvals_mask, vyvals_mask,
    #                                       name="V field", figure=fig)
    #vectors = mlab.pipeline.vectors(vdirfield, scale_factor=15.)
    ##vectors.glyph.glyph.orient = False
    #vectors.glyph.color_mode = 'no_coloring'
    
    if show_vel_top == True:
        vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top, np.zeros_like(vxvals_mask_top),
                                                    vyvals_mask_top, name="V field top",
                                                    figure=fig)
        vector_cut_plane_top = mlab.pipeline.vector_cut_plane(vdirfield_top, 
                                                              scale_factor=8.)
        vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
        vector_cut_plane_top.glyph.color_mode = 'no_coloring'
        vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,100., 50.5])
        vector_cut_plane_top.implicit_plane.widget.enabled = False
        vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
    
    if show_vel_front == True:
        vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front, vzvals_mask_front,
                                                     vyvals_mask_front, name="V field front",
                                                     figure=fig)
        vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                                scale_factor=8.)
        vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
        vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, 100])
        vector_cut_plane_front.glyph.color_mode = 'no_coloring'
        vector_cut_plane_front.implicit_plane.widget.enabled = False
        vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
        vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'

    if show_axes == True:
        axes = mlab.axes(field, nb_labels=1, line_width=3)
        axes.axes.x_label = ''
        axes.axes.y_label = ''
        axes.axes.z_label = ''
        axes.axes.label_format = ''
    
    #Set viewing angle
    if view == 'front parallel':
        field.scene.parallel_projection = True
        field.scene.z_plus_view()
        field.scene.camera.view_angle = 21
    if view == 'front':
        field.scene.z_plus_view()
        field.scene.camera.view_angle = 21
    if view == 'top':
        field.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
        field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
        field.scene.camera.view_angle = 21.0
        field.scene.camera.view_up = [-0, 0, -1]
        field.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
    if view == 'top parallel':
        field.scene.parallel_projection = True
        field.scene.camera.position = [53.107781380642741, 523.35670183503294, 50.948508989758153]
        field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
        field.scene.camera.view_angle = 21.0
        field.scene.camera.view_up = [-0, 0, -1]
        field.scene.camera.clipping_range = [368.83220888718552, 605.15289607145894]
    if view == 'front top':
    field.scene.camera.position = [48.764852970361503, 223.64895482756552, 498.62216293273576]
    field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
    field.scene.camera.view_angle = 21.0
    field.scene.camera.view_up = [-0.002418791139063777, 0.93281530024654913, -0.36034672896443193]
    field.scene.camera.clipping_range = [345.97885880654962, 650.71850659694883]
        
#    field.scene.camera.position = [249.04443094558394, 46.430453162279896, 318.42236533253191]
#    field.scene.camera.focal_point = [50.5, 50.5, 50.5]
#    field.scene.camera.view_angle = 23
#    field.scene.camera.view_up = [-0.0015851917479672694, 0.99986487801935642, 0.016361933579490465]
#    field.scene.camera.clipping_range = [192.05928737949483, 513.00739425263941]

    #Black background
    field.scene.background = (0., 0., 0.)
    
#    mlab.savefig(u'/media/matthew/W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/'
#                 + '{0:02d}'.format(t_ind+1) + '.png')
#    mlab.savefig(u'/media/matthew/W7_DATA/my_work/projects/Asymmetric_slab/Python/visualisations/ffmpeg/'
#                 + str(t_ind+1) + '.png')
    mlab.savefig('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations\\ffmpeg\\'
                 + str(t_ind+1) + '.png')

    mlab.close()
    t = t + (t_end - t_start)/nt

#i2v.img2vid(prefix='%01d', output_name='video', fps=10, n_loops=4)