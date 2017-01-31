# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 17:12:38 2016

@author: Matt
"""

import numpy as np
from pysac.plot.mayavi_seed_streamlines import SeedStreamline

from mayavi import mlab
from mayavi.tools.sources import vector_field, scalar_field
from mayavi.modules.scalar_cut_plane import ScalarCutPlane
from mayavi.modules.image_plane_widget import ImagePlaneWidget

# Define the sound speeds and alfven speeds.
ce = 1.
c0 = 1.
vA = 2.
cT = np.sqrt(c0**2 * vA**2*(c0**2 + vA**2)**(-1))

B0=10

constC = 0.1

# Specify oscillation parameters
W = 0.7
K = 2.

# Dependent variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

m0 = np.sqrt((c0**2 - W**2)*(vA**2 - W**2) / ((c0**2 + vA**2)*(cT**2 - W**2)))

me = np.sqrt(1 - W**2*ce**(-2))

# Define density ratio using dispersion relation
Re = (vA**2-W**2)*me / W**2*m0*(np.tanh(m0*K))  # rho_e/rho_0

lamb0 = -(vA**2-W**2)/(m0*W)
lamb00 = W/m0

lambe = Re*W/me

constB = 0. #constC*((lamb0*np.cosh(m0*K)+lambe*np.sinh(m0*K)) /
#                 (lambe*np.cosh(m0*K)+lamb0*np.sinh(m0*K)))

constA = ((constB*np.cosh(m0*K) - constC*np.sinh(m0*K)) /
            (np.cosh(me*K) - np.sinh(me*K)))

constD = ((constB*np.cosh(m0*K) + constC*np.sinh(m0*K)) /
            (np.cosh(me*K) - np.sinh(me*K)))

def vxhat(x):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhatfunction = np.zeros(len(x))
    for i in indices:
        vxhatfunction[i] = constB*np.cosh(m0*x[i]) + constC*np.sinh(m0*x[i])
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhatfunction[i] = constA*(np.cosh(me*x[i]) + np.sinh(me*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhatfunction[i] = constD*(np.cosh(me*x[i]) - np.sinh(me*x[i]))
    return vxhatfunction
    
def vxhat_dash(x):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    vxhat_dashfunction = np.zeros(len(x))
    for i in indices:
        vxhat_dashfunction[i] = m0*constC*np.cosh(m0*x[i])
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        vxhat_dashfunction[i] = me*constA*(np.sinh(me*x[i]) + np.cosh(me*x[i]))
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        vxhat_dashfunction[i] = me*constD*(np.sinh(me*x[i]) - np.cosh(me*x[i]))
    return vxhat_dashfunction
    
def bz_eq(x):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    bz_eqfunction = np.zeros(len(x))
    for i in indices:
        bz_eqfunction[i] = B0
    return bz_eqfunction
    
def p_tothat(x):
    truth = np.array(np.abs(x) <= K*np.ones(len(x)))
    indices = np.where(truth == True)
    p_tothatfunction = np.zeros(len(x))
    for i in indices:
        p_tothatfunction[i] = m0*constC*np.cosh(m0*x[i]) * lamb0
    truth2 = np.array(x < -K*np.ones(len(x)))
    indices2 = np.where(truth2 == True)
    for i in indices2:
        p_tothatfunction[i] = me*constA*(np.sinh(me*x[i]) + np.cosh(me*x[i])) * lambe
    truth3 = np.array(x > K*np.ones(len(x)))
    indices3 = np.where(truth3 == True)
    for i in indices3:
        p_tothatfunction[i] = me*constD*(np.sinh(me*x[i]) - np.cosh(me*x[i])) * lambe
    return p_tothatfunction

def vx(x, z, t):
    return np.outer(vxhat(x), np.cos(z-t))

def vx_dash(x, z, t):
    return np.outer(vxhat_dash(x), np.sin(z-t))
    
def p_tot(x, z, t):
    return np.outer(p_tothat(x), np.sin(z-t))

xmin = -4.
xmax = 4.
ymin = 0.
ymax = 4.
zmin = 0.
zmax = 8.

nx = 99
ny = 99
nz = 99
nt = 10

#xvals, zvals, tvals = np.meshgrid(np.linspace(xmin, xmax, nx+1),
#                                  np.linspace(zmin, zmax, nz+1),
#                                  np.linspace(tmin, tmax, nt+1))
#                            

t_start = 0.
t_end = zmax*0.7

t = t_start
for t_ind in range(nt+1):
    xvals = np.linspace(xmin, xmax, nx)
    yvals = np.linspace(ymin, ymax, ny+1)
    zvals = np.linspace(zmin, zmax, nz+1)
    
    vxvals = np.repeat(vx(xvals, zvals, t)[:, :, np.newaxis], ny+1, axis=2)
    
    bxvals = -B0*vxvals / W
    bz_eq2d = np.repeat(bz_eq(xvals)[:, np.newaxis], nz+1, axis=1)
    bz_eq3d = np.repeat(bz_eq2d[:, :, np.newaxis], ny+1, axis=2)
    bzvals = ((-B0 / W)*np.repeat(vx_dash(xvals, zvals, t)[:, :, np.newaxis], ny+1, axis=2) +
                bz_eq3d)
    byvals = np.zeros_like(bxvals)
    
    p_totvals = np.repeat(p_tot(xvals, zvals, t)[:, :, np.newaxis], ny+1, axis=2)
    
    xvals, zvals, yvals = np.mgrid[xmin:xmax:(nx)*1j,
                                   zmin:zmax:(nz+1)*1j,
                                   ymin:ymax:(ny+1)*1j]
    
    #
    ##
    
    
    fig = mlab.figure(size=(1000,1000))
    
    # Scalar field p_tot    
    sca = scalar_field(p_totvals, name="p_tot", figure=fig)
    
#    min = p_totvals.min()
#    max = p_totvals.max()
#    vol = mlab.pipeline.volume(sca, vmin=min + 0.65 * (max - min),
#                                   vmax=min + 0.9 * (max - min))
    
    image_plane_widget = ImagePlaneWidget()
    fig.parent.add_filter(image_plane_widget, sca)
    image_plane_widget.ipw.plane_orientation = 'y_axes'
    
    image_plane_widget2 = ImagePlaneWidget()
    fig.parent.add_filter(image_plane_widget2, sca)
    image_plane_widget2.ipw.plane_orientation = 'z_axes'
    
    module_manager = image_plane_widget.parent
    module_manager.scalar_lut_manager.lut_mode = 'RdBu'
    module_manager.scalar_lut_manager.reverse_lut = True
    
    
    # Vector field bxvals, bzvals, byvals
    field = mlab.pipeline.vector_field(bxvals, bzvals, byvals, name="B field", figure=fig)
    
    field.scene.camera.position = [140.24563591382554, 134.74789058452146, 360.44654358910702]
    field.scene.camera.focal_point = [50.5, 50.5, 50.5]
    field.scene.camera.view_angle = 30.0
    field.scene.camera.view_up = [-0.20792470843751543, 0.95746946474763595, -0.20004884327845979]
    field.scene.camera.clipping_range = [186.02517744318791, 519.93329980779686]
    
#    field.scene.camera.position = [50.0, 50.5, 382.87955413300307]
#    field.scene.camera.focal_point = [50.0, 50.5, 50.5]
#    field.scene.camera.view_angle = 30.0
#    field.scene.camera.view_up = [0.0, 1.0, 0.0]
#    field.scene.camera.clipping_range = [229.55575859167305, 462.61524744499809]    
        
#    field.scene.z_plus_view()
    
    field.scene.background = (0.12156862745098039, 0.12156862745098039, 0.12156862745098039)
    
    
    
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
    start_x = 35
    end_x = nx - start_x
    start_y = 0
    end_y = ny+1
    seeds=[]
    dx_res = (end_x - start_x) / nx_seed
    dy_res = (end_y - start_y) / ny_seed
    for i in range(0,nx_seed+2):
        for j in range(0,ny_seed+1):
            x= start_x + (i * dx_res)
            y= start_y + (j * dy_res)
            z= (1+t_start + t_ind*(t_end - t_start)/nt)/zmax * nz
            seeds.append((x,z,y))
                                                 
    field_lines = SeedStreamline(seed_points=seeds)
    field_lines.stream_tracer.integration_direction='both'
    
    magnitude2 = mlab.pipeline.extract_vector_norm(field)
    magnitude2.add_child(field_lines)
    module_manager = field_lines.parent
    module_manager.scalar_lut_manager.lut_mode = 'Reds'
    module_manager.scalar_lut_manager.data_range=[0,14]
    
    #field_lines.streamline_type = 'tube'
    field_lines.stream_tracer.maximum_propagation = 500.0
    field_lines.tube_filter.number_of_sides = 5
    field_lines.tube_filter.radius = 0.7
    
    ax = mlab.axes(field)
    ax.axes.x_label = ''
    ax.axes.y_label = ''
    ax.axes.z_label = ''

    mlab.savefig(u'D:/my_work/visualisations/sym_field_pressure/new/xz' + str(t_ind+1) + '.png')
    t = t + (t_end - t_start)/nt
    