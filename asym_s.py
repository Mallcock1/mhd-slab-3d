# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 11:07:22 2016

@author: Matt
"""

import slab_functions as sf
import solver
import numpy as np
from pysac.plot.mayavi_seed_streamlines import SeedStreamline

from mayavi import mlab
from mayavi.tools.sources import vector_field, scalar_field
from mayavi.modules.image_plane_widget import ImagePlaneWidget

# oscillation parameters
K = 2.
R1 = 1.5

# Dependent variables:
# x = k*x
# y = k*y
# z = k*z
# W = omega/k
# K = k*x_0
# t = omega*t

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)
    
W_init = 0.1
W_fin = sf.cT
W = solver.solver_forwards(disp_rel_asym_1var, W_init, W_fin, 10)

xmin = -2.*K
xmax = 2.*K
ymin = 0.
ymax = 4.
zmin = 0.
zmax = 2.*np.pi

nx = 100
ny = 100
nz = 100
nt = 20

#xvals, zvals, tvals = np.meshgrid(np.linspace(xmin, xmax, nx+1),
#                                  np.linspace(zmin, zmax, nz+1),
#                                  np.linspace(tmin, tmax, nt+1))
#                            


t = 0

xvals = np.linspace(xmin, xmax, nx)
yvals = np.linspace(ymin, ymax, ny+1)
zvals = np.linspace(zmin, zmax, nz+1)

vxvals = np.repeat(sf.vx_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2)

bxvals = -sf.B0*vxvals / W
bz_eq2d = np.repeat(sf.bz_eq(xvals, K)[:, np.newaxis], nz+1, axis=1)
bz_eq3d = np.repeat(bz_eq2d[:, :, np.newaxis], ny+1, axis=2)
bzvals = ((-sf.B0 / W)*np.repeat(sf.vx_dash_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2) +
            bz_eq3d)
byvals = np.zeros_like(bxvals)

p_totvals = np.repeat(sf.p_tot_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2)


fig = mlab.figure()

# Scalar field p_tot    
sca = scalar_field(p_totvals, name="p_tot", figure=fig)
#scalar_cut_plane = ScalarCutPlane()
#fig.parent.add_filter(scalar_cut_plane, sca)


min = p_totvals.min()
max = p_totvals.max()
#vol1 = mlab.pipeline.volume(sca, vmin=min + 0.65 * (max - min),
#                                   vmax=min + 0.9 * (max - min))
#vol2 = mlab.pipeline.volume(sca, vmin=min + 0.1 * (max - min),
#                                   vmax=min + 0.35 * (max - min))
#vol.parent.scalar_lut_manager.lut_mode = 'RdBu'
#vol = Volume()
#fig.parent.add_filter(vol, sca)

                                   
image_plane_widget = ImagePlaneWidget()
fig.parent.add_filter(image_plane_widget, sca)
image_plane_widget.ipw.plane_orientation = 'y_axes'

image_plane_widget2 = ImagePlaneWidget()
fig.parent.add_filter(image_plane_widget2, sca)
image_plane_widget2.ipw.plane_orientation = 'z_axes'

module_manager = image_plane_widget.parent
module_manager.scalar_lut_manager.lut_mode = 'RdBu'
module_manager.scalar_lut_manager.reverse_lut = True
#image_plane_widget = ImagePlaneWidget()
#fig.parent.add_filter(image_plane_widget, sca)
#image_plane_widget.ipw.plane_orientation = 'y_axes'
#module_manager = image_plane_widget.parent
#module_manager.scalar_lut_manager.lut_mode = 'RdBu'


# Vector field bxvals, bzvals, byvals
field = mlab.pipeline.vector_field(bxvals, bzvals, byvals, name="B field", figure=fig)
#field.scene.camera.position = [936.7613264352542, 841.53008573935188, 2193.2985708634601]
#field.scene.camera.focal_point = [251.0, 388.1923828125, 250.5]
#field.scene.camera.view_angle = 30.0
#field.scene.camera.view_up = [-0.058034716603158412, 0.97653691165056677, -0.20738281474790629]
#field.scene.camera.clipping_range = [1380.0778003475418, 3098.5523108029097]


field.scene.camera.position = [140.24563591382554, 134.74789058452146, 360.44654358910702]
field.scene.camera.focal_point = [50.5, 50.5, 50.5]
field.scene.camera.view_angle = 30.0
field.scene.camera.view_up = [-0.20792470843751543, 0.95746946474763595, -0.20004884327845979]
field.scene.camera.clipping_range = [186.02517744318791, 519.93329980779686]

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
start_x = 38
end_x = nx - start_x
start_y = 0
end_y = ny +1
seeds=[]
dx_res = (end_x - start_x) / nx_seed
dy_res = (end_y - start_y) / ny_seed
for i in range(0,nx_seed+2):
    for j in range(0,ny_seed+1):
        x= start_x + (i * dx_res)
        y= start_y + (j * dy_res)
        z= 1.
        seeds.append((x,z,y))
                                             
field_lines = SeedStreamline(seed_points=seeds)
field_lines.stream_tracer.integration_direction='both'
field_lines.streamline_type = 'tube'

magnitude = mlab.pipeline.extract_vector_norm(field)
magnitude.add_child(field_lines)
module_manager = field_lines.parent
module_manager.scalar_lut_manager.lut_mode = 'Reds'
module_manager.scalar_lut_manager.data_range=[0,14]
#field_lines.streamline_type = 'tube'
field_lines.stream_tracer.maximum_propagation = 500.
field_lines.tube_filter.number_of_sides = 5
field_lines.tube_filter.radius = 0.7

axes = mlab.axes(field)
axes.axes.x_label = ''
axes.axes.y_label = ''
axes.axes.z_label = ''
mlab.show()