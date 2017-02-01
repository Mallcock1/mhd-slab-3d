
import numpy as np
import scipy as sc

import sys
sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\visualisations')
sys.path.append('D:\\my_work\\projects\\Asymmetric_slab\\Python\\Mihai')

import Toolbox as tool

from scipy.optimize import newton
import solver

import slab_functions_perturbed_boundary as sf
from pysac.plot.mayavi_seed_streamlines import SeedStreamline

from mayavi import mlab
from mayavi.modules.image_plane_widget import ImagePlaneWidget

###############################################################################

# What mode do you want? OPTIONS:
mode_options = ['slow kink surf', 'slow saus surf', 'slow kink body 1',
                'slow saus body 1', 'slow kink body 2', 'slow saus body 2']
mode = mode_options[2]


# Which angle shall we view from?
#view = 'front'
#view = 'front parallel'
#view = 'top'
#view = 'top parallel'
#view = 'front top'
view = 'front side'

# Uniform lighting?
#uniform_light = True
uniform_light = False

show_density = False
show_density_lagrang = False
show_density_pert = False
show_density_pert2 = False
show_pressure = False
show_mag = False
show_mag_scale = False
show_mag_vec = False
show_vel_front = False
show_vel_front_pert = False
show_vel_top = False
show_vel_top_pert = False
show_disp_top = False
show_disp_front = False
show_axes = False
show_boundary = False

show_density = True
#show_density_lagrang = True
#show_density_pert = True
#show_density_pert2 = True #THIS ONE WORKS
#show_pressure = True
show_mag = True
#show_mag_scale = True
#show_mag_vec = True
show_vel_front = True
#show_vel_front_pert = True
#show_vel_top = True
#show_vel_top_pert = True
#show_disp_top = True
#show_disp_front = True
show_axes = True
show_boundary = True

# Specify oscillation parameters
K = 2.
R1 = 1.8 #1.8 #2.

def disp_rel_asym_1var(W):
    return sf.disp_rel_asym(W, K, R1)
    
def disp_rel_asym_2var(W, K):
    return sf.disp_rel_asym(W, K, R1)


Wrange = np.linspace(0., sf.c2, 51)
Woptions = tool.point_finder_scipy(disp_rel_asym_2var, np.array(K), Wrange, args=None).transpose()
Woptions = np.real(Woptions)
tol = 1e-2
indices_to_rm = []
for i in range(len(Woptions)):
    w = Woptions[i]
    print(type(Woptions[i]))
    print(i)
    if min(abs(np.array([w, w - sf.c0, w - sf.c1(R1), w - sf.c2, w - sf.vA]))) < tol or w < 0 or w > max(sf.c0, sf.vA, sf.c1, sf.c2):
        indices_to_rm.append(i)
Woptions = np.delete(Woptions, indices_to_rm)
Woptions = Woptions.sort()

#Something is wrong below. Fix on 2/2/17 pls
for i in range(len(mode_options)):
    if mode == mode_options[i]:
        W = Woptions[i]
        break

# Find a solution for W
ntries = 50
#

# Slow kink Surface
#W_init = 0.2
#W_end = sf.cT
#W = W_init
#
#W = solver.solver_forwards(disp_rel_asym_1var, W_init, W_end, ntries)

## Slow kink body 1
#W_init = 0.760431
#W_end = 0.760531
#W = W_init
#
#W = solver.solver_forwards(disp_rel_asym_1var, W_init, W_end, ntries)

## Slow kink body 2
W_init = 0.68336
W_end = 0.68338
W = W_init

W = solver.solver_forwards(disp_rel_asym_1var, W_init, W_end, ntries)


#
#W, nothing = tool.point_finder_scipy(disp_rel_asym_2var, (2., 2.), (0.75, 0.8),
#                                    args=(None))
#W = W[0]

#W, nothing = tool.point_finder_scipy(disp_rel_asym_2var, (2., 2.), (0.683, 0.6834),
#                                    args=(None))
#W = W[0]

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

# You can change ny but be careful changing nx, nz.
nx = 100
ny = 20
nz = 100
nt = 20                          

t = 0.

xvals = np.linspace(xmin, xmax, nx)
yvals = np.linspace(ymin, ymax, ny)
zvals = np.linspace(zmin, zmax, nz)

x_spacing = max(nx, ny, nz) / nx
y_spacing = max(nx, ny, nz) / ny
z_spacing = max(nx, ny, nz) / nz

if show_disp_top==True or show_disp_front==True:
    xixvals = np.real(np.repeat(sf.xix_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    xizvals = np.real(np.repeat(sf.xiz_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    xiyvals = np.zeros_like(xixvals)
    
    if show_disp_top==True:    
        mod = 4
        mod_top = np.ceil(4 / y_spacing)
                                
        xixvals_mask_top = np.copy(xixvals)
        xiyvals_mask_top = np.copy(xiyvals)
        xizvals_mask_top = np.copy(xizvals)
        
        for i in range(xixvals.shape[0]):
            for j in range(xixvals.shape[1]):
                for k in range(xixvals.shape[2]):
                    if (i%mod)!=0 or (k%1)!=0:
                        xixvals_mask_top[i,j,k] = 0
                        xizvals_mask_top[i,j,k] = 0
    if show_disp_front==True:
        xixvals_mask_front = np.copy(xixvals)
        xiyvals_mask_front = np.copy(xiyvals)
        xizvals_mask_front = np.copy(xizvals)
        
        for i in range(xixvals.shape[0]):
            for j in range(xixvals.shape[1]):
                for k in range(xixvals.shape[2]):
                    if (i%mod)!=0 or (j%mod)!=0:
                        xixvals_mask_front[i,j,k] = 0
                        xizvals_mask_front[i,j,k] = 0
                        
if show_vel_front == True or show_vel_top == True:
    vxvals = np.real(np.repeat(sf.vx_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    #vzvals = np.real((1j * sf.c0**2 / (sf.c0**2 - W**2)) * 
    #                 np.repeat(sf.vx_dash_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    vzvals = np.real(np.repeat(sf.vz_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    
    vyvals = np.zeros_like(vxvals)
    
    vxvals_mask = np.copy(vxvals)
    vyvals_mask = np.copy(vyvals)
    vzvals_mask = np.copy(vzvals)
    
if show_vel_front_pert == True or show_vel_top_pert == True:
    vxvals = np.real(np.repeat(sf.vx_kink_pert(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))

    vzvals = np.real(np.repeat(sf.vz_kink_pert(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    
    vyvals = np.zeros_like(vxvals)
    
    vxvals_mask = np.copy(vxvals)
    vyvals_mask = np.copy(vyvals)
    vzvals_mask = np.copy(vzvals)

if np.array([show_vel_front, show_vel_top, show_vel_front_pert, show_vel_top_pert]).any() == True:
    mod = 4
    mod_top = np.ceil(4. / y_spacing)
    
    for i in range(vxvals.shape[0]):
        for j in range(vxvals.shape[1]):
            for k in range(vxvals.shape[2]):
                if (i%mod)!=0 or (j%mod)!=0 or (k%mod)!=0:
                    vxvals_mask[i,j,k] = 0
                    vzvals_mask[i,j,k] = 0
                    
    vxvals_mask_top = np.copy(vxvals)
    vyvals_mask_top = np.copy(vyvals)
    vzvals_mask_top = np.copy(vzvals)
    
    for i in range(vxvals.shape[0]):
        for j in range(vxvals.shape[1]):
            for k in range(vxvals.shape[2]):
                if (i%mod)!=0 or (k%mod_top)!=0:
                    vxvals_mask_top[i,j,k] = 0
                    vzvals_mask_top[i,j,k] = 0
                    
    vxvals_mask_front = np.copy(vxvals)
    vyvals_mask_front = np.copy(vyvals)
    vzvals_mask_front = np.copy(vzvals)
    
    for i in range(vxvals.shape[0]):
        for j in range(vxvals.shape[1]):
            for k in range(vxvals.shape[2]):
                if (i%mod)!=0 or (j%mod)!=0:
                    vxvals_mask_front[i,j,k] = 0
                    vzvals_mask_front[i,j,k] = 0

#bxvals = sf.B0*vxvals / W
#bz_eq2d = np.repeat(sf.bz_eq(xvals, K)[:, np.newaxis], nz+1, axis=1)
#bz_eq3d = np.repeat(bz_eq2d[:, :, np.newaxis], ny+1, axis=2)
#bzvals = np.real(((-1j*sf.B0 / W)*np.repeat(sf.vx_dash_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny+1, axis=2) +
#            bz_eq3d))
# Maybe there should be a - in bzvals? Seems to work without though.

bxvals = np.real(np.repeat(sf.bx_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
#bz_eq2d = np.repeat(sf.bz_eq(xvals, K)[:, np.newaxis], nz, axis=1)
bz_eq3d = np.repeat(sf.bz_eq(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2)
#bz_eq3d = sf.B0*np.ones_like(bxvals)
bzvals = np.real(np.repeat(-sf.bz_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2) +
                 bz_eq3d)
                 
                 
byvals = np.zeros_like(bxvals)

bxvals_mask_front = np.copy(bxvals)
byvals_mask_front = np.copy(byvals)
bzvals_mask_front = np.copy(bzvals)

for i in range(bxvals.shape[0]):
    for j in range(bxvals.shape[1]):
        for k in range(bxvals.shape[2]):
            if (i%mod)!=0 or (j%mod)!=0:
                bxvals_mask_front[i,j,k] = 0
                bzvals_mask_front[i,j,k] = 0

if show_pressure == True:
    p_totvals = np.real(np.repeat(sf.p_tot_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))

if show_boundary == True:
    xi_boundary_r_vals = np.real(np.repeat(K + sf.xi_boundary_kink(zvals, t, W, K, R1, boundary='r')[:, np.newaxis], ny, axis=1))
    xi_boundary_l_vals = np.real(np.repeat(-K + sf.xi_boundary_kink(zvals, t, W, K, R1, boundary='l')[:, np.newaxis], ny, axis=1))

#xi_boundary_r_vals = np.real(np.repeat(K*np.ones_like(sf.xi_boundary_kink(zvals, t, W, K, R1, boundary='r'))[:, np.newaxis], ny, axis=1))
#xi_boundary_l_vals = np.real(np.repeat(-K*np.ones_like(sf.xi_boundary_kink(zvals, t, W, K, R1, boundary='r'))[:, np.newaxis], ny, axis=1))
if show_density == True:
    rho_vals = np.real(np.repeat(sf.rho_kink(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))

if show_density_lagrang == True:
    rho_vals_lagrang = np.real(np.repeat(sf.rho_kink_lagrang(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    
if show_density_pert == True:
    rho_vals_pert = np.real(np.repeat(sf.rho_kink_pert(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))
    
if show_density_pert2 == True:
    rho_vals_pert2 = np.real(np.repeat(sf.rho_kink_pert2(xvals, zvals, t, W, K, R1)[:, :, np.newaxis], ny, axis=2))

#zvals, yvals = np.mgrid[zmin:zmax:(nx)*1j,
#                        ymin:ymax:(ny)*1j]

zvals, yvals = np.mgrid[0:nz:(nz)*1j,
                        0:ny:(ny)*1j]

#
##
###
####
#####
####
###
##
#

fig = mlab.figure(size=tuple(100 * np.array((16,9)))) #16:9 aspect ratio for video upload

spacing =  np.array([x_spacing, z_spacing, y_spacing])

if show_boundary == True:
    ext_min_r = ((nx+1) * (xi_boundary_r_vals.min() - xmin) / (xmax - xmin)) * x_spacing
    ext_max_r = ((nx+1) * (xi_boundary_r_vals.max() - xmin) / (xmax - xmin)) * x_spacing
    
    ext_min_l = ((nx+1) * (xi_boundary_l_vals.min() - xmin) / (xmax - xmin)) * x_spacing #plus 2 after (xmax-xmin)?
    ext_max_l = ((nx+1) * (xi_boundary_l_vals.max() - xmin) / (xmax - xmin)) * x_spacing #plus 2 after (xmax-xmin)?
    
    
    boundary_r = mlab.mesh(xi_boundary_r_vals, zvals, yvals,
                           extent=[ext_min_r, ext_max_r, 0, nz, 0, (ny-1) * y_spacing],
                           color=(0.5,0.5,0.5), opacity=0.7)
    
    boundary_l = mlab.mesh(xi_boundary_l_vals, zvals, yvals,
                           extent=[ext_min_l, ext_max_l, 0, nz, 0, (ny-1) * y_spacing],
                           color=(0.5,0.5,0.5), opacity=0.7)

if show_density == True:
    # Scalar field density   
    rho = mlab.pipeline.scalar_field(rho_vals, name="density", figure=fig)
    #scalar_cut_plane = ScalarCutPlane()
    #fig.parent.add_filter(scalar_cut_plane, sca)
    rho.spacing = spacing
    minr = rho_vals.min()
    maxr = rho_vals.max()
    
    #Volume for high pressure
    rvmin1 = minr + 0.55 * (maxr - minr)
    rvmax1 = minr + 1. * (maxr - minr)
    rvol1 = mlab.pipeline.volume(rho, vmin=rvmin1, vmax=rvmax1)
    
    # Changing the ctf:
    from tvtk.util.ctf import ColorTransferFunction
    ctf1 = ColorTransferFunction()
    ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    ctf1.add_rgb_point(rvmin1 + 0.5 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    # ...
    rvol1._volume_property.set_color(ctf1)
    rvol1._ctf = ctf1
    rvol1.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(rvmin1, 0)
    otf.add_point(rvmax1, 0.10)
    ##vol1._otf = otf
    rvol1._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol1.volume_property.shade = False
    rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    
    
    #Volume for low pressure
    rvmin2 = minr + 0. * (maxr - minr)
    rvmax2 = minr + 0.45 * (maxr - minr)
    rvol2 = mlab.pipeline.volume(rho, vmin=rvmin2, vmax=rvmax2)
    
    # Changing the ctf:
    ctf2 = ColorTransferFunction()
    ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    ctf2.add_rgb_point(rvmin2 + 0.5 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    # ...
    rvol2._volume_property.set_color(ctf2)
    rvol2._ctf = ctf2
    rvol2.update_ctf = True

    #Changing the opacity of the volume vol1
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(rvmax2, 0)
    otf.add_point(rvmin2, 0.10)
    ##vol1._otf = otf
    rvol2._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol2.volume_property.shade = False
    rvol2.volume_property.scalar_opacity_unit_distance = 2.0

if show_density_lagrang == True:
    # Scalar field density   
    rho_lagrang = mlab.pipeline.scalar_field(rho_vals_lagrang, name="density", figure=fig)
    #scalar_cut_plane = ScalarCutPlane()
    #fig.parent.add_filter(scalar_cut_plane, sca)
    rho_lagrang.spacing = spacing
    minr = rho_vals_lagrang.min()
    maxr = rho_vals_lagrang.max()
    
    #Volume for high pressure
    rvmin1 = 0.01 #minr + 0.501 * (maxr - minr)
    rvmax1 = 0.5 #minr + 0.51 * (maxr - minr)
    rvol1 = mlab.pipeline.volume(rho_lagrang, vmin=rvmin1, vmax=rvmax1)
    
    # Changing the ctf:
    from tvtk.util.ctf import ColorTransferFunction
    ctf1 = ColorTransferFunction()
    ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    ctf1.add_rgb_point(rvmin1 + 0.5 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    # ...
    rvol1._volume_property.set_color(ctf1)
    rvol1._ctf = ctf1
    rvol1.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(rvmin1, 0)
    otf.add_point(rvmax1, 0.10)
    ##vol1._otf = otf
    rvol1._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol1.volume_property.shade = False
    rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    
    
    #Volume for low pressure
    rvmin2 = -0.5 #minr + 0.49 * (maxr - minr)
    rvmax2 = -0.01 #minr + 0.499 * (maxr - minr)
    rvol2 = mlab.pipeline.volume(rho_lagrang, vmin=rvmin2, vmax=rvmax2)
    
    # Changing the ctf:
    ctf2 = ColorTransferFunction()
    ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    ctf2.add_rgb_point(rvmin2 + 0.5 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    # ...
    rvol2._volume_property.set_color(ctf2)
    rvol2._ctf = ctf2
    rvol2.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(rvmax2, 0)
    otf.add_point(rvmin2, 0.10)
    ##vol1._otf = otf
    rvol2._volume_property.set_scalar_opacity(otf)

    # exempt volume from shading and improve overall look by increasing opacity
    rvol2.volume_property.shade = False
    rvol2.volume_property.scalar_opacity_unit_distance = 2.0
    
if show_density_pert == True:
    # Scalar field density   
    rho_pert = mlab.pipeline.scalar_field(rho_vals_pert, name="density", figure=fig)
    #scalar_cut_plane = ScalarCutPlane()
    #fig.parent.add_filter(scalar_cut_plane, sca)
    rho_pert.spacing = spacing
    minr = rho_vals_pert.min()
    maxr = rho_vals_pert.max()
    
    #Volume for high pressure
    rvmin1 = minr + 0.55 * (maxr - minr)
    rvmax1 = minr + 1. * (maxr - minr)
    rvol1 = mlab.pipeline.volume(rho_pert, vmin=rvmin1, vmax=rvmax1)
    
    # Changing the ctf:
    from tvtk.util.ctf import ColorTransferFunction
    ctf1 = ColorTransferFunction()
    ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    ctf1.add_rgb_point(rvmin1 + 0.5 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    # ...
    rvol1._volume_property.set_color(ctf1)
    rvol1._ctf = ctf1
    rvol1.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(rvmin1, 0)
    otf.add_point(rvmax1, 0.10)
    ##vol1._otf = otf
    rvol1._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol1.volume_property.shade = False
    rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    
    
    #Volume for low pressure
    rvmin2 = minr + 0. * (maxr - minr)
    rvmax2 = minr + 0.45 * (maxr - minr)
    rvol2 = mlab.pipeline.volume(rho_pert, vmin=rvmin2, vmax=rvmax2)
    
    # Changing the ctf:
    ctf2 = ColorTransferFunction()
    ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    ctf2.add_rgb_point(rvmin2 + 0.5 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    # ...
    rvol2._volume_property.set_color(ctf2)
    rvol2._ctf = ctf2
    rvol2.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(rvmax2, 0)
    otf.add_point(rvmin2, 0.10)
    ##vol1._otf = otf
    rvol2._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol2.volume_property.shade = False
    rvol2.volume_property.scalar_opacity_unit_distance = 2.0
    
if show_density_pert2 == True:
    # Scalar field density   
    rho_pert2 = mlab.pipeline.scalar_field(rho_vals_pert2, name="density", figure=fig)
    #scalar_cut_plane = ScalarCutPlane()
    #fig.parent.add_filter(scalar_cut_plane, sca)
    rho_pert2.spacing = spacing
    minr = rho_vals_pert2.min()
    maxr = rho_vals_pert2.max()
    
    #Volume for high pressure
    rvmin1 = minr + 0.55 * (maxr - minr)
    rvmax1 = minr + 1. * (maxr - minr)
    rvol1 = mlab.pipeline.volume(rho_pert2, vmin=rvmin1, vmax=rvmax1)
    
    # Changing the ctf:
    from tvtk.util.ctf import ColorTransferFunction
    ctf1 = ColorTransferFunction()
    ctf1.add_rgb_point(rvmin1, 1., 1., 0.5)
    ctf1.add_rgb_point(rvmin1 + 0.5 * (rvmax1 - rvmin1), 1, 0.3, 0.1)
    ctf1.add_rgb_point(rvmax1, 1., 0., 0.)
    # ...
    rvol1._volume_property.set_color(ctf1)
    rvol1._ctf = ctf1
    rvol1.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(rvmin1, 0)
    otf.add_point(rvmax1, 0.10)
    ##vol1._otf = otf
    rvol1._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol1.volume_property.shade = False
    rvol1.volume_property.scalar_opacity_unit_distance = 2.0
    
    
    #Volume for low pressure
    rvmin2 = minr + 0. * (maxr - minr)
    rvmax2 = minr + 0.45 * (maxr - minr)
    rvol2 = mlab.pipeline.volume(rho_pert2, vmin=rvmin2, vmax=rvmax2)
    
    # Changing the ctf:
    ctf2 = ColorTransferFunction()
    ctf2.add_rgb_point(rvmin2, 0., 0.5, 1.)
    ctf2.add_rgb_point(rvmin2 + 0.5 * (rvmax2 - rvmin2), 0.1, 0.7, 1.)
    ctf2.add_rgb_point(rvmax2, 0.5, 1., 1.)
    # ...
    rvol2._volume_property.set_color(ctf2)
    rvol2._ctf = ctf2
    rvol2.update_ctf = True
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(rvmax2, 0)
    otf.add_point(rvmin2, 0.10)
    ##vol1._otf = otf
    rvol2._volume_property.set_scalar_opacity(otf)
    
    # exempt volume from shading and improve overall look by increasing opacity
    rvol2.volume_property.shade = False
    rvol2.volume_property.scalar_opacity_unit_distance = 2.0

if show_pressure == True:
    # Scalar field p_tot    
    sca = mlab.pipeline.scalar_field(p_totvals, name="p_tot", figure=fig)
    #scalar_cut_plane = ScalarCutPlane()
    #fig.parent.add_filter(scalar_cut_plane, sca)
    sca.spacing = spacing
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
    otf = PiecewiseFunction()
    otf.add_point(pvmin1, 0)
    otf.add_point(pvmax1, 0.1)
    ##vol1._otf = otf
    vol1._volume_property.set_scalar_opacity(otf)
    
    
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
    
    #Changing the opacity of the volume vol1
    ## Changing the otf:
    otf = PiecewiseFunction()
    otf.add_point(pvmax2, 0)
    otf.add_point(pvmin2, 0.10)
    ##vol1._otf = otf
    vol2._volume_property.set_scalar_opacity(otf)



                                   
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
    

# Vector field bxvals, bzvals, byvals
field = mlab.pipeline.vector_field(bxvals, bzvals, byvals, name="B field", 
                                       figure=fig)
field.spacing = spacing
    
if show_mag == True:
    magnitude = mlab.pipeline.extract_vector_norm(field)
    #contours = mlab.pipeline.iso_surface(magnitude,
    #                                        contours=range(2, 14, 3),
    #                                        transparent=True,
    #                                        opacity=0.4,
    #                                        colormap='YlGnBu',
    #                                        vmin=0, vmax=14)
    
    # Create an array of points for which we want mag field seeds
    nx_seed = 9 #7
    ny_seed = 13 #10
    start_x = 30. #38
    end_x = nx+1 - start_x
    start_y = 1.
    end_y = ny-1 #ny-2 for ny = 100
    seeds=[]
    dx_res = (end_x - start_x) / (nx_seed-1)
    dy_res = (end_y - start_y) / (ny_seed-1)
    for i in range(0,nx_seed):
        for j in range(0,ny_seed):
            x = start_x + (i * dx_res) * x_spacing
            y = start_y + (j * dy_res) * y_spacing
            z = 1.
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
    

if show_mag_vec == True:
    bdirfield_front = mlab.pipeline.vector_field(bxvals_mask_front, bzvals_mask_front,
                                                 byvals_mask_front, name="B field front",
                                                 figure=fig)
    bdirfield_front.spacing = spacing
    vector_cut_plane_front = mlab.pipeline.vector_cut_plane(bdirfield_front, 
                                                      scale_factor=8.)
    vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
    vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
    vector_cut_plane_front.glyph.color_mode = 'no_coloring'
    vector_cut_plane_front.implicit_plane.widget.enabled = False
    vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
    vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'


#vdirfield = mlab.pipeline.vector_field(vxvals_mask, vzvals_mask, vyvals_mask,
#                                       name="V field", figure=fig)
#vectors = mlab.pipeline.vectors(vdirfield, scale_factor=15.)
##vectors.glyph.glyph.orient = False
#vectors.glyph.color_mode = 'no_coloring'

if show_vel_top == True or show_vel_top_pert == True:
    vdirfield_top = mlab.pipeline.vector_field(vxvals_mask_top, np.zeros_like(vxvals_mask_top),
                                                vyvals_mask_top, name="V field top",
                                                figure=fig)
    vdirfield_top.spacing = spacing
    vector_cut_plane_top = mlab.pipeline.vector_cut_plane(vdirfield_top, 
                                                      scale_factor=8.)
    vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
    vector_cut_plane_top.glyph.color_mode = 'no_coloring'
    vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-1, 50.5])
    vector_cut_plane_top.implicit_plane.widget.enabled = False
    vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
    vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
    
if show_vel_front == True or show_vel_front_pert == True:
    vdirfield_front = mlab.pipeline.vector_field(vxvals_mask_front, vzvals_mask_front,
                                                 vyvals_mask_front, name="V field front",
                                                 figure=fig)
    vdirfield_front.spacing = spacing
    vector_cut_plane_front = mlab.pipeline.vector_cut_plane(vdirfield_front, 
                                                      scale_factor=8.)
    vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
    vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
    vector_cut_plane_front.glyph.color_mode = 'no_coloring'
    vector_cut_plane_front.implicit_plane.widget.enabled = False
    vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
    vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'

if show_disp_top == True:
    xidirfield_top = mlab.pipeline.vector_field(xixvals_mask_top, np.zeros_like(xixvals_mask_top),
                                                xiyvals_mask_top, name="Xi field top",
                                                figure=fig)
    xidirfield_top.spacing = spacing
    vector_cut_plane_top = mlab.pipeline.vector_cut_plane(xidirfield_top, 
                                                      scale_factor=8.)
    vector_cut_plane_top.implicit_plane.widget.normal_to_y_axis = True
    vector_cut_plane_top.glyph.color_mode = 'no_coloring'
    vector_cut_plane_top.implicit_plane.widget.origin = np.array([ 50.,nz-1, 50.5])
    vector_cut_plane_top.implicit_plane.widget.enabled = False
    vector_cut_plane_top.glyph.glyph_source.glyph_source = vector_cut_plane_top.glyph.glyph_source.glyph_dict['arrow_source']
    vector_cut_plane_top.glyph.glyph_source.glyph_position = 'center'
    
if show_disp_front == True:
    xidirfield_front = mlab.pipeline.vector_field(xixvals_mask_front, xizvals_mask_front,
                                                 xiyvals_mask_front, name="Xi field front",
                                                 figure=fig)
    xidirfield_front.spacing = spacing
    vector_cut_plane_front = mlab.pipeline.vector_cut_plane(xidirfield_front, 
                                                      scale_factor=8.)
    vector_cut_plane_front.implicit_plane.widget.normal_to_z_axis = True
    vector_cut_plane_front.implicit_plane.widget.origin = np.array([ 50., 25.91140784, (ny-1)*y_spacing])
    vector_cut_plane_front.glyph.color_mode = 'no_coloring'
    vector_cut_plane_front.implicit_plane.widget.enabled = False
    vector_cut_plane_front.glyph.glyph_source.glyph_source = vector_cut_plane_front.glyph.glyph_source.glyph_dict['arrow_source']
    vector_cut_plane_front.glyph.glyph_source.glyph_position = 'center'

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
    field.scene.camera.view_angle = 17.0
    field.scene.camera.view_up = [-0.002418791139063777, 0.93281530024654913, -0.36034672896443193]
    field.scene.camera.clipping_range = [345.97885880654962, 650.71850659694883]
 
if view == 'front side':
    field.scene.camera.position = [219.75240522356913, 79.115530658939363, 499.350055117431]
    field.scene.camera.focal_point = [50.821544647216797, 50.413210511207581, 50.159849926829338]
    field.scene.camera.view_angle = 17.0
    field.scene.camera.view_up = [-0.044827368757468969, 0.99789300973432205, -0.046904670701947544]
    field.scene.camera.clipping_range = [348.0134523780402, 654.26836207642464]

if show_axes == True:
    axes = mlab.axes(field, nb_labels=1, line_width=3)
    axes.axes.x_label = ''
    axes.axes.y_label = ''
    axes.axes.z_label = ''
    axes.axes.label_format = ''

#if uniform_light == True:
#    #uniform lighting
#    field.scene.light_manager.number_of_lights = 6
#    
#    camera_light1 = field.scene.light_manager.lights[0]
#    camera_light1.activate = True
#    camera_light1.intensity = 0.7
#    camera_light1.elevation = 90.
#    camera_light1.azimuth = 0.
#
#    camera_light2 = field.scene.light_manager.lights[1]
#    camera_light2.activate = True
#    camera_light2.intensity = 0.7
#    camera_light2.elevation = -90.
#    camera_light2.azimuth = 0.
#
#    camera_light3 = field.scene.light_manager.lights[2]
#    camera_light3.activate = True
#    camera_light3.intensity = 0.7
#    camera_light3.elevation = 0.
#    camera_light3.azimuth = -90
#
#    camera_light4 = field.scene.light_manager.lights[3]
#    camera_light4.activate = True
#    camera_light4.intensity = 0.7
#    camera_light4.elevation = 0.
#    camera_light4.azimuth = 0.
#
#    camera_light5 = field.scene.light_manager.lights[4]
#    camera_light5.activate = True
#    camera_light5.intensity = 0.7
#    camera_light5.elevation = 0.
#    camera_light5.azimuth = 90.
#
#    camera_light6 = field.scene.light_manager.lights[5]
#    camera_light6.activate = True
#    camera_light6.intensity = 0.7
#    camera_light6.elevation = 0.
#    camera_light6.azimuth = 180.

#Black background
field.scene.background = (0., 0., 0.)

mlab.show()