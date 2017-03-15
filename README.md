# mhd-slab-3d
3D visualisations of MHD waves in asymmetric an asymmetric slab

Prerequisites:
Python 2 (https://www.python.org/download/releases/2.7.2/),
Mayavi (http://mayavi.sourceforge.net/).

3D_slab_modes.py: Builds an animation of the magnetoacoustic eigenmodes of an isolated asymmetric magnetic slab with possible visualisations of:
magnetic field lines,
magentic field strength,
magnetic field vectors,
velocity (eularian and lagrangian),
displacement,
density (eularian and lagrangian),
slab boundary perturbations.
Also, can produce dispersion diagrams.

alfven_mixed_driver.py: Builds an animation of a shear alfven wave in a plasma with spacially variable alfven speed (otherwise uniform m plasma

We use the analytical solutions given in our paper: http://link.springer.com/article/10.1007/s11207-017-1054-y
