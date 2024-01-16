# QG_spec
3D Quasi-geostrophic pseudo-spectral code

To compile:

gfortran -O3 -o qgps commons.f90  omp_stafft.f90 fft3d.f90 spectral.f90  init.f90 qgps.f90

qgps.f90 uses RK4 for time integration
qgps_semi.f90 is semi-implicit with a Leapfrog time integration

The current version uses static arrays, so use -mcmodel=medium for high resolution (defined in commons.f90)



