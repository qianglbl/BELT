# EBLT: Electron Beam Longitudinal Tracking (forward or backward)
*** Copyright Notice ***

Electron Beam Longitudinal Tracking  (EBLT) Copyright (c) 2024, The
Regents of the University of California, through Lawrence Berkeley National Labo
ratory (subject to receipt of any required approvals from the U.S. Dept. of Ener
gy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit others to do so.

-------------------------
Contact: Ji Qiang (jqiang@lbl.gov), LBNL
-------------------------

This is a serial/parallel code for fast electron beam longitudinal
forward or backward tracking through an electron linear accelerator.
It includes longitudinal space-charge effect, structure
and resistive wall wakefields, longitudinal coherent synchrotron radiation effect.
To use it as a parallel code, please
comment out the "use mpistub" in BeamBunch.f90, rename mpif.h as mpif.hh,
and use the Makefile_parallel as Makefile.
Each macroparticle has three attributes: delta z (m), delta gamma, and weight (total charge/Np)
The physical model is in:

------------------------
Citation:

J. Qiang, "Fast longitudinal beam dynamics optimization in x-ray free electron laser
linear accelerators," Phys. Rev. Accel. Beams 22, 094401 (2019).

