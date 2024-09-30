!-------------------------------------------
!BELT - BEam Longitudinal Tracking (forward or backward)
!-------------------------------------------
*** Copyright Notice ***

Electron Beam Longitudinal Tracking  (EBLT) Copyright (c) 2024, The
Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.
!-------------------------
!Author: Ji Qiang, LBNL
!-------------------------
!This is a serial/parallel code for fast electron beam longitudinal
!forward or backward tracking through an electron linear accelerator. 
!It includes longitudinal space-charge effect, structure 
!and resistive wall wakefields, longitudinal coherent synchrotron radiation effect.
!To use it as a parallel code, please
!comment out the "use mpistub" in BeamBunch.f90, rename mpif.h as mpif.hh,
!and use the Makefile_parallel as Makefile.
!Each macroparticle has three attributes: delta z (m), delta gamma, and weight (total charge/Np)
!The physical model is in: 
!J. Qiang, "Fast longitudinal beam dynamics optimization in x-ray free electron laser 
!linear accelerators," Phys. Rev. Accel. Beams 22, 094401 (2019).

!----------------
!Input file: eblt.in
!-----------------------------------
!"1024" # of longitudinal particles for beam distribution, 
!"256 # of grid points for collective effect (LSC,CSR,wakefield) calculation
!"-0.7910726E-02" z_min (m) and "0.7105112E-02" z_max (m)
!"1" switch flag for forward tracking, otherwise backward tracking.
1024 256 -0.7910726E-02  0.7105112E-02 1 200
!--------
!"-200" switch for input distribution:
!1-use the polynominal coefficients of current and longitudinal phase space,
!2-use the Gaussian current profile and poly. longitudinal phase space,
!100.0-readin from EBLT particle output,
!200.0-readin from ImpactZ slice output,
!300.0-readin from ImpactZ particle output,
!those are polynominal coefficients a0 to a9 for current (A) profile
0.0 0.0 0.0 /
!---------
!polynominal coefficients b0 to b9 delta gamma-z profile
1.0 0.0 0.0  /
!--
!"0.13" current (I=Q*freq), "0.92706472e+08" e beam energy, 
!"0.511005d+06" e mass (eV), "-1.0" e charge, "1.3d9" reference RF freq.
0.13 0.92706472e+08 0.511005d+06 -1.0 1.3d9
!--------------------------------------------------------------------
!lattice input infor.
!length of ele., # of steps, # of maps (not used in most elements), 
!element type code ("0" for drift), beam trans. radius in SC model (m)
74.9744159418 2 10 0 0.3e-3 /
!
!"-41" longitudinal wake field, "1.0" scaling of wakefunction,
!"41" wake function stored in file rfdata41, "1.0" turn-on swith 
0.0 0 1 -41 1.0 41 1.0 /
!
!lumped RF cavity model
!length, "2" # of steps, "1" not used, "103" RF type code, "0.3e-3" trans. radius,
!"0.1576186516E+08" accel. gradient (V/m), "1.3e9" rf freq., "-1.01" not used, "1.0" not used.
!"-0.1789850476E+02" rf. phase, 
16.603888 2 1 103 0.3e-3 0.1576186516E+08 1.3e9 -0.1789850476E+02 -1.01 1.0 /
!
!"-41" longitudinal wake field, "1.0" scaling of wakefunction,
!"41" wake function stored in file rfdata41, "-1.0" turn-off swith 
0.0 0 1 -41 1.0 41 -1.0 /
!
!user specified current profile output in fort.50 and particle distribution in fort.51
0.0 0 50 -2 1.0 /
!
!length, "10" # of steps, if # of maps is "1" calculate r56, t566, u5666 for 4 dipole C chicane using
!bending angle "0.09962502193", drift between b1-2, b3-4 "2.44781923989" and 
!bend length "0.203558243263",
!if # of maps is "-1", "2.44781923989"->r56, "-0.0773102006514754"->t566, "0.103080267535301"->u5666. "1.01"-> IGF CSR including A and B, "2.01" -> IGF steady state CSR"0.0"-> no CSR.
0.203558243263 10 1 4 0.3e-3 2.44781923989 -0.0773102006514754 0.103080267535301 0.09962502193 1.01 0.0 / 1.005x

!--------------------------------------------------------
!Output files:
!-------------------------------------------------------
!fort.2: 3 columns - distance, kin. energy, gamma , <dz>, rms dz, <delta gamma>, rms delta gamma
!fort.100: initial current profile (3 columns) - bunch length (m), charge per cell, current (A)
!fort.101: initial particle distribution (4 columns) - z coordinate (m), del_gamma, ptcl. weight, dE/E0
!fort.50: user specified (-2) current profile (3 columns) - bunch length (m), charge per cell, current (A)
!fort.51: user specified (-2) particle distribution (4 columns) - z coordinate (m), del_gamma, ptcl. weight, dE/E0
!fort.200: final current profile (3 columns) - bunch length (m), charge per cell, current (A)
!fort.201: final particle distribution (4 columns) - z coordinate (m), del_gamma, ptcl. weight, dE/E0
