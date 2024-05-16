!----------------------------------------------------------------
!EBLT -  Electron Beam Longitudinal Tracking (forward and backward)
!*** Copyright Notice ***
!
!Electron Beam Longitudinal Tracking  (EBLT) Copyright (c) 2024, The
!Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
!
!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Intellectual Property Office at
!IPO@lbl.gov.
!NOTICE.  This Software was developed under funding from the U.S. Department
!of Energy and the U.S. Government consequently retains certain rights.  As
!such, the U.S. Government has been granted for itself and others acting on
!its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
!Software to reproduce, distribute copies to the public, prepare derivative 
!works, and perform publicly and display publicly, and to permit others to do so.
!****************************
! Description: v.1.0
! Comments: This is a parallel and serial code for electron beam longitudinal
!           forward or backward tracking. To use it as a parallel code, please
!           comment out the "use mpistub" in BeamBunch.f90 and use the Makefile_parallel.
!----------------------------------------------------------------
      program main
      use AccSimulatorclass
      implicit none

      call construct_AccSimulator()
      call run_AccSimulator()
      call destruct_AccSimulator()

      end program main
