!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BPMclass: Beam position monitor class in Lattice module of APPLICATION 
!           layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the different beam diagnostics at given
!              beam position.
! Comments:
!  1) Itype = -1, shift the transverse centroid position to 0.
!  2) Itype = -2, shift the transverse centroid position and angle to 0.
!                 (this one not work yet due to conflict of definition)
!  3) Itype = -10, mismatch the beam distribution by the amount given in
!                  Param(3) - Param(10).  
!  4) Itype = -13, collimator slit
!  5) Itype = -21, shift the beam centroid in 6D phase space by the amount
!                  given in Param(3) - Param(10).
!----------------------------------------------------------------
      module BPMclass
        integer, private, parameter :: Nparam = 8
        type BPM
          !Itype < 0
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
          ! Param(1) : zedge
          !      (2) : radius
          !      (3) : xmax
          !      (4) : pxmax
          !      (5) : ymax
          !      (6) : pymax
          !      (7) : zmax
          !      (8) : pzmax
        end type BPM
        interface getparam_BPM
          module procedure getparam1_BPM,  &
                          getparam2_BPM,   &
                          getparam3_BPM
        end interface
        interface setparam_BPM
          module procedure setparam1_BPM,  &
                           setparam2_BPM, setparam3_BPM
        end interface
      contains
        subroutine construct_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_BPM
   
        subroutine setparam1_BPM(this,i,value)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_BPM

        subroutine setparam2_BPM(this,values)
        implicit none
        type (BPM), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_BPM

        subroutine setparam3_BPM(this,numseg,nmpstp,type,blength)
        implicit none
        type (BPM), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength
        
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_BPM
   
        subroutine getparam1_BPM(this,i,blparam) 
        implicit none 
        type (BPM), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_BPM
  
        subroutine getparam2_BPM(this,blparams)
        implicit none
        type (BPM), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_BPM

        subroutine getparam3_BPM(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (BPM), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_BPM

       !"-55" introduce an instant R56 to the beam
        subroutine chicane_BPM(Pts1,innp,gamma0,r56,t566,u5666,g0)
        implicit none
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: gamma0,r56,t566,u5666
        real*8 :: de,dzz
        integer :: i
        real*8 :: z0,g0,tmp1

!        t566 = -1.5d0*r56
!        u5666 = 2*r56
!----------------
! find the centroid and shift to the centroid.
        tmp1 = sum(Pts1(3,:))
        z0 = 0.0d0
        g0 = 0.0d0
!        do i = 1, innp
!          z0 = z0 + Pts1(1,i)*Pts1(3,i)/tmp1
!          g0 = g0 + Pts1(2,i)*Pts1(3,i)/tmp1
!        enddo
!        z0 = z0
!        g0 = g0
        Pts1(1,:) = Pts1(1,:) - z0
        Pts1(2,:) = Pts1(2,:) - g0
!-------------------
!        print*,"z0,g0: ",z0,g0

        do i = 1, innp
          de = Pts1(2,i)/gamma0
          dzz = r56*de+t566*de**2+u5666*de**3
          Pts1(1,i) = Pts1(1,i) + dzz
        enddo

        end subroutine chicane_BPM

        subroutine chicanebkwd_BPM(Pts1,innp,gamma0,r56,t566,u5666,g0)
        implicit none
        integer, intent(in) :: innp
        double precision, pointer, dimension(:,:) :: Pts1
        double precision, intent(in) :: gamma0,r56,t566,u5666
        real*8 :: de,dzz
        integer :: i
        real*8 :: z0,g0,tmp1

!        t566 = -1.5d0*r56
!        u5666 = 2*r56
!----------------
! find the centroid and shift to the centroid.
        tmp1 = sum(Pts1(3,:))
        z0 = 0.0d0
        g0 = 0.0d0
!        do i = 1, innp
!          z0 = z0 + Pts1(1,i)*Pts1(3,i)/tmp1
!          g0 = g0 + Pts1(2,i)*Pts1(3,i)/tmp1
!        enddo
!        z0 = z0
!        g0 = g0
        Pts1(1,:) = Pts1(1,:) - z0
        Pts1(2,:) = Pts1(2,:) - g0
!-------------------
!        print*,"z0,g0: ",z0,g0

        do i = 1, innp
          de = Pts1(2,i)/gamma0
          dzz = r56*de+t566*de**2+u5666*de**3
          Pts1(1,i) = Pts1(1,i) - dzz
        enddo

        end subroutine chicanebkwd_BPM

      end module BPMclass
