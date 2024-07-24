!----------------------------------------------------------------
!EBLT -  Electron Beam Longitudinal Tracking (forward and backward)
!*** Copyright Notice ***

!Electron Beam Longitudinal Tracking  (EBLT) Copyright (c) 2024, The
!Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

!If you have questions about your rights to use or distribute this software,
!please contact Berkeley Lab's Intellectual Property Office at
!IPO@lbl.gov.

!NOTICE.  This Software was developed under funding from the U.S. Department
!of Energy and the U.S. Government consequently retains certain rights.  As
!such, the U.S. Government has been granted for itself and others acting on
!its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
!Software to reproduce, distribute copies to the public, prepare derivative 
!works, and perform publicly and display publicly, and to permit others to do so.
!-------------------------
! AccSimulatorclass: Linear accelerator simulator class.
! Author: Ji Qiang, LBNL
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments:
!----------------------------------------------------------------
      module AccSimulatorclass
        use BeamBunchclass
        use Inputclass
        use BPMclass
        implicit none

        !# of num. total x, total and local y mesh pts., type of BC, 
        !# of beam elems, type of integrator.
        integer, private :: Nz,Nblem,Np,Nplc
        integer, private :: nproc,myid

        !communicator used in the eblo
        integer, private :: commeblt

        !beam current, kin. energy, part. mass, and charge.
        double precision, private :: Bcurr,Bkenergy,Bmass,Bcharge,Bfreq,&
                                     Perdlen

        !conts. in init. dist.
        double precision, private, dimension(20) :: distparam

        !beam particle object and array.
        type (BeamBunch), private :: Bpts

        !beam line element array.
        type (BPM),target,dimension(20000) :: beamln

        real*8 :: clite

        integer :: flagfwd
        real*8 :: Bchgtot

        !beam line element period.
        interface construct_AccSimulator
          module procedure init_AccSimulator
        end interface

      contains
        !set up objects and parameters.
        subroutine init_AccSimulator()
        implicit none
        include 'mpif.h'
        integer :: i,j
        integer, allocatable, dimension(:) :: bnseg,bmpstp,bitype
        double precision, allocatable, dimension(:) :: blength,val1,&
        val2, val3,val4,val5,val6,val7,val8
        double precision :: z,zmin,zmax
        double precision, dimension(8) :: tmpbpm
        double precision, dimension(9) :: tmp1
        integer :: ibpm,isamp
        real*8 :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
        real*8 :: b0,b1,b2,b3,b4,b5,b6,b7,b8,b9
        real*8 :: hz,phsini,gam0,gam0tmp
        integer :: ierr,npl,npu,ii
        real*8 :: Pi,scxl,beta0,chglc,sqpi

        call MPI_INIT(ierr)

        commeblt = MPI_COMM_WORLD
        call MPI_COMM_SIZE(commeblt,nproc,ierr)
        call MPI_COMM_RANK(commeblt,myid,ierr)

        flagfwd = 1
        print*,"EBLT - Electron Beam Longitudinal Tracking (forward or backward): Vs1.0beta"
        clite = 299792458.0d0
!-------------------------------------------------------------------
! get all global input parameters.
        distparam = 0.0d0
        call in_Input(Np,Nz,distparam,20,Bcurr,Bkenergy,Bmass,Bcharge,&
        Bfreq,zmin,zmax,Nblem,flagfwd)
 
        Nplc = Np/nproc
        if(mod(Np,nproc).ne.0) then
          print*,"# of particles not divisible by # of processors!!!!"
          stop
        endif
        npl = myid*Nplc+1
        npu = (myid+1)*Nplc

        print*,"flagfwd:",flagfwd
        if(flagfwd.eq.1) then
          print*,"forward tracking..."
        else
          print*,"backward tracking..."
        endif

!-------------------------------------------------------------------
! construct BeamBunch class.
        phsini = 0.0
        call construct_BeamBunch(Bpts,Bcurr,Bkenergy,Bmass,Bcharge,&
                            Np,Nplc,phsini)

!-------------------------------------------------------------------
! sample initial particle distribution.

        gam0 = Bkenergy/Bmass+1.0d0
        beta0 = sqrt(1.0d0-1.0d0/gam0**2)
        Pi = 2*asin(1.0d0)
        sqpi = sqrt(2*Pi)
        scxl = clite*beta0/(2*Pi*Bfreq)

        allocate(Bpts%Pts1(3,Nplc))

        hz = (zmax-zmin)/(Np-1)
        
        a0 = distparam(1)
        a1 = distparam(2)
        a2 = distparam(3)
        a3 = distparam(4)
        a4 = distparam(5)
        a5 = distparam(6)
        a6 = distparam(7)
        a7 = distparam(8)
        a8 = distparam(9)
        a9 = distparam(10)
        b0 = distparam(11)
        b1 = distparam(12)
        b2 = distparam(13)
        b3 = distparam(14)
        b4 = distparam(15)
        b5 = distparam(16)
        b6 = distparam(17)
        b7 = distparam(18)
        b8 = distparam(19)
        b9 = distparam(20)

        if(a0.lt.0) then
!          isamp = 0
          isamp = abs(a0)+0.001 
        else
          isamp = a0
        endif
        print*,"isamp: ",isamp
        !print*,"hz:",hz,clite,a0

        if(isamp.eq.1) then
          do i = 1, Np
            z = zmin + (i-1)*hz
            if((i.ge.npl).and.(i.le.npu)) then
              ii = i - myid*Nplc
              Bpts%Pts1(1,ii) = z 
              Bpts%Pts1(2,ii) = b0+b1*z+b2*z**2+b3*z**3+b4*z**4+b5*z**5+&
                           b6*z**6+b7*z**7+b8*z**8+b9*z**9
              Bpts%Pts1(3,ii) = (a0+a1*z+a2*z**2+a3*z**3+a4*z**4+a5*z**5+&
                           a6*z**6+a7*z**7+a8*z**8+a9*z**9)*hz/clite
            endif
          enddo
        
        else if(isamp.eq.2) then
          do i = 1, Np
            z = zmin + (i-1)*hz
            if((i.ge.npl).and.(i.le.npu)) then
              ii = i - myid*Nplc
              Bpts%Pts1(1,ii) = z
              Bpts%Pts1(2,ii) = b0+b1*z+b2*z**2+b3*z**3+b4*z**4+b5*z**5+&
                           b6*z**6+b7*z**7+b8*z**8+b9*z**9
              Bpts%Pts1(3,ii) = a0*(a1*exp(-((z-a2)/a3)**2/2)/a3+&
                      a4*exp(-((z-a5)/a6)**2/2)/a6 + &
                      a7*exp(-((z-a8)/a9)**2/2)/a9)*hz/clite/sqpi
            endif
          enddo
        else if(isamp.eq.100) then !read in from the EBLT output particle dist.
          open(1,file="pts.in",status="old") 
          do i = 1, Np
            read(1,*)tmp1(1:4)
            if((i.ge.npl).and.(i.le.npu)) then
              ii = i - myid*Nplc
              Bpts%Pts1(1,ii) = tmp1(1)
              Bpts%Pts1(2,ii) = tmp1(2)
              Bpts%Pts1(3,ii) = tmp1(3)
            endif
          enddo
          close(1)
        else if(isamp.eq.200) then !readin
          open(1,file="pts.in",status="old") !read in from the Impact slice output file
          do i = 1, Np
            read(1,*)tmp1(1:6)
            if((i.ge.npl).and.(i.le.npu)) then
              ii = i - myid*Nplc
              Bpts%Pts1(1,ii) = tmp1(1)
              Bpts%Pts1(2,ii) = tmp1(6)*gam0
              Bpts%Pts1(3,ii) = tmp1(3)/clite*hz
            endif
          enddo
          close(1)
        else if(isamp.eq.300) then !readin
          open(1,file="pts.in",status="old") !read in from the Impact particle output file
          do i = 1, Np
            read(1,*)tmp1(1:9)
            if((i.ge.npl).and.(i.le.npu)) then
              ii = i - myid*Nplc
              Bpts%Pts1(1,ii) = -tmp1(5)*scxl
              Bpts%Pts1(2,ii) = -tmp1(6)
              Bpts%Pts1(3,ii) = abs(tmp1(8))
            endif
          enddo
          close(1)
        else
          print*,"wrong initial distribution!!"
          stop
        endif


        chglc = sum(Bpts%Pts1(3,:))

        print*,"chglc:",chglc,a0,b0
        call MPI_ALLREDUCE(chglc,Bchgtot,1,MPI_DOUBLE_PRECISION,&
             MPI_SUM,commeblt,ierr)

        if(abs(Bchgtot).lt.1.0d-16) then
            Bchgtot = 1.0d-17
            Bpts%Pts1(3,:) = Bchgtot/Np
        endif
        print*,"total charge (C):",Bchgtot

!-------------------------------------------------------------------
! construct beam line elements.
        allocate(blength(Nblem),bnseg(Nblem),bmpstp(Nblem))
        allocate(bitype(Nblem))
        allocate(val1(Nblem),val2(Nblem),val3(Nblem),val4(Nblem))
        allocate(val5(Nblem),val6(Nblem),val7(Nblem),val8(Nblem))

        call in_Input(Nblem,blength,bnseg,bmpstp,bitype,val1,val2,val3,&
        val4,val5,val6,val7,val8)

        ibpm = 0
        do i = 1, Nblem
            ibpm = ibpm + 1
            call construct_BPM(beamln(ibpm),bnseg(i),bmpstp(i),&
                 bitype(i),blength(i))
            tmpbpm(1) = 0.0
            tmpbpm(2) = val1(i)
            tmpbpm(3) = val2(i)
            tmpbpm(4) = val3(i)
            tmpbpm(5) = val4(i)
            tmpbpm(6) = val5(i)
            tmpbpm(7) = val6(i)
            tmpbpm(8) = val7(i)
            call setparam_BPM(beamln(ibpm),tmpbpm)
        enddo
!-------------------------------------------------------------------
        print*,"pass setting up lattice..."

        deallocate(blength,bnseg,bmpstp,bitype)
        deallocate(val1,val2,val3,val4,val5,val6,val7,val8)
        call MPI_BARRIER(commeblt,ierr)

        end subroutine init_AccSimulator

        !Run beam dynamics simulation through accelerator.
        subroutine run_AccSimulator()
        implicit none
        include 'mpif.h'
        integer, parameter :: nwkmax = 10000
        integer :: i,j,bnseg,bmpstp,bitype,nstep,bmap
        double precision :: z,tau1,tau2,blength
        double precision :: hz,zmin,zedge,radius,rfile
        double precision :: zminpt,zmaxpt,zmax
        !parameters for stripper modeling
        double precision :: gamma0,r56,t566,u5666,bangle
        double precision, allocatable, dimension(:) :: &
            densz,ezwake,ezlsc
        integer :: flagwake,flagcsr,flagsc
        real*8 :: zwkmin,bendlen,zbleng,scwk
        integer :: ihlf,ifile,ndatawk
        real*8, dimension(nwkmax) :: wklong
        real*8 :: tmpwk,tmp1,lenwk,r0,deleng,g0
        real*8 :: driftbtbd
        integer :: sample
        real*8 :: zmingl,zmaxgl
        integer :: ierr,ipt
        real*8, dimension(4) :: ztmplc,ztmpgl
        real*8 :: zavg,zsig,deavg,desig

!-------------------------------------------------------------------
! prepare initial parameters, allocate temporary array.
        scwk = 1.0d0

!-------------------------------------------------------------------
! prepare for round pipe, open longitudinal

        !assign initial storage used for wakefield calculation
        allocate(densz(Nz))
        allocate(ezlsc(Nz))
        allocate(ezwake(Nz))
        ezlsc = 0.0d0
        ezwake = 0.0d0
        z = 0.0d0

        gamma0 = Bpts%refptcl(6)
        sample = 1
        call output(myid,commeblt,nproc,Nplc,sample,Nz,100,Bpts%Pts1,gamma0)
        ztmplc(1) = 0.0d0
        ztmplc(2) = 0.0d0
        ztmplc(3) = 0.0d0
        ztmplc(4) = 0.0d0
        do ipt = 1, Nplc
          ztmplc(1) = ztmplc(1) + Bpts%Pts1(1,ipt)*Bpts%Pts1(3,ipt)/Bchgtot
          ztmplc(2) = ztmplc(2) + Bpts%Pts1(1,ipt)**2*Bpts%Pts1(3,ipt)/Bchgtot
          ztmplc(3) = ztmplc(3) + Bpts%Pts1(2,ipt)*Bpts%Pts1(3,ipt)/Bchgtot
          ztmplc(4) = ztmplc(4) + Bpts%Pts1(2,ipt)**2*Bpts%Pts1(3,ipt)/Bchgtot
        enddo
        call MPI_ALLREDUCE(ztmplc,ztmpgl,4,MPI_DOUBLE_PRECISION,&
                    MPI_SUM,commeblt,ierr)
        zavg = ztmpgl(1)
        zsig = sqrt(ztmpgl(2)-zavg**2)
        deavg = ztmpgl(3)
        desig = sqrt(ztmpgl(4)-deavg**2)

        if(myid.eq.0) then
        write(2,1011)z,(Bpts%refptcl(6)-1)*Bmass,Bpts%refptcl(6),zavg,zsig,deavg,desig
        endif

        flagcsr = 0
        flagsc  = 1
        flagwake = 0
        zbleng = 0.0d0
        nstep = 0
!-------------------------------------------------------------------
! start looping through 'Nblem' beam line elements.
        do i = 1, Nblem

          bitype = beamln(i)%Itype
          blength = beamln(i)%Length
          bnseg = beamln(i)%Nseg
          bmap = beamln(i)%Mapstp

          if(myid.eq.0) print*,"enter elment (type code): ",i,bitype

          tau1 = 0.0d0
          if(bitype.ge.0) tau1 = 0.5d0*blength/bnseg
          tau2 = 2.0d0*tau1

          gamma0 = Bpts%refptcl(6)
          if(bitype.eq.-2)then
            ifile = beamln(i)%Mapstp
            sample = beamln(i)%Param(2)
            if(sample.eq.0) sample = 1
            call output(myid,commeblt,nproc,Nplc,sample,Nz,ifile,Bpts%Pts1,gamma0)
          endif
          if(bitype.eq.-39)then
            deleng = beamln(i)%Param(2)
            Bpts%refptcl(6) = Bpts%refptcl(6) + deleng/Bmass
!            print*,"refpt eng:",(Bpts%refptcl(6)-1.d0)*Bmass
          endif
          if(bitype.eq.-41)then
            tmpwk = beamln(i)%Param(2)
            rfile = beamln(i)%Param(3)
            tmp1 = beamln(i)%Param(4)
            if(rfile.gt.0.0d0) then
              ifile = int(rfile + 0.1)
              call read1wk_Data(ifile,ndatawk,lenwk,wklong)
              if(ndatawk.gt.nwkmax) then
                print*,"too many data points in readin wakefield:"
                stop
              endif
            endif
            !print*,"-41flagbc: ",tmpwk,rfile,tmp1
            if(tmp1.gt.0.0d0) then !turn on the read-in wakefield.
              flagwake = 1
              scwk = tmpwk
            else !turn off
              flagwake = 0
              scwk = 1.0d0
            endif
          endif
          if(bitype.eq.-99) then
            exit
          endif

          zedge = z
          Beamln(i)%Param(1) = zedge
          if(myid.eq.0) print*,"zedge: ",zedge
          radius = Beamln(i)%Param(2)

          !/bend using Transport transfer map
          if(bitype.eq.4) then
              gamma0 = Bpts%refptcl(6)
!              r56 = beamln(i)%Param(3)
!              t566 = beamln(i)%Param(4)
!              u5666 = beamln(i)%Param(5)
! one can calculate those parameters from bending angle, bend length,
! and drift between 1 and 2 (3 and 4).
              bangle = beamln(i)%Param(6)
              if(bmap.gt.0) then
                driftbtbd = beamln(i)%Param(3)
                r56 = 2*bangle**2*(driftbtbd+2.0d0/3*abs(blength))*1.005
                t566 = -1.5d0*r56
                u5666 = 2*r56
              else
                r56 = beamln(i)%Param(3)
                t566 = beamln(i)%Param(4)
                u5666 = beamln(i)%Param(5)
              endif
              flagcsr = int(beamln(i)%Param(7))
              flagsc = int(beamln(i)%Param(8))
              if(flagfwd.eq.1) then !forward
                call chicane_BPM(Bpts%Pts1,Nplc,gamma0,r56,t566,u5666,g0)
              endif
              r0 = abs(blength)/bangle
              !this is due longitudinal shift
              Bpts%refptcl(6) = Bpts%refptcl(6) + g0
          else
             flagsc = 1
             flagcsr = 0
          endif
!-------------------------------------------------------------------
! loop through 'bnseg' numerical segments in each beam element
! using 2 step symplectic integeration (ie. leap frog).
          ihlf = 0
          do j = 1, bnseg
!-------------------------------------------------------------------
! use linear map or nonlinear Lorentz integrator to advance particles.
              ! spatial drift.
              !linear map integrator
             call map1(Bpts,beamln(i),z,tau1,bitype,&
                                    bnseg,j,ihlf,flagfwd)
             gamma0 = Bpts%refptcl(6)
!-------------------------------------------------------------------
!-------------------------------------------------------------------
             zminpt = minval(Bpts%Pts1(1,:))
             zmaxpt = maxval(Bpts%Pts1(1,:))

             call MPI_ALLREDUCE(zminpt,zmingl,1,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,commeblt,ierr)
             call MPI_ALLREDUCE(zmaxpt,zmaxgl,1,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,commeblt,ierr)
             zmin = zmingl - 0.1d0*(zmaxgl-zmingl)/Nz
             zmax = zmaxgl + 0.1d0*(zmaxgl-zmingl)/Nz
             
             hz = (zmax-zmin)/(Nz-1)

!             print*,"zmin,zmax",zmin,zmax,flagsc,flagwake,flagcsr
             
             ! deposit particles onto grid to obtain charge density.
             call deposit(Nplc,Nz,Bpts%Pts1,densz,hz,zmin,commeblt)
!             if(myid.eq.0) print*,"sum densz: ",sum(densz),hz
             ezlsc = 0.0d0
             if(flagsc.eq.1) then
                call Lsc(Nz,densz,ezlsc,hz,radius,gamma0)
             endif
!             if(myid.eq.0) print*,"sum ezlsc: ",sum(ezlsc),hz

             ezwake = 0.0d0
!-------------------------------------------------------------------------
             if(flagwake.eq.1) then
                !includes wakefield
                call wakefieldread(Nz,densz,ezwake,hz,lenwk,ndatawk,wklong)
                ezwake = scwk*ezwake
             endif

             !if(flagcsr.eq.1) then
             if(flagcsr.ge.1) then
                if(bitype.eq.4) then
                   !zwkmin = zminpt + (z-zbleng)
                   zwkmin = zmingl + (z-zbleng)
                   bendlen = abs(blength) !inside the bend
                else
                  print*,"wrong csr element"
                endif

                ezwake = 0.0d0
                if(flagcsr.eq.1) then
!using IGF for the s-s csr wake
                  call csrwakeTrIGF_FieldQuant(Nz,r0,zwkmin,hz,&
                              bendlen,densz,gamma0,ezwake)
                else if(flagcsr.eq.2) then
                  call csrwakeSS_FieldQuant(Nz,r0,hz,&
                              densz,gamma0,ezwake)
                else if(flagcsr.eq.3) then
                  call csrwakeSS2_FieldQuant(Nz,r0,hz,&
                              densz,gamma0,ezwake)
                endif
             endif

             ezlsc = ezlsc + ezwake
!             ezlsc = ezwake
!             if(myid.eq.0) print*,"sum ezlsc2: ",sum(ezlsc)

            call scatter(Nplc,Nz,Bpts%Pts1,ezlsc,hz,zmin,tau2,Bmass)

            call map1(Bpts,beamln(i),z,tau1,bitype,&
                                    bnseg,j,ihlf,flagfwd)

            nstep = nstep + 1

            ztmplc(1) = 0.0d0
            ztmplc(2) = 0.0d0
            ztmplc(3) = 0.0d0
            ztmplc(4) = 0.0d0
            do ipt = 1, Nplc
              ztmplc(1) = ztmplc(1) + Bpts%Pts1(1,ipt)*Bpts%Pts1(3,ipt)/Bchgtot
              ztmplc(2) = ztmplc(2) + Bpts%Pts1(1,ipt)**2*Bpts%Pts1(3,ipt)/Bchgtot
              ztmplc(3) = ztmplc(3) + Bpts%Pts1(2,ipt)*Bpts%Pts1(3,ipt)/Bchgtot
              ztmplc(4) = ztmplc(4) + Bpts%Pts1(2,ipt)**2*Bpts%Pts1(3,ipt)/Bchgtot
            enddo 
            call MPI_ALLREDUCE(ztmplc,ztmpgl,4,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,commeblt,ierr)
            zavg = ztmpgl(1)
            zsig = sqrt(ztmpgl(2)-zavg**2)
            deavg = ztmpgl(3)
            desig = sqrt(ztmpgl(4)-deavg**2)

            if(myid.eq.0) then
              write(2,1011)z,(Bpts%refptcl(6)-1)*Bmass,Bpts%refptcl(6),zavg,zsig,deavg,desig
            endif
            gamma0 = Bpts%refptcl(6)
          end do
!end loop bnseg
!---------------------------------------------------------------------
          if(bitype.eq.4 .and. flagfwd.ne.1) then !backward
            call chicanebkwd_BPM(Bpts%Pts1,Nplc,gamma0,r56,t566,u5666,g0)
          endif

          zbleng = zbleng + blength
        enddo
!end loop through nbeam line element
!------------------------------------------------

! final output.
        sample = 1
        call output(myid,commeblt,nproc,Nplc,sample,Nz,200,Bpts%Pts1,gamma0)

        deallocate(densz)
        deallocate(ezlsc)
        deallocate(ezwake)
1011    format(7(1x,e18.9))


        end subroutine run_AccSimulator

        subroutine destruct_AccSimulator()
        implicit none
        include "mpif.h"
        double precision :: time
        integer :: ierr
 
        deallocate(Bpts%Pts1)
        call MPI_Finalize(ierr)

        end subroutine destruct_AccSimulator

        !read in the discrete wake field data in z, x, y.
        !distribution along axis zdat from files "rfdatax or rfdataxx 
        !or rfdataxxx".
        subroutine read1wk_Data(ifile,nwake,lenwk,wkfld)
        implicit none
        integer, intent(in) :: ifile
        integer, intent(out) :: nwake
        real*8, intent(out) :: lenwk
        real*8, intent(inout), dimension(:) :: wkfld
        integer :: myrank,ierr,i,ii,jj,kk,ll,n
        double precision :: tmp1,tmp2,tmp3,tmp4,zdat1
        character*7 name1
        character*8 name2
        character*9 name3

        name1 = 'rfdatax'
        name2 = 'rfdataxx'
        name3 = 'rfdataxxx'

          if((ifile.ge.1).and.(ifile.le.9)) then
            name1(7:7) = char(ifile+48)
            open(14,file=name1,status='old')
!            open(15,file=name1//"out",status='unknown')
          else if((ifile.ge.10).and.(ifile.le.99)) then
            ii = ifile/10
            jj = ifile - ii*10
            name2(7:7) = char(ii+48)
            name2(8:8) = char(jj+48)
            open(14,file=name2,status='old')
!            open(15,file=name2//"out",status='unknown')
          else if((ifile.ge.100).and.(ifile.le.999)) then
            ii = ifile/100
            jj = ifile - 100*ii
            kk = jj/10
            ll = jj - 10*kk
            name3(7:7) = char(ii+48)
            name3(8:8) = char(kk+48)
            name3(9:9) = char(ll+48)
            open(14,file=name3,status='old')
!            open(15,file=name3//"out",status='unknown')
          else
            print*,"out of the range of maximum 999 files!!!!"
          endif

          n = 0
50        continue
            read(14,*,end=77)tmp1,tmp2,tmp3,tmp4
            n = n + 1
            wkfld(n) = tmp2
          goto 50
77        continue
          close(14)
          nwake = n
          lenwk = tmp1

        end subroutine read1wk_Data
  
        ! rays(1,:) -> z
        ! rays(2,:) -> delta E
        ! rays(3,:) -> weight (C)
        ! deposit particles onto grid.
        subroutine output(myid,commin,npr,innplc,samplePeriod,innz,ifile,rays,gamma0)
        implicit none
        include "mpif.h"
        integer, intent(in) :: myid,commin,npr
        integer, intent(in) :: innplc,innz,ifile,samplePeriod
        double precision, intent (in), dimension (3,innplc) :: rays
        double precision, dimension (innz) :: rho,rhogl
        real*8, intent(in) :: gamma0
        integer :: kx,kx1,i,np3,j,ierr
        double precision :: ef
        double precision :: hz,hzi,zminpt,zmaxpt,zmin,zmax,zmingl,zmaxgl
        real*8, allocatable, dimension(:,:) :: recvbuf
        integer status(MPI_STATUS_SIZE)
        integer :: ifile1


        zminpt = minval(rays(1,:))
        zmaxpt = maxval(rays(1,:))
        call MPI_ALLREDUCE(zminpt,zmingl,1,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,commin,ierr)
        call MPI_ALLREDUCE(zmaxpt,zmaxgl,1,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,commin,ierr)
        zmin = zmingl - 0.1d0*(zmaxgl-zmingl)/innz
        zmax = zmaxgl + 0.1d0*(zmaxgl-zmingl)/innz

        hz = (zmax-zmin)/(innz-1)

        !print*,"zmin2: ",zminpt,zmaxpt,zmin,zmax,hz,innz,innp

        hzi = 1.0d0/hz

        rho=0.0d0
        rhogl=0.0d0
        do i = 1, innplc
          kx=(rays(1,i)-zmin)*hzi + 1 
          ef=((zmin-rays(1,i))+kx*hz)*hzi
          kx1=kx+1

          rho(kx) = rho(kx) + ef*rays(3,i)
          rho(kx1) = rho(kx1)+(1.0d0-ef)*rays(3,i)
        enddo

        call MPI_ALLREDUCE(rho,rhogl,innz,MPI_DOUBLE_PRECISION,&
                        MPI_SUM,commin,ierr)
        rho = rhogl
        rho = rho/hz

        if(myid.eq.0) then
          do i = 1, innz
            write(ifile,110)zmin+(i-1)*hz,rho(i),rho(i)*clite
          enddo
          close(ifile)
        endif
110     format(3(1x,e21.12))

        np3 = 3*innplc
        allocate(recvbuf(3,innplc))
        recvbuf = 0.0d0

        ifile1 = ifile + 1

        if(myid.eq.0) then
          do i = 1, innplc, samplePeriod
            write(ifile1,111)rays(1,i),rays(2,i),rays(3,i),rays(2,i)/gamma0
          enddo
          do i = 1, npr-1
            call MPI_RECV(recvbuf(1,1),np3,MPI_DOUBLE_PRECISION,&
                          i,1,commin,status,ierr)
            do j = 1, innplc, samplePeriod
              write(ifile1,111)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),recvbuf(2,j)/gamma0
            enddo
          enddo
          call flush(ifile1)
          close(ifile1)
        else
          call MPI_SEND(rays(1,1),np3,MPI_DOUBLE_PRECISION,0,1,&
                        commin,ierr)
        endif

111     format(4(1x,e21.12))

        deallocate(recvbuf)
        call MPI_BARRIER(commin,ierr)
        
        end subroutine output

      end module AccSimulatorclass
