!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! BeamBunchclass: Charged beam bunch class in Beam module of APPLICATION 
!                 layer.
! Version: 2.0
! Author: Ji Qiang, LBNL
! Description: This class defines the charged particle beam bunch 
!              information in the accelerator.
! Comments: 1) I have added the 3 attributes to the particle array:
!           x,px,y,py,t,pt,charge/mass,charge weight,id. We have moved
!           the charge*curr/freq/Ntot into the charge density calculation,
!           which is represented by the "charge weight" of a particle.
!           2) The map integrator does not work for multiple species, only
!              the Lorenze integrator works for the multiple species.
!----------------------------------------------------------------
      module BeamBunchclass
        use BPMclass
        use FFTclass
        use mpistub
        type BeamBunch
          !beam freq, current, part. mass and charge.
          double precision :: Current,Mass,Charge
          !# of total global macroparticles and local particles
          integer :: Npt,Nptlocal
          !particles type one.
          !Pts1(1,:) - z (m)
          !Pts1(2,:) - delta gamma 
          !Pts1(3,:) - longitudinal charge (C)
          double precision, pointer, dimension(:,:) :: Pts1
          !reference particle
          double precision, dimension(6) :: refptcl
        end type BeamBunch
      contains
        subroutine construct_BeamBunch(this,incurr,inkin,inmass,incharge,innp,&
                                       innplc,phasini)
        implicit none
        type (BeamBunch), intent(inout) :: this
        double precision, intent(in) :: incurr,inkin,inmass,&
                                        incharge,phasini
        integer, intent(in) :: innp,innplc
        integer :: myid, myidx, myidy,comm2d,commrow,commcol,ierr
        integer :: nptot,nprocrow,nproccol
   
        this%Current = incurr
        this%Mass = inmass
        this%Charge = incharge
        this%Npt = innp
        this%Nptlocal = innplc

        this%refptcl = 0.0d0
        this%refptcl(5) = phasini
        this%refptcl(6) = (inkin/this%Mass + 1.0)

        end subroutine construct_BeamBunch

        ! Drift beam half step using linear map for external field.
        ! this%Pts1(1,:) is z
        ! this%Pts1(2,:) is delta gamma
        subroutine map1(this,beamln,z,tau,bitype,nseg,nst,ihlf,flagfwd)
        implicit none
        type (BeamBunch), intent(inout) :: this
        type (BPM), intent(inout) :: beamln
        double precision, intent(inout) :: z
        double precision, intent (in) :: tau
        integer, intent(in) :: bitype,nseg,nst,flagfwd
        integer, intent(inout) :: ihlf
        integer :: i,j,k
        real*8 :: beta0,gam0,gambet0,tmppt,tmph,qmass
!for ideal RF cavity model
        real*8 :: vtmp,phi0lc,phi,gam,gambet,&
                  rk,pi,clite
        real*8, dimension(12) :: drange

        pi = 2*asin(1.0d0)
        clite = 2.99792458d8
        qmass = this%Charge/this%Mass

!        print*,"tau: ",tau

        if(bitype.eq.0) then
          gam0 = this%refptcl(6)
          gambet0 = sqrt(gam0**2-1.0d0)
          beta0 = gambet0/gam0
!          print*,"beta0: ",beta0
          do i = 1, this%Nptlocal
            tmppt = this%Pts1(2,i)+gam0
            tmph = sqrt(tmppt**2-1.0d0)
!            print*,"ipt00: ",i,this%Pts1(1,i),(1.d0/beta0-tmppt/tmph)
            this%Pts1(1,i) = this%Pts1(1,i)+(1.d0/beta0-tmppt/tmph)*tau
!            print*,"ipt: ",i,this%Pts1(1,i),(1.d0/beta0-tmppt/tmph)
          enddo
          z=z+tau
        else if(bitype.eq.103) then
          !ideal RF cavity model with entrance and exit focusing
            !max. accelerating gradient (V/m)
            vtmp = beamln%Param(3)/this%Mass
            !wave number
            rk = (2*pi*beamln%Param(4))/clite
            !synchronous phase
            phi0lc = beamln%Param(5)*pi/180

            !drift under constant acceleration

            if(flagfwd.eq.1) then
              gam0 = this%refptcl(6)
              gambet0 = sqrt(gam0**2-1.0d0)
              do i = 1, this%Nptlocal
                this%Pts1(1,i) = this%Pts1(1,i) + &
                               tau/gambet0**3*this%Pts1(2,i)
              enddo
              !update the reference particle energy
              this%refptcl(6) = this%refptcl(6) + tau*vtmp*cos(phi0lc)
            else
              !update the reference particle energy
              this%refptcl(6) = this%refptcl(6) + tau*vtmp*cos(phi0lc)
              gam0 = this%refptcl(6)
              gambet0 = sqrt(gam0**2-1.0d0)
              do i = 1, this%Nptlocal
                this%Pts1(1,i) = this%Pts1(1,i) + &
                               tau/gambet0**3*this%Pts1(2,i)
              enddo
            endif
!            print*,"vtmp2: ",vtmp,phi0lc,tau,this%refptcl(6)

            !if(nst.eq.nseg .and. mod(ihlf,2).eq.1) then
            !if(nst.eq.1 .and. mod(ihlf,2).eq.0) then
            if(mod(ihlf,2).eq.0) then
              do i = 1, this%Nptlocal
                phi = -rk*this%Pts1(1,i)+phi0lc
                !this%Pts1(2,i) = this%Pts1(2,i)+2*nseg*tau*vtmp*(cos(phi)-cos(phi0lc))
                this%Pts1(2,i) = this%Pts1(2,i)+2*tau*vtmp*(cos(phi)-cos(phi0lc))
              enddo
            endif
            z=z+tau
            ihlf = ihlf + 1
        else if(bitype.eq.4) then
            z=z+tau
        endif

        end subroutine map1

        !inputs are real units
        ! scatter field onto particles from grid.
        subroutine scatter(innp,innz,rays,fld,hz,zmin,tau,qmass)
        implicit none
        integer, intent(in) :: innp,innz
        real*8, intent(in) :: qmass
        double precision, intent (inout), dimension (3,innp) :: rays
        double precision, intent (in), dimension (innz) :: fld
        integer :: kx,kx1,i
        double precision :: ef
        double precision :: hz,hzi,zmin,tmpfld,tau

        hzi = 1.0d0/hz

        do i = 1, innp
          kx=(rays(1,i)-zmin)*hzi + 1 
          ef=(zmin+kx*hz-rays(1,i))*hzi
          kx1=kx+1

          tmpfld = fld(kx)*ef+fld(kx1)*(1.0d0-ef)
!          print*,"iez:",i,tmpfld,rays(1,i)
          rays(2,i) = rays(2,i) + tmpfld/qmass*tau
        enddo

!        print*,"zmm: ",minval(rays(1,:))
        
        end subroutine scatter

        ! rays(1,:) -> z
        ! rays(2,:) -> delta E
        ! rays(3,:) -> density weight (C/m)
        ! deposit particles onto grid.
        subroutine deposit(innp,innz,rays,rho,hz,zmin,commin)
        implicit none
        include "mpif.h"
        integer, intent(in) :: innp,innz
        integer, intent(in) :: commin
        double precision, intent (in), dimension (3,innp) :: rays
        double precision, intent (out), dimension (innz) :: rho
        double precision, dimension (innz) :: rhogl
        integer :: kx,kx1
        double precision :: ef
        double precision :: hz,hzi,zmin
        integer :: ierr,i

        hzi = 1.0d0/hz

        rho=0.0d0
        rhogl=0.0d0
        do i = 1, innp
          kx=(rays(1,i)-zmin)*hzi + 1 
          ef=((zmin-rays(1,i))+kx*hz)*hzi
          kx1=kx+1

          rho(kx) = rho(kx) + ef*rays(3,i)
          rho(kx1) = rho(kx1)+(1.0d0-ef)*rays(3,i)
        enddo

        call MPI_ALLREDUCE(rho,rhogl,innz,MPI_DOUBLE_PRECISION,&
             MPI_SUM,commin,ierr)
        rho = rhogl

!        print*,"sum rho1:",sum(rho)

        rho = rho/hz

        end subroutine deposit

       !longitudinal space-charge on a beam 
       !All inputs are in real units
       !recvdensz is rho(z) in unit of C/m
       !aa is the radius of the beam in m,
       !hz is in m
       !Output Ez also in real units V/m
       subroutine Lsc(Nz,recvdensz,ezwake,hz,aa,gam)
       implicit none
       integer, intent(in) :: Nz
       double precision, intent(in) :: hz,aa,gam
       double precision, dimension(Nz), intent(in) :: recvdensz
       double precision, dimension(Nz), intent(out) :: ezwake
       double precision, dimension(2*Nz) :: densz2n,densz2nout,&
                         greenwake,greenwakeout,tmpwk
       integer :: kz,twonz,one,ksign,kkzz,i,iz,iz1
       double precision :: scale,zz,zz1,dz
       real*8 :: Epsilon0,pi

       Epsilon0 = 8.854187817d-12
       pi = 2*asin(1.0d0)

!       print*,"aa: ",aa,gam,hz
!       print*,"sumaa: ",sum(recvdensz),hz,sum(recvdensz)*hz
  
       do kz = 1, Nz
          densz2n(kz) = recvdensz(kz) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       !longitudinal wakefield function
       do kz = 1, Nz+1
         zz = (kz-1)*hz
         greenwake(kz) = gam*zz/sqrt((gam*zz)**2+aa**2)
       enddo

       do kz = Nz+2, twonz
         greenwake(kz) = -greenwake(twonz-kz+2)
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       do kz = 1, 2
          greenwake(kz) = densz2nout(kz)*greenwakeout(kz)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1) = densz2nout(2*kz-1)*greenwakeout(2*kz-1)-&
                            densz2nout(2*kz)*greenwakeout(2*kz)
          greenwake(2*kz) = densz2nout(2*kz-1)*greenwakeout(2*kz)+&
                            densz2nout(2*kz)*greenwakeout(2*kz-1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

!-----------------------------------------
!direct brutal force
!
!       do kz = 1, Nz
!         tmpwk(kz) = 0.0d0
!         zz = (kz-1)*hz
!         do i = 1, Nz
!           zz1 = (i-1)*hz
!           dz = zz-zz1
!           tmpwk(kz) = tmpwk(kz)+recvdensz(i)*&
!                       gam*dz/sqrt((gam*dz)**2+aa**2)
!         enddo   
!       enddo
!-----------------------------------------

!       do kz = 1, Nz
!         write(1111,101)(kz-1)*hz,greenwakeout(kz),tmpwk(kz)
!       enddo
!101    format(3(1x,e16.8))

       do kz = 1, Nz
         ezwake(kz) = sum(recvdensz(1:kz))-sum(recvdensz(kz:Nz))-&
                      greenwakeout(kz)
       enddo
!       do kz = 1, Nz
!         write(1112,*)(kz-1)*hz,ezwake(kz)
!       enddo

       ezwake = ezwake*2/aa**2*hz/(4*pi*Epsilon0)
!       do kz = 1, Nz
!         write(1113,*)(kz-1)*hz,ezwake(kz)
!       enddo

!       print*,"ezwake:",sum(ezwake)

!------------------------------------------------------------------------------
       end subroutine Lsc

       !Inputs are real units.
       !Outputs are V/m.
       !longitudinal and transverse wakefield on a beam using the 
       !readin longitudinal and transverse wake functions
       subroutine wakefieldread(Nz,recvdensz,ezwake,&
                          hz,leng,ndatawk,wklong)
       implicit none
       integer, intent(in) :: Nz,ndatawk
       double precision, intent(in) :: hz, leng
       double precision, dimension(Nz), intent(in) :: recvdensz
       double precision, dimension(ndatawk), intent(in) :: wklong
       double precision, dimension(Nz), intent(out) :: ezwake
       double precision, dimension(2*Nz) :: densz2n,densz2nout,&
                         greenwake,greenwakeout
       integer :: kz,twonz,one,ksign,kkzz,i,iz,iz1
       double precision :: scale,zz,zziz
       real*8 :: hzwake,Scxlt

       Scxlt = 1.0d0

       hzwake = leng/(ndatawk-1)
  
       do kz = 1, Nz
          densz2n(kz) = recvdensz(kz) 
       enddo
       do kz = Nz + 1, 2*Nz
          densz2n(kz) = 0.0d0
       enddo

       twonz = 2*Nz
       one = 1
       ksign = 1
       scale = 1.d0
       call fftrclocal2_FFT(ksign,scale,densz2n,twonz,one,densz2nout)

       if(Nz*hz*Scxlt.gt.leng) then
         print*,"warning: total bunch length is greater than the readin wake field range!!"
       endif

       !longitudinal wakefield function
       greenwake(1) = 0.5d0*wklong(1)
       do kz = 2, Nz+1
         zz = (kz-1)*hz*Scxlt
         iz = zz/hzwake + 1
         iz1 = iz + 1
         if(iz1.gt.ndatawk) then
           iz = ndatawk - 1 
           iz1 = ndatawk  
         endif
         zziz = (iz-1)*hzwake
         greenwake(kz) = wklong(iz)+(wklong(iz1)-wklong(iz))*(zz-zziz)/hzwake
       enddo

       do kz = Nz+2, twonz
         greenwake(kz) = greenwake(twonz-kz+2)
       enddo
       do kz = 2, Nz
         greenwake(kz) = 0.0
       enddo

       call fftrclocal2_FFT(ksign,scale,greenwake,twonz,one,&
             greenwakeout)

       do kz = 1, 2
          greenwake(kz) = densz2nout(kz)*greenwakeout(kz)
       enddo
       do kz = 2, twonz/2
          greenwake(2*kz-1) = densz2nout(2*kz-1)*greenwakeout(2*kz-1)-&
                            densz2nout(2*kz)*greenwakeout(2*kz)
          greenwake(2*kz) = densz2nout(2*kz-1)*greenwakeout(2*kz)+&
                            densz2nout(2*kz)*greenwakeout(2*kz-1)
       enddo

       scale = 1.0d0/twonz
       ksign = -1
       call fftcrlocal2_FFT(ksign,scale,greenwake,twonz,one,&
            greenwakeout)

       !the "-" sign is from definition of longitudinal wake function, see Wangler's book
       do kz = 1, Nz
         ezwake(kz) = -greenwakeout(kz)*hz*Scxlt
       enddo
       
!------------------------------------------------------------------------------

       end subroutine wakefieldread

        !This subroutine calculates the 1d csr wakefield including
        !entrance, stead-state, and transitions effects.
        !Here, hx, rho,...are in real units. The return ezwake is in V/m.
        !The current version uses IGF corresponding to the four cases of Saldin et al.
        subroutine csrwakeTrIGF_FieldQuant(Nx,r0,ptmin,hx,blength,rhonew,gam,ezwake)
        implicit none
        integer, intent(in) :: Nx
        real*8 :: r0,ptmin,hx,blength,gam
        real*8, dimension(Nx) :: rhonew,ezwake
        real*8 :: xx,xxl,xx2,dx,tmprho,epstol,deltas,deltasmax,xxbar,psi,&
                  phim,pilc,yy,hx2,xconst,xk,psitmp,xslpN
        integer :: Ni,il,i,Nmax,j,islpN,islpN1
        integer :: myidlc,ierr,islp,islp1,islp0,jstart
        real*8 :: aa,uuh,uuh24,xslp,csrss,xxp,ssh,xk2,tcoef
        real*8 :: bb,cc,yh,phih,phpy2,csrtr,xblg,xxbarh,phimh,psimax
        real*8 :: psip2x2,psipx2,csrdrm1,csrdr1,phpxy2,phpy,csrss1,csrss2,&
                  csrtr1,csrtr2,csrdrm2,csrdr2
!        real*8 :: IcsrCaseA,IcsrCaseB,IcsrCaseC,IcsrCaseD

!        epstol = 1.0d-10 !tolerance for root finding
        epstol = 2.0d-9
        Nmax = 100 !maximum # of iteration in root finding
        pilc = 2*asin(1.0d0)
        phim = blength/r0
        xconst = 1./(4*pilc*8.854187817d-12)
        xk2 = gam/r0

        psitmp = 0.0
        ezwake(1) = 0.0 
        do i = 2, Nx
           xx = ptmin + (i-1)*hx 
           xxl = (xx/r0)**3*r0/24

           Ni = i
           ezwake(i) = 0.0

           xblg = (i-1)*hx

           !if(myidlc.eq.0) print*,"ii:",i,xxl,xx

           if(xx.le.0) then
           else if(xx.le.blength) then

             xslp = (xx/r0)*r0/2/gam/gam + xxl
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1

             !Case A 
             phih = xx/r0*gam
             ! IGF integral over the interior sample points for Case A:
             do j = 2, islp1-1
                 !write(*,*) 'Inside Case A!'
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrtr1-csrtr2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case A
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             !Add the lower end IGF integral for case A
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrtr1-csrtr2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr2 = IcsrCaseA(phih,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrtr1 = IcsrCaseA(phih,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrtr1-csrtr2)
             endif

             ! Case B (steady-state regime)
             jstart = max(islp1,1)
             ! IGF integral over the interior sample points for Case B:
             do j = jstart+2,Ni-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             enddo
             !add end integrals
             if(jstart.le.Ni-2) then
             !Add the upper end IGF integral for case B
                 xxp = ptmin + (Ni-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = ptmin + (Ni-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(Ni)*(csrss1-csrss2)
             !Add the lower end IGF integral for case B
                 if(islp1.lt.1) then
                   j = 1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)  
                   xxp = ptmin + (j-1)*hx 
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 else
                   j = islp1+1
                   xxp = ptmin + (j-1)*hx + hx/2
                   ssh = (xx-xxp)*gam**3/r0
                   csrss2 = IcsrCaseB(ssh,xk2)
                   xxp = xx-xslp
                   ssh = (xx-xxp)*gam**3/r0
                   csrss1 = IcsrCaseB(ssh,xk2)
                   ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
                 endif
             else if(jstart.eq.Ni-1) then
                 j = Ni
                 xxp = ptmin + (j-1)*hx 
                 ssh = (xx-xxp)*gam**3/r0
                 csrss2 = IcsrCaseB(ssh,xk2)
                 xxp = xx-xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrss1 = IcsrCaseB(ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrss1-csrss2)
             endif

           !  ezwake(i) = ezwake(i)*xconst

             !if(myidlc.eq.0) print*,"xxl: ",xxl,xx,ezwake(i),r0,ptmin,phim
             !if(myidlc.eq.0) print*,"xxl: ",xx,ezwake(i),xxl,xx-ptmin,xslp,hx

           else

             !if(myidlc.eq.0) print*,"CD:",i,islp,xslp,hx,xx,ezwake(i)

             xxbar = (xx - blength)
             xslp= (r0*phim+xxbar)/2/gam/gam + r0*phim**3/24*&
                   (r0*phim+4*xxbar)/(r0*phim+xxbar)
             islp0 = (xx-xslp-ptmin)/hx
             islp1 = islp0 + 1
             xslpN = xxbar/2.d0/gam/gam
             islpN = (xx-xslpN-ptmin)/hx
             islpN1 = islpN + 1

             ! Case C
             xxbarh = xxbar*gam/r0
             phimh = phim*gam
             ! IGF integral over the interior sample points for Case C:
             do j = 2, islp1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdrm1-csrdrm2)
             enddo
             if(islp1.ge.2) then
             !Add the upper end IGF integral for case C
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin + (islp1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             !Add the lower end IGF integral for case C
                 xxp = ptmin + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(1)*(csrdrm1-csrdrm2)
             ! Special case
             else if(islp1.eq.1) then
                 xxp = xx - xslp
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm2 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 xxp = ptmin 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdrm1 = IcsrCaseC(phimh,xxbarh,ssh,xk2)
                 ezwake(i) = ezwake(i) + rhonew(islp1)*(csrdrm1-csrdrm2)
             endif


             !Case D 
             psimax = phim*gam
             xxbarh = xxbar*gam/r0
             jstart = max(islp1,0)
             !write(13,*) 'xslp and xslpN',xslp,xslpN
             ! IGF integral over the interior sample points for Case D:
             do j = jstart+2,islpN1-1
                 xxp = ptmin + (j-1)*hx + hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (j-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(j)*(csrdr1-csrdr2)
             enddo
             if(islpN1.ge.(jstart+2)) then
             !Add the upper end IGF integral for case D
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 xxp = ptmin + (islpN1-1)*hx - hx/2
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(islpN1)*(csrdr1-csrdr2)
             !Add the lower end IGF integral for case D
                 xxp = ptmin + (jstart)*hx + hx/2 
                 ssh = (xx-xxp)*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             else if(islpN1.eq.(jstart+1)) then
                 ssh = xslpN*gam**3/r0
                 csrdr2 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ssh = xslp*gam**3/r0
                 csrdr1 = IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
                 ezwake(i) = ezwake(i) + rhonew(jstart+1)*(csrdr1-csrdr2)
             endif


           endif
           ezwake(i) = ezwake(i)*xconst
        enddo

        end subroutine csrwakeTrIGF_FieldQuant

           function IcsrCaseA(phih,ssh,xk2)
           implicit none
           double precision:: phih,ssh,xk2,IcsrCaseA
           double precision:: bb,cc,yh,phpy

           bb = 2*phih+phih**3/3-2*ssh
           cc = phih**2+phih**4/12-2*ssh*phih
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpy = (phih+yh)
           IcsrCaseA = xk2*(-(2*phpy+phih**3)/&
                    (phpy**2+phih**4/4)+1.0d0/ssh)
           end function IcsrCaseA


           function IcsrCaseB(ssh,xk2)
           implicit none
           double precision:: ssh,xk2,IcsrCaseB
           double precision:: aa,uuh
           aa = sqrt(64.0+144.0d0*ssh**2)
           uuh = (aa+12*ssh)**(1.0d0/3.0d0)-(aa-12*ssh)**(1.0d0/3.0d0)
           IcsrCaseB = xk2*(-4*uuh*(uuh**2+8)/ &
                    ((uuh**2+4.0d0)*(uuh**2+12.0d0)))
           end function IcsrCaseB


           function IcsrCaseC(phimh,xxbarh,ssh,xk2)
           implicit none
           double precision:: phimh,xxbarh,ssh,xk2,IcsrCaseC
           double precision:: bb,cc,yh,phpxya,phpxyb
           bb = 2*(phimh+xxbarh)-2*ssh+phimh**3/3+phimh**2*xxbarh
           cc = (phimh+xxbarh)**2+phimh**2*(phimh**2+4*phimh*xxbarh)/12-&
                      2*ssh*(phimh+xxbarh)
           yh = (-bb+sqrt(bb*bb-4*cc))/2
           phpxya = (phimh+xxbarh+yh)
           phpxyb = phimh*xxbarh+phimh*phimh/2.d0
           IcsrCaseC = xk2*(-2.d0*(phpxya+phimh*phpxyb)/ &
                          (phpxya**2+phpxyb**2)+1.0d0/ssh)
           end function IcsrCaseC

           function IcsrCaseD(xxbarh,psimax,ssh,xk2,epstol,Nmax)
           implicit none
           double precision:: xxbarh,ssh,xk2,epstol,IcsrCaseD
           double precision:: psi,psipxa,psipxb,psimax,x1,x2
           integer:: Nmax
           x1 = -epstol
           x2 = psimax*(1.d0+epstol)
           call root(x1,x2,epstol,Nmax,ssh,xxbarh,psi)
           psipxa = (psi+xxbarh)
           psipxb = psi*(xxbarh+psi/2.d0)
           IcsrCaseD = xk2*(-2.d0*(psipxa+psi*psipxb)/&
                    (psipxa**2+psipxb**2)+1.0d0/ssh)
           end function IcsrCaseD

     subroutine root(x1,x2,xacc,maxit,zeta,xxh,rtsafe)
!************************************************************************
!  This routine computes the root of the function that is evaluated in 
!  the subroutine 'funcd'. It is based on the subroutine 'root' of 
!  Numerical Recipes 9.4, which makes use of a Newton-Raphson method 
!  with root bracketing.  It has been modified to handle the two bracket 
!  endpoints carefully. The routine searches for a root in the interval
!  [x1,x2] with a tolerance given by 'xacc', and returns this value
!  as 'rtsafe'.  The maximum number of iterations allowed is 'maxit'.
!  C.E.M.
!***********************************************************************
     implicit none
     double precision:: rtsafe,x1,x2,xacc
     double precision:: xxh,zeta
     integer:: j,maxit
     double precision:: df,dx,dxold,f,fh,fl,temp,xh,xl
     call funcd(x1,xxh,zeta,fl,df)
     call funcd(x2,xxh,zeta,fh,df)
     if((fl>0.d0.and.fh>0.d0).or.(fl<0.d0.and.fh<0.d0)) then
           pause 'root must be bracketed in rtsafe'
           write(*,*) 'psimax,fl,fh = ',x2,fl,fh
     endif
     if(dabs(fl)< xacc) then
       rtsafe=x1
       return
     else if(dabs(fh)< xacc) then
       rtsafe=x2
       return
     else if(fl<0.d0) then
       xl=x1
       xh=x2
     else
       xh=x1
       xl=x2
     endif
     rtsafe=0.5d0*(x1+x2)
     dxold=dabs(x2-x1)
     dx=dxold
     call funcd(rtsafe,xxh,zeta,f,df)
     do j=1,maxit
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0.d0.or. &
&         dabs(2.d0*f)>dabs(dxold*df)) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl==rtsafe) return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp==rtsafe) return
        endif
        if(abs(dx)<xacc) return
        call funcd(rtsafe,xxh,zeta,f,df)
        if(f<0.d0) then
           xl=rtsafe
        else
           xh=rtsafe
        endif
     enddo
     pause 'root finding exceeding maximum iterations'
     return
     end subroutine

     subroutine funcd(psi,xxh,deltas,f,derivf)
!**********************************************************
!  This routine evaluates the function whose root produces
!  the retarded angle psi that is required for evaluating
!  the CSR kernel in Case D of Saldin et al.  The value
!  of the function is output as 'f', and its derivative
!  is output as 'derivf'.  C.E.M.
!*********************************************************
     implicit none
     double precision:: deltas,psi,term1,xxh,gamma,f,derivf
     double precision:: alpha,kappa,tau,theta
     f = psi**4/12+xxh*psi**3/3+psi**2+(2*xxh-2*deltas)*psi-&
               2*deltas*xxh+xxh*xxh
     derivf = psi**3/3+xxh*psi**2+2*psi+2*xxh-2*deltas
     end subroutine

      end module BeamBunchclass
