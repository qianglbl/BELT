!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! FFTclass: Fourier function class in Math Function module of FUNCTION layer.
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines the 3d FFT transformation subject to
!              open or periodic conditions, Fourier Sine transformation,
!              Complex-Complex, Complex-Real, and Real-Complex FFT.
! Comments:
!----------------------------------------------------------------
      module FFTclass
      contains
!----------------------------------------------------------------
      subroutine fftrclocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(in) :: x
      double precision, dimension(ny,nsizex), intent(out) :: y
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
         tempi = x(:,j)
         call realft(tempi,ny,ksign)
         y(:,j) = tempi*scale
         !y(:,j) = tempi
      end do

      end subroutine fftrclocal2_FFT

      subroutine fftcrlocal2_FFT(ksign,scale,x,ny,nsizex,y)
      implicit none
      integer, intent(in) :: ksign,ny,nsizex
      double precision, intent(in) :: scale
      double precision, dimension(ny,nsizex), intent(out) :: y
      double precision, dimension(ny,nsizex), intent(in) :: x
      real*8, dimension(ny) :: tempi
      integer :: i,j
      double precision :: t0

      ! Perform multiple FFTs with scaling:
      do j = 1, nsizex
           tempi = x(:,j)
           call realft(tempi,ny,ksign)
           y(:,j) = tempi*scale*2
           !y(:,j) = tempi
      end do

      end subroutine fftcrlocal2_FFT

      subroutine realft(data,n,isign)
      integer isign,n
      real*8 data(n)
      integer i,i1,i2,i3,i4,n2p3
      real*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=wr
        wis=wi
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      end subroutine realft

      subroutine four1(data,nn,isign)
      integer isign,nn
      real*8 data(2*nn)
      integer i,istep,j,m,mmax,n
      real*8 tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=wr*data(j)-wi*data(j+1)
            tempi=wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end subroutine four1

      subroutine sinft(y,n)
      integer n
      real*8 y(n)
      integer j
      real*8 sum,y1,y2
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0
      do 11 j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
11    continue
      call realft(y,n,+1)
      sum=0.0
      y(1)=0.5*y(1)
      y(2)=0.0
      do 12 j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
12    continue
      return
      END subroutine sinft
    end module FFTclass
