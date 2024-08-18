!----------------------------------------------------------------
! (c) Copyright, 2023 by the Regents of the University of California.
! Inputclass: Input class 
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: This class defines functions to input the global
!              beam and computational parameters and the lattice input
!              parameters in the accelerator.
! Comments: J. Q. modified the source code so that the user can put
!           coments line starting with "!" for each number line in the
!           input file "eblt.in".
!----------------------------------------------------------------
      module Inputclass
        interface in_Input
          module procedure in1_Input, in2_Input
        end interface
      contains
        subroutine init_Input(time)
        implicit none
        real*8 :: time
        integer :: ierr

        end subroutine init_Input
 
        ! Input all parameters except beam line element parameters.
        subroutine in1_Input(onp,onz,distparam,nparam,obcurr,obkenergy,obmass,&
        obcharge,obfreq,oxrad,oyrad,onblem,oflagfwd,oflagdist)

        implicit none
        integer, intent(out) :: onp,onz
        integer, intent(out) :: onblem,oflagfwd,oflagdist
        integer, intent(in) :: nparam
        double precision, dimension(nparam), intent(out) :: distparam
        double precision, intent(out) :: obcurr,obkenergy,obmass
        double precision, intent(out) :: obcharge,obfreq,&
                                         oxrad,oyrad
        double precision :: xjunk
        integer :: my_rank,nproc,ierr,np,itot,njunk1,njunk2,njunk3
        character*1 comst
        integer :: ii,jj,i

          print*,"Read input data from file - eblt.in:"
          open(unit=13,file='eblt.in',status='old')

          ii = 0
          jj = 0 
30        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 30
          else
            backspace(13,err=789)
            read(13,*)onp,onz,oxrad,oyrad,oflagfwd,oflagdist
            ii = ii+1
          endif 
          distparam = 0.0
80        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 80
          else
            backspace(13,err=789)
            read(13,*)distparam(1:10)
            ii = ii+1
          endif 
90        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 90
          else
            backspace(13,err=789)
            read(13,*)distparam(11:20)
            ii = ii+1
          endif
102       continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 102
          else
            backspace(13,err=789)
            read(13,*)obcurr,obkenergy,obmass,obcharge,obfreq
            ii = ii+1
          endif

          !count the # of beam line elements.
          itot=0
          njunk3 = 0
123       continue
            read(13,*,end=789)comst
            if(comst.ne."!") then
              backspace(13,err=789)
              read(13,*,end=789)xjunk,njunk1,njunk2,njunk3
              itot = itot + 1
            endif
            if(njunk3.eq.-99)then
              goto 789
            endif
          goto 123
  789     continue
          onblem=itot
!          write(6,*)'onblem = ',onblem
          rewind(13)

          do i = 1, jj
            read(13,*)comst
          enddo

        end subroutine in1_Input

        ! Input beam line element parameters.
        subroutine in2_Input(onblem,oblength,obnseg,obmpstp,&
        obtype,value1,value2,value3,value4,value5,value6,value7,value8)
        implicit none
        integer,intent(in) :: onblem
        integer,intent(out) :: obnseg(onblem)
        integer,intent(out) :: obmpstp(onblem)
        integer,intent(out) :: obtype(onblem)
        double precision,intent(out) :: oblength(onblem)
        double precision,dimension(onblem),intent(out) :: value1,value2,&
        value3,value4,value5,value6,value7,value8
        integer :: i,irf
        integer :: myrank,ierr
        character*1 comst

        value1 = 0.0
        value2 = 0.0
        value3 = 0.0
        value4 = 0.0
        value5 = 0.0
        value6 = 0.0
        value7 = 0.0
        value8 = 0.0

          i=0
123       continue
            read(13,*,end=789)comst
            if(comst.ne."!") then
              backspace(13,err=789)
              i = i + 1
              read(13,*)oblength(i),obnseg(i),obmpstp(i),obtype(i),&
              value1(i),value2(i),value3(i),value4(i),value5(i),value6(i),&
              value7(i),value8(i)
            endif
            if(obtype(i).eq.-99)then
              goto 789
            endif
          goto 123
789       continue

          do i = 1, onblem
!          print*,"value11: ",value11(i),value12(i)
          end do
          close(13)


        end subroutine in2_Input

      end module Inputclass

