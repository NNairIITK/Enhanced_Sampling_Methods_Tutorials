! 
! Fortran program to re-construct the 2D free energy profile for REST2 Simulation
!
! Authour: ANJI BABU KAPAKAYALA
!          C/O Prof. Nisanth Nair
!          Dept. of Chemistry
!	   IIT Kanpur, India.
!          anjibabu480@gmail.com
!            (  24/02/2020 )
!
!
! Definitions of variables used in this code
! 
! T          = Physical system temperature (300 K)
! kb          = 1.9872041*10^-3  (kcal K-1 mol-1)
! kbT        = kb*T
! md_steps    = The total MD steps 
! mtd_steps   = The total MetaD steps 
! gridmin1    = The minimum value of CV1
! gridmax1    = The maximum value of CV1
! griddiff1   = The witdh of bin along CV1
! gridmin2    = The minimum value of CV2
! gridmax2    = The maximum value of CV2
! griddiff2   = The witdh of bin along CV2
!
!
!  USAGE: COLVAR.0                    #Use ground replica to reconstruct FES
!         sed -i '/^#/d' COLVAR.0     #Remove all lines starts with '#'
!         gfortran FES_Reconstruction_REST2.f90 -o fes_rest2.x
!         ./fes_rest2.x
!
PROGRAM fes_rest2
IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddiff1, gridmin2, gridmax2, griddiff2
REAL*8 :: T, den, kbT
REAL*8 :: diff_s2, ds2, ss, hh, dum, num
REAL*8, ALLOCATABLE :: prob(:,:)
REAL*8, ALLOCATABLE :: cv1(:), cv2(:), vbias(:),s1,s2, s3, s4
INTEGER ::  md_steps, i, t_min, t_max
INTEGER :: i_s1, i_s2, i_s3, i_s4, nbin1, nbin2
INTEGER :: index1, index2,i_md
CHARACTER(20)::periodicity

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: kj_to_kcal = 0.239006

OPEN(1,FILE='input',STATUS='old')
OPEN(2,FILE='free_energy',STATUS='replace')
OPEN(11,FILE='COLVAR.0',STATUS='old')


CALL get_steps(11,md_steps)


read(1,*) T, t_min, t_max

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'

read(1,*) gridmin1, gridmax1, griddiff1
read(1,*) gridmin2, gridmax2, griddiff2
read(1,*) periodicity

! Keeping no of steps as maximum steps mentioned
md_steps = t_max

! Temperature in energy units
 kbT = kb*T

WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',T
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
WRITE(*,'(A,F9.2)')'Reweigtht: Grid min1    =',gridmin1
WRITE(*,'(A,F9.2)')'Reweigtht: Grid max1    =',gridmax1
WRITE(*,'(A,F9.2)')'Reweigtht: Grid width1  =',griddiff1
WRITE(*,*)
WRITE(*,'(A,F9.2)')'Reweigtht: Grid min2    =',gridmin2
WRITE(*,'(A,F9.2)')'Reweigtht: Grid max2    =',gridmax2
WRITE(*,'(A,F9.2)')'Reweigtht: Grid width2  =',griddiff2
WRITE(*,*)


ALLOCATE(cv1(md_steps),cv2(md_steps))

DO i_md=1,md_steps
 READ(11,*) dum,cv1(i_md),cv2(i_md)

   IF(periodicity .eq. "periodic") THEN
     PRINT*, "Given CVS are Periodic"
     IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
     IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
     IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
     IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0
   ELSE
   PRINT*, "Given CVS are Non-Periodic"
   ENDIF

END DO

nbin1 = NINT((gridmax1-gridmin1)/griddiff1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddiff2)+1

!WRITE(*,*) nbin1, nbin2
WRITE(*,'(A,I10)')'Reweigtht: nbin 1 ',nbin1
WRITE(*,'(A,I10)')'Reweigtht: nbin 2 ',nbin2

ALLOCATE(prob(nbin1,nbin2))

WRITE(*,*) 'Re-constructing Free energy'

den=0.d0
prob=0.d0

DO i_md=1,md_steps
   IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
      index1 = nint((cv1(i_md)-gridmin1)/griddiff1) +1
      index2 = nint((cv2(i_md)-gridmin2)/griddiff2) +1

         IF(index1.gt.0.and.index2.gt.0.and.index1.le.nbin1.and.index2.le.nbin2) then
            prob(index1,index2) = prob(index1,index2) + 1.d0
            den=den+1.d0
         END IF
  
   END IF
END DO

  dum=den*griddiff1*griddiff2
  den=1.d0/dum

DO i_s1=1,nbin1
s1=DFLOAT(i_s1-1)*griddiff1+gridmin1
   DO i_s2=1,nbin2
   s2=DFLOAT(i_s2-1)*griddiff2+gridmin2

           prob(i_s1,i_s2)=prob(i_s1,i_s2)*den
          
           WRITE(2,'(3E16.8)')s1, s2, -kbT*DLOG(prob(i_s1,i_s2))
   END DO
      WRITE(2,*)
END DO
 WRITE(*,'(A)')'Free energy has been written in file (free_energy)'

close(1)
close(2)
close(11)

DEALLOCATE(cv1, cv2)
DEALLOCATE(prob)

END PROGRAM fes_rest2
!-------------------------------------------------!
SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER :: iunit, nsteps
INTEGER :: ios
nsteps=0
REWIND(iunit)
  Read_Loop: DO
  READ(iunit,*,IOSTAT=ios)
    IF(ios.ne.0)EXIT Read_Loop
    nsteps=nsteps+1
  END DO Read_Loop
REWIND(iunit)
END SUBROUTINE
