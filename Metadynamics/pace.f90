PROGRAM check_pace
IMPLICIT NONE

REAL*8, ALLOCATABLE :: cv(:,:)
REAL*8, ALLOCATABLE :: step(:)
REAL*8 :: dum3, init_pos1, diff1, width

INTEGER :: md_steps, pace,number
INTEGER :: i,j,k,m,avg_pace,steps,sum
REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)
INTEGER:: column
OPEN(1,file='input',status='old')
READ(1,*) width
!READ(1,*) column

column=1
OPEN(4,file='COLVAR',status='old')
OPEN(10,file='PACE.dat',status='replace')

CALL step_count(4,md_steps)

!PRINT *, 'md steps from count =', md_steps
REWIND(4)

md_steps=md_steps-5

ALLOCATE(cv(2,md_steps))

READ(4,*)
READ(4,*)
READ(4,*)
READ(4,*)
READ(4,*)

DO k = 1, md_steps
     READ(4,*) dum3, cv(1,k), cv(2,k), dum3
END DO

init_pos1=cv(column,1)
diff1=0.d0
width=width*1.5d0
j=0

DO i=1,md_steps
  diff1= cv(column,i) - init_pos1

  IF((diff1).GT.pi) diff1=diff1 -2.d0*pi
  IF((diff1).LT.-pi) diff1=diff1 +2.d0*pi

   j=j+1
   IF(ABS(diff1).GT.width)THEN
     pace=j
     j=0
     init_pos1=cv(column,i)
     WRITE(10,*)  pace
   END IF
END DO
REWIND(10)

CALL step_count(10,steps)
REWIND(10)

sum=0

DO m=1,steps
READ(10,*) number
sum = sum+number 
END DO

avg_pace = sum/steps

PRINT *, 'AVERAGE PACE =', avg_pace

end program

subroutine step_count(file_number,steps)
integer :: file_number, steps, ios,i
steps=0
do
 read(file_number,*,iostat=ios)
 if(ios.ne.0) exit
  steps=steps+1
end do
end subroutine
