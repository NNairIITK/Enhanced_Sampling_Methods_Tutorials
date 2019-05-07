PROGRAM ste_dev
IMPLICIT NONE

INTEGER:: steps,ndim
REAL*8,ALLOCATABLE::mean(:),stdev(:)

ndim=2
OPEN(1,file="COLVAR")

ALLOCATE(mean(ndim))
ALLOCATE(stdev(ndim))

CALL get_steps(1,steps)

CALL cal_mean(steps,mean,1,ndim)

CALL standard_dev(1,steps,stdev,mean,ndim)

DEALLOCATE(mean,stdev)

END PROGRAM

SUBROUTINE get_steps(file_number,steps)
IMPLICIT NONE

INTEGER :: file_number, steps, ios

steps=0
DO
 READ(file_number,*,iostat=ios)
 IF(ios.NE.0)EXIT
  steps=steps+1
END DO
REWIND(file_number)
END SUBROUTINE get_steps

SUBROUTINE cal_mean(steps,mean,file_number,ndim)
IMPLICIT NONE

INTEGER::steps,file_number,i,ndim
REAL*8::mean(ndim)
REAL*8:: dummy,cv(ndim)

mean(1:ndim)=0.d0

DO i=1,steps
 READ(file_number,*) dummy,cv(1:ndim)
 mean(1:ndim)=mean(1:ndim)+cv(1:ndim)
END DO

REWIND(file_number)
mean(1:ndim)=mean(1:ndim)/steps
PRINT*,mean(1:ndim)

END SUBROUTINE

SUBROUTINE standard_dev(file_number,steps,stdev,mean,ndim)
IMPLICIT NONE

INTEGER:: i,file_number,steps,ndim
REAL*8:: stdev(ndim),mean(ndim)
REAL*8::cv(ndim),dummy

stdev(1:ndim)=0.d0
DO i=1,steps
 READ(file_number,*) dummy,cv(1:ndim)
 stdev(1:ndim)=stdev(1:ndim)+(cv(1:ndim)-mean(1:ndim))**2
END DO
stdev(1:ndim)=dsqrt(stdev(1:ndim)/steps)
PRINT*,stdev(1:ndim)
END SUBROUTINE