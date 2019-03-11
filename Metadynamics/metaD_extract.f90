      program extract
      implicit none
      integer:: nmtd,i,j,nbin,bla
      real*8, allocatable :: cv(:)
      real*8 :: gridmin,gridmax,griddif,s1,height,width,pot,dum,bla1
      nmtd=1000
      allocate (cv(nmtd))
      height=3.0             ! gaussian hight in KJ/Mol
      width=0.1              ! width in rad
      gridmin=0.0
      gridmax=7.4
      griddif=0.01

      nbin=nint((gridmax-gridmin)/griddif)+1

      do i=1,nmtd
      read(*,*)bla1,cv(i)
      end do

      do j=1,nbin
      pot=0.0
      s1=(dfloat(j-1))*griddif+gridmin
      do i=1,nmtd
      dum=s1-cv(i)
      pot=pot-height*exp(-0.5*(dum/width)**2)    !Sum of Gaussians = F.E
      end do
      write(*,*)s1,pot
      enddo


      end program
