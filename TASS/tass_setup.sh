#  Bash Script to perform Temperature Accelerated Sliced Sampling (TASS) Simulation on Alanine dipeptide in vacuum using gromacs
#
#            Authour:  Anji Babu Kapakayala
#                      C/O Prof. Nisanth N. Nair
#                      Dept. of Chemistry
#                      IIT Kanpur, India.
#                      
#             USAGE:  sh tass_setup.sh
#             
#             Requirements:  Plumed and Gromacs should be installed with mpi
#                            gfortran, gnuplot, mpif90
#
#	      INPUTS:  ala_di-pep.gro and ala_di-pep.top
#                          
#             
#             Description:  This script will generate required input files and runs the umbrella sampling simulation
#                           for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns to obtain the
#			    initial structures for TASS simulations. And performs TASS then constructs the 2D free energy 
#         		    along Phi and Psi collective variables
#             
#             On the fly, It will create following directories and stores related files in those directories:
#             
#             MIN : write necessary input files and does the energy minimization.
#             NVT : Writes the required input files and performs equilibration here.
#	      US  : Performs the Umbrella Sampling Simulations here in respective directory for every CV value
#	      TASS:  Runs TASS Simulations in respective sub directories
#	      TASS/INPUTS: Stores all required input files here
#	      TASS/TASS_*: Respective directory for umbrella window        
#	      TASS/ANALYSIS: Writes all required inputs here to do analysis 
#	      TASS/PROB: Stores all the unbiased probabilities after reweighting metadynamics
#	      TASS/WHAM: Performs WHAM run here 
#
#
#
# 	     Cheers !!
#	    Anji Babu
#
#!/bin/bash
function write_mdp() {
case "$1" in
  min)
cat >> min.mdp<<EOF
; For Vacuum condition
; minim.mdp - used as input into grompp to generate em.tpr
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps		= 50000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 10	    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = group
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
coulombtype	    = cutoff
rcoulomb	    = 2.0		; Short-range electrostatic cut-off
rvdw		    = 2.0		; Short-range Van der Waals cut-off
pbc		    = no 		; Periodic Boundary Conditions (yes/no)
rlist               = 2.0
EOF
;;
 nvt)
cat >> nvt.mdp <<EOF
title                    = Alanine dipeptide in vacuum
;Preprocessor
cpp                      = /lib/cpp
;Directories to include in the topology format
;include                 = -I../top
;Run control: A leap-frog algorithm for integrating Newton's equations. 
integrator               = md
;Total simulation time: 100 ps
;time step in femtoseconds 
dt                       = 0.002
;number of steps
nsteps                   = 1000000
;frequency to write coordinates to output trajectory file
nstxout                  = 100
nstvout                  = 100
nstfout                  = 100
;frequency to write energies to log file
nstlog                   = 100
;frequency to write energies to energy file
nstenergy                = 100
;frequency to write coordinates to xtc trajectory 
nstxtcout                = 100
;group(s) to write to xtc trajectory
xtc_grps                 = System
;group(s) to write to energy file 
energygrps               = Protein
;Frequency to update the neighbor list (and the long-range forces, 
;when using twin-range cut-off's). 
nstlist                  = 10
;Make a grid in the box and only check atoms in neighboring grid cells 
;when constructing a new neighbor list every nstlist steps. 
ns_type                  = grid
;cut-off distance for the short-range neighbor list
;treatment of electrostatic interactions
coulombtype = cutoff
;treatment of van der waals interactions
rvdw = 2.0
rlist = 2.0
rcoulomb = 2.0 
cutoff-scheme=group
comm-mode=Angular
; Periodic boudary conditions in all the directions 
pbc                      = no
;Temperature coupling
tcoupl                   = v-rescale
tc-grps                  = Protein
tau_t                    = 0.1
ref_t                    = 300
;Velocity generation
gen_vel                  = yes 
gen_temp                 = 300
gen_seed                 = 173529
;Constrain all bonds
constraints              = all-bonds
EOF
;;
*) echo " Invalid argument with write_mdp ";;
esac
}
#--------------------------------------------------------#
# Function to write plumed input file
function write_plumed_input() {
case "$1" in
us)
cat >>plumed.dat<<EOF
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

# Umbrella Bias
us: RESTRAINT ARG=phi KAPPA=500 AT=$2
# monitor the two variables
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR
EOF
;;
#
tass)
cat >>plumed.dat<<EOF
#-------------------------------------------------------#
# PLUMED INPUT TO PERFORM TASS
#-------------------------------------------------------#
#set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
#-------------------------------------------------------#
# EXTENDED LAGRANGIAN ACTIVATED
ex: EXTENDED_LAGRANGIAN ARG=phi,psi KAPPA=5260,5260 TAU=0.062,0.062 FRICTION=100,100 TEMP=900
#-------------------------------------------------------#
# UMBRELLA SAMPLING ACTIVATED along Phi
restraint-ex.phi_fict: RESTRAINT ARG=ex.phi_fict KAPPA=500 AT=$2
#-------------------------------------------------------#
# METAD ACTIVATED along Psi
# with height equal to 2.4 kJoule/mol,
# Well-tempered metadynamics is activated,
metad: METAD ARG=ex.psi_fict PACE=500 HEIGHT=2.4 SIGMA=0.05 FILE=HILLS BIASFACTOR=2.0 TEMP=900.0
#-------------------------------------------------------#
# monitor the auxiliary variables
PRINT STRIDE=10 ARG=ex.phi_fict,ex.psi_fict FILE=COLVAR
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR_REAL
PRINT ARG=ex.phi_vfict,ex.psi_vfict FILE=TEMP
#------------------------------------------------------#
EOF
;;
*)echo " Invalid plumed option";;
esac
}
#================================================================#
#Prepare Initial structure using gromacs & AMBERFF14SB forcefields
#gmx pdb2gmx -f ala2.pdb -ter -ignh
function Prepare_initial_structure() {
gmx_mpi pdb2gmx -f ala_di-peptide.pdb -ter -ignh <<EOF
5 
7
EOF
mkdir INITIAL_STRUCTURES
mv *.gro *.top *.itp INITIAL_STRUCTURES/
}
#--------------------------------------------------------------#
# Energy Minimisation in vacuum
#gmx grompp -f min.mdp -c conf.gro -p topol.top -o min.tpr
function Perform_Minimisation() {
mkdir MIN
cd MIN
cp ../*.gro .
cp ../*.top .
#cp ../INITIAL_STRUCTURES/* .
write_mdp min
gmx_mpi grompp -f min.mdp -c *.gro -p *.top -o min.tpr -maxwarn 10
gmx_mpi mdrun -deffnm min -v
cd ..
}
#--------------------------------------------------------------#
# Equilibration NVT
#gmx grompp -f nvt.mdp -c min.gro -p topol.top -o nvt.tpr
function Perform_Equilibration() {
mkdir NVT
cd NVT
#cp ../INITIAL_STRUCTURES/* .
cp ../MIN/min.gro .
cp ../MIN/*.top .
write_mdp nvt
gmx_mpi grompp -f nvt.mdp -c min.gro -p *.top -o nvt.tpr -maxwarn 10
gmx_mpi mdrun -deffnm nvt -v 
cd ..
}
#--------------------------------------------------------------#
# Perform TASS simulations
function Perform_TASS() {
mkdir TASS
cd TASS
mkdir INPUTS
sed -i "s/nsteps                   = 1000000/nsteps                   = 500000/g" nvt.mdp
write_mdp nvt 
mv nvt.mdp INPUTS
#cp ../INITIAL_STRUCTURES/topol.top INPUTS
#cp ../INITIAL_STRUCTURES/*.itp INPUTS
#cp ../NVT/nvt.gro INPUTS
cp ../NVT/*.top INPUTS
#------------------------------------------#
for i in `seq -3.2 0.2 3.2`;do
mkdir TASS_$i
cd TASS_$i
cp ../INPUTS/* .
write_plumed_input tass $i
cp ../../US/US_$i/us.gro .
gmx_mpi grompp -f nvt.mdp -c us.gro -p *.top -o tass.tpr -maxwarn 10
gmx_mpi mdrun -v -deffnm tass -plumed plumed.dat
cd ../
done
cd ..
}
#--------------------------------------------------------------#
# Perform US simulations
function Perform_US() {
mkdir US
cd US
mkdir INPUTS
write_mdp nvt 
sed -i "s/nsteps                   = 1000000/nsteps                   = 500000/g" nvt.mdp
mv nvt.mdp INPUTS
#cp ../INITIAL_STRUCTURES/topol.top INPUTS
#cp ../INITIAL_STRUCTURES/*.itp INPUTS
cp ../NVT/nvt.gro INPUTS
cp ../NVT/*.top INPUTS
#------------------------------------------#
firstTime="yes"
# Submit US Runs
#for i in `seq 1.4 0.2 3.2` `seq 1.2 -0.2 -3.2`;do
for i in `seq -3.2 0.2 3.2`;do

mkdir US_$i
cd US_$i
cp ../INPUTS/* .
write_plumed_input us $i
#if [ "$firstTime"  == "yes" ];then
gmx_mpi grompp -f nvt.mdp -c nvt.gro -p *.top -o us.tpr -maxwarn 10
#firstTime="no"
#else
#if [ "$i" -gt "0.0" ];then
#j=$( bc <<<  "$i-0.2")
#rm nvt.gro
#echo "$j" >> LOG
#else
#j=`echo "$i+0.2"|bc`
#fi
#cp ../US_$j/us.gro .
#gmx_mpi grompp -f nvt.mdp -c us.gro -p *.top -o us.tpr -maxwarn 10
#fi
gmx_mpi mdrun -v -deffnm us -plumed plumed.dat
cd ..


done
cd ..
}
#-----------------------------------------------------------------#
function gnuplot_2d_script() {

cat >>plot_fes.gnp<<EOF
reset
set term png #output terminal and file
set output "FES_2D.png"
set border lw 2
set pm3d interpolate 0,0
set cntrp bspline
set cntrp linear
set cntrp points 10
set cntrp order 10
set contour both
set cntrp levels 100
set xrange [-3.14:3.08]
set yrange [-3.14:3.09]
set xlabel "Phi"
set ylabel "Psi"
set key right outside
set view 0,0
set palette color
set palette model HSV
set palette rgbformulae -3,-2,2
#set palette defined ( -30 "green", -24 "blue", -12 "red", -0.5 "orange", 0 "white" )
set colorbox vertical user origin 0.08,0.26 size 0.02,0.51
#set view equal xy
unset ztics
unset clabel
unset key
#
#set zrange [-50:54.1]
set cntrp lev incr -30,2.05,200
splot 'free_energy' us 1:2:3 with pm3d
EOF
#
}
#-------------------------------------------------------------------#

function write_gnuplot_script() {
cat >>plot.gnp<<EOF
reset
set term png #output terminal and file
set output "All_histograms.png"
set border lw 2
#set tics font' ,15'
set ylabel " Probability Distribution" 
set xlabel " Potential Energy" 
set title "P.E Overlap REMD" font' ,15'
#Plot
plot "pe_hist_0.dat" u 1:2 w l lw 2 title'Replica_0', "pe_hist_1.dat" u 1:2 w l lw 2 title'Replica_1',"pe_hist_2.dat" u 1:2 w l lw 2 title'Replica_2',"pe_hist_3.dat" u 1:2 w l lw 2 title'Replica_3'                                                                                                                                                                                                     
EOF
}
#===============================================#
# Check for gromacs & plumed before starting simulation

function check_plumed() {
plumed --is-installed
if [ $? -eq 0 ];then
echo "Checking plumed: Found"
gmx_mpi mdrun -h &> temp
grep " -plumed" temp &>/dev/null
  if [ $? -eq 0 ];then
    echo "Checking plumed: Plumed patched with Gromacs"
  else
    echo "Checking plumed: Plumed did NOT patched with Gromacs"
    echo "Exiting..!"
    exit
  fi
rm temp
else
echo "Checking plumed: Plumed NOT found"
echo "Thank you !"
exit
fi
}
#-----------------------------------------------#
function check_gromacs() {
gmx_mpi mdrun -h &> /dev/null
if [ $? -eq 0 ];then
echo "Checking Gromacs: Found"
else
echo "Checking Gromacs: Gromacs mpi NOT found"
echo "Thank you !"
exit
fi
}
#----------------------------------------------#
# function to check for initial input ala_di-peptide.pdb
function check_for_pdb() {
if [  -f ./ala_di-peptide.pdb ]; then
 echo "Checking PDB: OK"
else    
echo "Checking PDB: pdb file not found! Sorry."
exit
fi
}
#----------------------------------------------#
# function to check for initial input gro & top 
function check_for_gro_top() {
if [  -f ./*.gro ]; then
 echo "Checking GRO: OK"
if [  -f ./*.top ]; then
 echo "Checking TOP: OK"
else
 echo "Checking TOP: NOT FOUND"
exit
fi
else    
 echo "Checking GRO: NOT FOUND"
exit
fi
}
#------------------------------------------------#
# Function to write WHAM Codes
function Write_WHAM_Code() {
cat >>wham.f90 <<EOF
!***********************************************!
!THE PROGRAM BEGINS.
!***********************************************!
program fes_calc
implicit none
integer i, iter, ncv, umbr_n, nbin1, nbin2
real*8, allocatable :: grid0(:,:), v(:), prob(:,:),biased_prob(:,:,:),grid(:,:)
real*8 :: kt, toler, dummy
logical :: parent
integer :: i_umbr, i_s1, i_s2
character (len=50) cvfile,outputfile,f1,f3
real*8 cnvg
real*8, allocatable :: umbr_mean(:), umbr_k(:),a(:)
integer, allocatable :: nmax(:)
real*8 :: dummy_1,dummy_2,dummy_3
integer :: rank, gleng1_min,gleng1_max,gleng2,ngrid
real*8, parameter :: kb=1.9872041E-3 !kcal K-1 mol-1
real*8, parameter :: au_to_kcal = 627.51
real*8, parameter :: kj_to_kcal = 0.239006

CALL MPI_Start
CALL Set_Parent(parent)


if(parent)then
  open (unit=1, file='input',status='old')
  !read number of cv, number of umbrella,kt in energy unit
  read(1,*)ncv, umbr_n, kt
end if
  kt=kb*kt

!broadcast ncv
CALL IBcast(ncv,1)
CALL IBcast(umbr_n,1)
CALL RBcast(kt,1)

allocate(grid(3,ncv))
allocate(grid0(3,ncv))
allocate(a(umbr_n))
allocate(umbr_mean(umbr_n))
allocate(umbr_k(umbr_n))
allocate(nmax(umbr_n))
a=1.d0

if (parent) then
 !read grid_min, grid_max, grid_bin
 do i=1,ncv
    read(1,*)grid0(1:3,i)
    write(*,'(i10,3f16.6)')i,grid0(1:3,i)
 end do
end if

if (parent) then
if (ncv.eq.2) then
nbin1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
nbin2=nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1
else if (ncv.ge.3)then
STOP '3 or more CVs not implemented'
end if
end if

!broadcast grids and bin info
CALL RBcast(grid0,3*ncv)
CALL IBcast(nbin1,1)
CALL IBcast(nbin2,1)

allocate(biased_prob(nbin1,nbin2,umbr_n))
allocate(prob(nbin1,nbin2))

if (parent) then
   open( unit =2, file= 'whaminput', status = 'old')
   read(2,*)toler
   do i_umbr=1, umbr_n
      write(*,*) 'umbrella simulation #', i_umbr
      !reads probability file name
      read(2,'(a)') cvfile
      write(*,*)'CVFILE',cvfile
      !reads force constant, r0, max points in that umbrella
      read(2,*) umbr_mean(i_umbr),umbr_k(i_umbr),nmax(i_umbr)
        umbr_k(i_umbr)=umbr_k(i_umbr)*kj_to_kcal !converted in kcal
       if(umbr_mean(i_umbr).gt. 3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)- 6.28d0
       if(umbr_mean(i_umbr).lt.-3.14d0)umbr_mean(i_umbr)=umbr_mean(i_umbr)+ 6.28d0

      f1 = cvfile
      open( unit=3, file=f1, status='old' )
      do i_s1=1,nbin1 !US
         do i_s2=1,nbin2 !MTD
         read(3,*)dummy,dummy,biased_prob(i_s1,i_s2,i_umbr)
         end do
      end do
   enddo
end if

call RBcast(toler,1)
call RBcast(umbr_mean,umbr_n)
call RBcast(umbr_k,umbr_n)
call IBcast(nmax,umbr_n)
call RBcast(biased_prob,nbin1*nbin2*umbr_n)

!PERFORMS WHAM.
 if (parent)  write(*,*) 'wham begins'
 CALL DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)

 gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1

 write(*,*)'new_grid', gleng1_min, gleng1_max,gleng2, rank, grid(1,1)
  iter=0
  scf_loop : do
       iter=iter+1
       call wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
       if (parent) write(*,*) 'iteration #', iter, 'convergence =',cnvg
       if (mod(iter,100) .eq. 0 ) then
       call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
       endif
       if((cnvg.lt.toler).or.(iter.ge.20000))then
       if (parent) write(*,*)'** convergence achieved **'
       exit scf_loop
       end if
   end do scf_loop

    !prints the free energy.
    call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
call MPI_Stop
end program

!***********************************************!


 subroutine print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
!prints pmf 
 implicit none
 integer:: i_s1, i_s2
 integer:: nbin1, nbin2, ncv
 real*8 :: prob(nbin1,nbin2), grid0(3,ncv)
 real*8 :: s1, s2, kt, dum
 logical:: parent
 real*8, allocatable :: prob1(:,:)
 character (len=50):: f2

 allocate(prob1(nbin1,nbin2))
 prob1=0.0
 CALL GlobSumR(prob,prob1,nbin1*nbin2)
 if (parent)then
 f2= 'free_energy'
 open( unit =7 , file = f2, status =  'unknown' )
 do i_s1=1,nbin1 !US cv
  s1=DFLOAT(i_s1-1)*grid0(3,1)+grid0(1,1)
  do i_s2=1,nbin2 !MTD cv
     s2=DFLOAT(i_s2-1)*grid0(3,2)+grid0(1,2)
     dum= -kt*DLOG(prob1(i_s1,i_s2))
     write(7,'(4E16.8)') s1, s2, dum, prob(i_s1,i_s2)
  enddo
  write(7,*)
 enddo
 write(*,*) 'free energy written in ',f2
 close(7)
 end if
 endsubroutine
!**************************************************************************************************!

subroutine wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
!performs wham scf.
implicit none
integer:: i_s1, i_s2, i_umbr, nbin1, nbin2, umbr_n,ncv
integer:: nmax(umbr_n)
real*8 :: umbr_k(umbr_n), umbr_mean(umbr_n),dum
real*8 :: num, den, dummy_v, dummy_s1, avg, del_s1, kt, dummy, cnvg
real*8 :: prob(nbin1,nbin2),a(umbr_n)
real*8 :: grid0(3,ncv), biased_prob(nbin1,nbin2,umbr_n),grid(3,ncv)
integer:: rank, gleng1_min,gleng1_max,gleng2,ngrid
real*8,allocatable :: dummy_a1(:),dummy_a(:)

allocate(dummy_a(umbr_n))
allocate (dummy_a1(umbr_n) )

dummy_a = 0.0d0

!calculates probability at each grid_point.
do i_s1 =gleng1_min,gleng1_max !over US cv
do i_s2 =1,gleng2 !over MTD cv
   num = 0.0d0
   den = 0.0d0
   dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)

   !calculates probability.
   do i_umbr=1,umbr_n
      del_s1=dummy_s1-umbr_mean(i_umbr)
      if ( del_s1 .gt. 3.14d0 ) del_s1 = del_s1 - 6.28d0
      if ( del_s1 .lt.-3.14d0 ) del_s1 = del_s1 + 6.28d0
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      num=num+dfloat(nmax(i_umbr))*biased_prob(i_s1,i_s2,i_umbr)
      den=den+dfloat(nmax(i_umbr))*a(i_umbr)*dummy_v
   enddo

   prob(i_s1,i_s2)=num/den
   if(prob(i_s1,i_s2).ne.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove NaN
   if(prob(i_s1,i_s2)+1.eq.prob(i_s1,i_s2)) prob(i_s1,i_s2)=1.0D-16 !remove infinity

!calculate a.
      dum=grid0(3,1)*grid0(3,2)
   do i_umbr=1,umbr_n
      del_s1=dummy_s1 - umbr_mean(i_umbr)
      if ( del_s1 .ge. 3.14d0 ) del_s1 = 6.28d0 - del_s1
      if ( del_s1 .le.-3.14d0 ) del_s1 = 6.28d0 + del_s1
      dummy_v=dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      dummy_a(i_umbr)=dummy_a(i_umbr) +dum*dummy_v*prob(i_s1,i_s2)
   enddo
enddo !end of MTD cv loop
enddo !end of US cv loop

!=======================================================================================
CALL GlobSumR(dummy_a,dummy_a1,umbr_n)
do i_umbr=1,umbr_n
dummy_a(i_umbr)=dummy_a1(i_umbr)
end do
!========================================================================================

!finds convergence and update a.
 avg = 0.0d0
 cnvg = 0.0d0
 do i_umbr=1,umbr_n
 dummy_a(i_umbr) = 1.0d0/dummy_a(i_umbr)
 cnvg = cnvg + dabs(dlog(dummy_a(i_umbr))-dlog(a(i_umbr)))
 a(i_umbr) = dummy_a(i_umbr)
 enddo
 cnvg = kt*cnvg
 end subroutine
!**************************************************************************************************!

SUBROUTINE DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)
!Distribute X grid over processors by mapping grid0 to grid 
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
!
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m,ngrids_y,ngrids_z,ngrids_o
INTEGER :: gleng1_min, gleng1_max

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

write(*,*)'NCPUS',ncpu
DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO
rank=0

IF(ncv.EQ.1)THEN
 ngrids_y=1
 ngrids_z=1
ELSE IF(ncv.EQ.2)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=1
ELSE IF(ncv.EQ.3)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=(nint((grid0(2,3)-grid0(1,3))/grid0(3,3))+1)
END IF

if(ncpu.eq.1) then
gleng1_min=1
gleng1_max=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

!Distribute X grids
if(icpu.eq.0)WRITE(*,'(3A12,3A16)') 'CPU','CV', 'GRID SIZE', 'GRID MIN', 'GRID MAX', 'GRID BIN'
CALL Sync_procs
IF(ncpu.GT.1)THEN
  ngrids=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  write(*,*)'NEW_GRID',ngrids,icpu,ncpu
  ngrids_o=ngrids
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids_o,ncpu)
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids_m-1)*grid0(3,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids-1)*grid0(3,1)
  END IF
  CALL Sync_procs
  WRITE(*,'(3I12,3F16.6)') icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  rank=ngrids_z*ngrids_y*ngrids*icpu
  gleng1_min=ngrids*icpu+1
  gleng1_max=ngrids*(icpu+1)
  if(icpu.eq.ncpu-1) gleng1_max=ngrids*(icpu+1)+mod(ngrids_o,ncpu)
END IF
END
!****************************************************************************************!
!****************************************************************************************!

SUBROUTINE MPI_Start()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
  call MPI_INIT(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE MPI_Stop()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
call MPI_FINALIZE(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: ncpu, i_err
!ncpu=1
!#if defined (_PARALLEL)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
END
!****************************************************************************************!
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: icpu, i_err
!icpu=0
!#if defined (_PARALLEL)
call MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint(*), i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: myreal(*)
INTEGER :: leng, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
!****************************************************************************************!

SUBROUTINE Sync_Procs
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER i_err
!#if defined (_PARALLEL)
call MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
IMPLICIT NONE
INCLUDE 'mpif.h'
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
!icpu=0
!#if defined (_PARALLEL)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal_in,myreal_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
REAL*8 :: myreal_in(*), myreal_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
!****************************************************************************************!

SUBROUTINE GlobSumI(myint_in,myint_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
INTEGER :: myint_in(*), myint_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
call MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
EOF
}
#---------------------------------------------------------------------------------------#
# Function to reweight Metadynamics Bias
function Write_Reweighting_Code() {
cat >>metad.f<<EOF
!---------------------------------------------------------------------!
!Written by Shalini Awasthi (ashalini@iitk.ac.in)
!---------------------------------------------------------------------!
      PROGRAM WSMTD_rw_2D
      IMPLICIT NONE
      REAL*8 gridmin1, gridmax1, griddif1, dummy2,v,
     &       gridmin2, gridmax2, griddif2,num,
     &       gridmin3, gridmax3, griddif3,
     &       gridmin4, gridmax4, griddif4,
     &       prob,den,alpha,fes,fes1,grid,Ro,
     &       cv1,cv2,ht,kt0,kt,ktb,bias_fact,ct,vbias,hill,width,
     &       diff_s2,ds2,ss,hh,dum,s1,s2,dummy11,cv3,cv4,s3,s4
      ALLOCATABLE cv1(:),cv2(:),ht(:),vbias(:),ct(:),hill(:),width(:),
     & prob(:,:),fes(:,:),fes1(:),grid(:),cv3(:),cv4(:)
      INTEGER mtd_steps,md_steps,dummy1,i,j,index1,index2,k,
     &        index3,index4,
     &        t_min,t_max,
     &        nbin1, nbin2,nbin3,nbin4,
     &        w_cv,w_hill,
     &        i_mtd, i_md, i_s2, i_s1,i_s3,i_s4,mtd_max,
     &        narg
      LOGICAL pmf,inpgrid
      CHARACTER*120 :: arg
      REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
      REAL*8, PARAMETER :: au_to_kcal = 627.51
      REAL*8, PARAMETER :: kj_to_kcal = 0.239006

!      character(len=50) f1

      OPEN(11,FILE='COLVAR',STATUS='unknown')
      OPEN(12,FILE='HILLS',STATUS='unknown')
!      OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')
!
      CALL get_steps(11,md_steps)
      CALL get_steps(12,mtd_steps)
!
      kt0=300.D0
      kt=300.D0
      bias_fact=1500.D0
      t_min=1
      t_max=md_steps
      pmf=.FALSE.
      inpgrid=.false.
      t_max=md_steps
      narg = IARGC()
      DO i=1,narg
        CALL GETARG(i,arg)
        IF(INDEX(arg,'-T0').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)kt0
        ELSEIF(INDEX(arg,'-T').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)kt
        ELSE IF(INDEX(arg,'-dT').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)bias_fact
        ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_min
        ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_max
           IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
        ELSE IF(INDEX(arg,'-pmf').NE.0)THEN
           pmf=.true.
        ELSE IF(INDEX(arg,'-grid').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)gridmin1
            CALL GETARG(i+2,arg)
            READ(arg,*)gridmax1
            CALL GETARG(i+3,arg)
            READ(arg,*)griddif1
            CALL GETARG(i+4,arg)
            READ(arg,*)gridmin2
            CALL GETARG(i+5,arg)
            READ(arg,*)gridmax2
            CALL GETARG(i+6,arg)
            READ(arg,*)griddif2
           inpgrid=.true.
        ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)w_cv
        ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)w_hill
        END IF
      END DO

      WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
      WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
      WRITE(*,'(A,F9.2)')'Physical Temp (K)      =',kt0
      WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
      WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
      WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
      WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
      WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
      WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
      IF(pmf)WRITE(*,'(A)')'PMF will be written in PMF.dat'
!
!      bias_fact=(kt+bias_fact)/kt
      bias_fact=(kt0+bias_fact)/kt0
      kt=kb*kt
      ktb=kt*bias_fact

!Make CV file containing CV1(t), CV2(t) from CONSTRAINT and cvmdck_mtd
      ALLOCATE(cv1(md_steps),cv2(md_steps),vbias(md_steps))
      ALLOCATE(cv3(md_steps),cv4(md_steps))
      ALLOCATE(ht(mtd_steps),ct(mtd_steps),hill(mtd_steps),
     &         width(mtd_steps))

      OPEN(16,file='cv.dat',status='unknown')
      DO i_md=1,md_steps
!        READ(11,*)dummy11,cv1(i_md),cv3(i_md),cv2(i_md),cv4(i_md)
        READ(11,*)dummy11,cv1(i_md),cv2(i_md)

          IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
          IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
          IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
          IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0
!          IF( cv3(i_md) .gt.  3.14d0)  cv3(i_md) = cv3(i_md) - 6.28d0
!          IF( cv3(i_md) .lt. -3.14d0 ) cv3(i_md) = cv3(i_md) + 6.28d0
!          IF( cv4(i_md) .gt.  3.14d0)  cv4(i_md) = cv4(i_md) - 6.28d0
!          IF( cv4(i_md) .lt. -3.14d0 ) cv4(i_md) = cv4(i_md) + 6.28d0


         WRITE(16,'(I10,4F16.6)')i_md,cv1(i_md),cv2(i_md)
!        WRITE(16,'(I10,4F16.6)')i_md,cv1(i_md),cv2(i_md),cv3(i_md),
!     &                         cv4(i_md)
      END DO
      WRITE(*,'(A)')'CV values written in cv.dat'
!
      IF(.NOT.inpgrid)
     &  CALL get_gridmin_max(16,gridmin1,gridmax1,griddif1,
     &                          gridmin2,gridmax2,griddif2)


      alpha=bias_fact/(bias_fact-1.D0)

      nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
      nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1
!      nbin3 = NINT((gridmax3-gridmin3)/griddif3)+1
!      nbin4 = NINT((gridmax4-gridmin4)/griddif4)+1
      WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
      WRITE(*,'(A10,3F8.4,I10)')'US  COORD:',
     &          gridmin1,gridmax1,griddif1,nbin1
      WRITE(*,'(A10,3F8.4,I10)')'MTD COORD:',
     &          gridmin2,gridmax2,griddif2,nbin2
!      WRITE(*,'(A10,3F8.4,I10)')'TAMD COORD:',
!     &          gridmin3,gridmax3,griddif3,nbin3
!      WRITE(*,'(A10,3F8.4,I10)')'TAMD COORD:',
!     &          gridmin4,gridmax4,griddif4,nbin4


      ALLOCATE(prob(nbin1,nbin2))
!      ALLOCATE(fes(nbin1,nbin2,nbin3,nbin4))
      ALLOCATE(fes(nbin1,nbin2))
      ALLOCATE(fes1(nbin2),grid(nbin2))

      DO i_mtd=1,mtd_steps
        READ(12,*) dummy11,hill(i_mtd),width(i_mtd),ht(i_mtd)
          IF( hill(i_mtd) .gt.  3.14d0) hill(i_mtd) = hill(i_mtd) - 6.28d0
          IF( hill(i_mtd) .lt. -3.14d0 )hill(i_mtd) = hill(i_mtd) + 6.28d0
        ht(i_mtd)=ht(i_mtd)*kj_to_kcal
      END DO

      CLOSE(11 12 13 14 16)

!calculate c(t)
      DO i_s2=1,nbin2 !grid over cv2 on which MTD is being done
        grid(i_s2)=gridmin2+DFLOAT(i_s2-1)*griddif2
      END DO

      OPEN(21,FILE='ct.dat',STATUS='unknown')
      DO i_mtd=1,mtd_steps
        ds2=width(i_mtd)*width(i_mtd)
        ss=hill(i_mtd)
        hh=ht(i_mtd)
        num=0.D0
        den=0.D0
        DO i_s2=1,nbin2
           diff_s2=grid(i_s2)-ss
            if (diff_s2 .gt. 3.14d0 ) diff_s2 =diff_s2 - 6.28d0
            if (diff_s2 .lt.-3.14d0 ) diff_s2 =diff_s2 + 6.28d0
           diff_s2=diff_s2*diff_s2*0.5D0
           fes1(i_s2)=fes1(i_s2)-hh*DEXP(-diff_s2/ds2)
           num=num+DEXP(-fes1(i_s2)/kt)
           den=den+DEXP(-fes1(i_s2)/ktb)
        END DO
        ct(i_mtd)=kt*DLOG(num/den)
        WRITE(21,'(I10,F16.8)')i_mtd,ct(i_mtd)
      END DO
      CLOSE(21)
      WRITE(*,'(A)')'CV values written in cv.dat'

!calculate v(s,t)
      DO i_md=1,md_steps
         mtd_max=(i_md*w_cv/w_hill)+1
         ss=cv2(i_md)
         dum=0.d0
         DO i_mtd=1,mtd_max
           ds2=width(i_mtd)*width(i_mtd)
           hh=ht(i_mtd)/alpha
           diff_s2=ss-hill(i_mtd)
            if (diff_s2 .gt. 3.14d0 ) diff_s2 =diff_s2 - 6.28d0
            if (diff_s2 .lt.-3.14d0 ) diff_s2 =diff_s2 + 6.28d0
           diff_s2=diff_s2*diff_s2*0.5D0
           dum=dum+hh*DEXP(-diff_s2/ds2)
         END DO
         vbias(i_md)=dum
      END DO

!calculate prob (unbiased from MTD potential)       
      den=0.d0
      prob=0.d0
      DO i_md=1,md_steps
        IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
          index1 = nint((cv1(i_md)-gridmin1)/griddif1) +1
          index2 = nint((cv2(i_md)-gridmin2)/griddif2) +1
!          index3 = nint((cv3(i_md)-gridmin3)/griddif3) +1
!          index4 = nint((cv4(i_md)-gridmin4)/griddif4) +1
          i_mtd=(i_md*w_cv/w_hill) + 1
          dum=vbias(i_md) - ct(i_mtd)
!          prob(index1,index2,index3,index4)=
!     &    prob(index1,index2,index3,index4)+DEXP(dum/kt)
          prob(index1,index2)=
     &    prob(index1,index2)+DEXP(dum/kt)
          den=den+DEXP(dum/kt)
        END IF
      END DO
!      dum=den*griddif1*griddif2*griddif3*griddif4
      dum=den*griddif1*griddif2
      den=1.D0/dum
      OPEN(2,FILE='Pu.dat',STATUS='unknown')
      DO i_s1=1,nbin1
        s1=DFLOAT(i_s1-1)*griddif1+gridmin1
        DO i_s2=1,nbin2
          s2=DFLOAT(i_s2-1)*griddif2+gridmin2
!          DO i_s3=1,nbin3
!            s3=DFLOAT(i_s3-1)*griddif3+gridmin3
!             DO i_s4=1,nbin4
!               s4=DFLOAT(i_s4-1)*griddif4+gridmin4
!          prob(i_s1,i_s2,i_s3,i_s4)=prob(i_s1,i_s2,i_s3,i_s4)*den
          prob(i_s1,i_s2)=prob(i_s1,i_s2)*den
!          WRITE(2,'(4E12.3,E16.8)')s1,s2,s3,s4,prob(i_s1,i_s2,i_s3,i_s4)
!          WRITE(2,'(E16.8)')prob(i_s1,i_s2,i_s3,i_s4)
          WRITE(2,'(3E16.8)')s1, s2, prob(i_s1,i_s2)
        END DO
        WRITE(2,*)
      END DO
!      END DO
!      END DO
      WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

      END PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
      SUBROUTINE get_steps(iunit,nsteps)
      IMPLICIT NONE
      INTEGER iunit, nsteps
      INTEGER ios
      nsteps=0
      REWIND(iunit)
      Read_Loop: DO
         READ(iunit,*,IOSTAT=ios)
         IF(ios.ne.0)EXIT Read_Loop
         nsteps=nsteps+1
      END DO Read_Loop
      REWIND(iunit)
      END
!---------------------------------------------------------------------!
EOF
}
#-----------------------------------------------------#
# Function for REwieting MEtaD
function Reweight_MetaD() {
mkdir TASS/ANALYSIS
cd TASS/ANALYSIS
mkdir PROB
Write_Reweighting_Code
gfortran metad.f -o metad.x
cat >> whaminput <<  EOF1
0.00001
EOF1

for i in `seq -3.2 0.2 3.2 `; do
cp ../TASS_$i/COLVAR .
cp ../TASS_$i/HILLS .
sed -i "/^#/d" COLVAR
sed -i "/^#/d" HILLS
a=`wc -l COLVAR |awk '{print $1}'`
b=`echo $(($a - 1))`
execute_histo $a
cat >>whaminput<<EOF
PROB/PROB_$i
$i  500  $b
EOF

cp Pu.dat PROB/PROB_$i
rm COLVAR HILLS cv.dat Pu.dat
done
mv whaminput PROB
rm metad.f 
cd ../..
}
#--------------------------------------------------#
# Run WHAM
function Run_WHAM() {
mkdir TASS/ANALYSIS/WHAM
cd TASS/ANALYSIS/WHAM
cp -rp ../PROB .
cp ../PROB/whaminput .
Write_WHAM_inputs
Write_WHAM_Code
mpif90  wham.f90 -o wham.x
rm wham.f90
./wham.x

gnuplot_2d_script
gnuplot plot_fes.gnp

rm -rf PROB
cd ../../..

}

#-----------------------------------------------------#
function Write_WHAM_inputs() {
cat >>input<<EOF
2 33 900
-3.14 3.14 0.10
-3.14 3.14 0.10
EOF
}
#------------------------------------------------------#
function execute_histo() {
./metad.x -T0 300 -T 900 -dT 900 -tmin 1 -tmax $1 -grid -3.14 3.14 0.10 -3.14 3.14 0.10 -pfrqMD 10 -dtMTD 500
}

#================================================================#
#     MAIN CODE STARTS FROM HERE  
#================================================================#
# Check for PDB file, gromacs & gmx_mpi and plumed
#check_for_pdb
check_for_gro_top
check_gromacs
check_plumed
#Prepare_initial_structure
#Perform_Minimisation
#Perform_Equilibration
#----------------
#Perform_US 
Perform_TASS
Reweight_MetaD
Run_WHAM

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
