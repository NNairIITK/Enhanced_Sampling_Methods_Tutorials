#  Tutorial on Replica Exchange Molecular Dynamics (t-REMD) on Alanine dipeptide in vacuum.
#
#            Authour:  Anji Babu Kapakayala
#                      C/O Prof. Nisanth N. Nair
#                      Dept. of Chemistry
#                      IIT Kanpur, India.
#                      
#             USAGE:  sh remd_setup.sh
#             
#             Requirements:  Plumed and Gromacs should be installed with mpi
#                            gfortran, gnuplot
#                          
#             
#             Description:  This script will generate required input files and run the remd simulation
#             for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns.
#             
#             It will create following directories and stores related files in those directories:
#             
#             INITIAL_STRUCTURES:   Build the initial structures( gro & top) files
#             MIN : write necessary input files and does the energy minimization.
#             REMD: Prepare the input files for running REMD and Runs simulation for 2ns with 4 replica.
#             ANALYSIS: This script also does basic analysis and stores the plots in this directory.
#             
#             
#             Note: This script uses the geometric progression to generate the required temperatures. 
#
#
#
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
cat >>plumed.dat<<EOF
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

# Umbrella Bias
us: RESTRAINT ARG=phi KAPPA=500 AT=$1
# monitor the two variables
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR
EOF
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
# Perform US simulations
function Perform_US() {
mkdir US
cd US
mkdir INPUTS
write_mdp nvt 
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
write_plumed_input $i
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
#-------------------------------------------------#
#  Function for Analysing basic results
function Analyse_runs() {
mkdir ANALYSIS
cd REMD
# Prints Avg. exchange probabilities into file
grep -A9 "average probabilities" *.log > ../ANALYSIS/Avg_Exchanges.dat
#
# Concatenate all lof files into single log file (REMD.log)
# Which produces replica_index.xvg and replica_temp.xvg
cat *.log > REMD.log
demux.pl REMD.log
cp replica_*.xvg ../ANALYSIS/
# 
# Produces the continuos coordinate trajectories
gmx_mpi trjcat -f *.xtc -demux replica_index.xvg 
cp *_trajout.xtc ../ANALYSIS/
#
# Produce the plots of  potential energy overlap
# Uses gmx energy, gnuplot, and inhouse fortran code for calculating distribution
write_histogram_code
gfortran histogram.f90 -o hist.x
for i in 0 1 2 3;do
gmx_mpi energy -f remd_${i}.edr -s remd_${i}.tpr -o pe_${i}.xvg <<EOF
8

EOF
sed -i "/^@/d" pe_${i}.xvg
sed -i "/^#/d" pe_${i}.xvg
mv pe_${i}.xvg pe_${i}.dat
awk '{print $2}' pe_${i}.dat > inputfile.dat
./hist.x
mv histogram.dat pe_hist_${i}.dat
done
write_gnuplot_script
#
case "$1" in
--not_CVs)
gnuplot plot.gnp;;
*)
gnuplot plot_cv1.gnp
gnuplot plot.gnp
gnuplot plot_cv2.gnp
gnuplot plot_cv12.gnp
# FES
plumed sum_hills --histo COLVAR.0 --idw psi --sigma 0.2 --kt 2.5 --outhisto fes_psi.dat
plumed sum_hills --histo COLVAR.0 --idw phi --sigma 0.2 --kt 2.5 --outhisto fes_phi.dat
gnuplot plot_fes_phi.gnp
gnuplot plot_fes_psi.gnp

;;
esac
cp *.png ../ANALYSIS/
cd ..
}
#-------------------------------------------------#
#Function to plot histogram using gnuplot
function Plot_histogram() {
cat >>histogram_gnuplot.gnp<<EOF
reset
n=100 #number of intervals
max=-250. #max value
min=-350. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set term png #output terminal and file
set output "histogram.png"
set xrange [min:max]
#set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"
#count and plot
plot "data.dat" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle
EOF
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
# Plotting COLVARS with Time
cat >>plot_cv1.gnp<<EOF
reset
set term png #output terminal and file
set output "phi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Time (ps)" 
set title "REMD Simulation"
set ylabel " Phi" 
#Plot
plot "COLVAR.0" u 1:2 w l lw 2 title'Phi'
EOF
#
cat >>plot_cv2.gnp<<EOF
reset
set term png #output terminal and file
set output "psi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Time (ps)" 
set title "REMD Simulation" 
set ylabel " Psi" 
#Plot
plot "COLVAR.0" u 1:3 w l lw 2 title'Psi'
EOF
#
#
cat >>plot_fes_phi.gnp<<EOF
reset
set term png #output terminal and file
set output "fes_phi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Phi" 
set title "REMD Simulation" 
set ylabel " Free Energy" 
#Plot
plot "fes_phi.dat" u 1:2 w l lw 2 notitle
EOF
#
cat >>plot_fes_psi.gnp<<EOF
reset
set term png #output terminal and file
set output "fes_psi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Psi" 
set title "REMD Simulation" 
set ylabel " Free Energy" 
#Plot
plot "fes_psi.dat" u 1:2 w l lw 2 notitle
EOF
cat >>plot_cv12.gnp<<EOF
reset
set term png #output terminal and file
set output "psi_phi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Phi" 
set title "REMD Simulation" 
set ylabel " Psi" 
#Plot
plot "COLVAR.0" u 2:3 w l lw 2 notitle
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
cat >>wham_us.f90 <<EOF
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

#-----------------------------------------------------#
# Fortran code to create probabilities
function Write_histogram_Code() {
cat >>make_histograms.f<<EOF
!---------------------------------------------------------------------!
!Written by Shalini Awasthi (ashalini@iitk.ac.in)
!---------------------------------------------------------------------!
      PROGRAM WSMTD_rw_2D
      IMPLICIT NONE
      REAL*8 gridmin1, gridmax1, griddif1, dummy2,dummy3,v,
     &       gridmin2, gridmax2, griddif2,
     &       prob,den,alpha,fes,fes1,grid,Ro,
     &       cv1,cv2,kt0,kt,ktb,
     &       dum,s1,s2,dummy11
      ALLOCATABLE cv1(:),cv2(:),prob(:,:),fes(:,:),fes1(:),grid(:)
      INTEGER md_steps,dummy1,i,j,index1,index2,k,
     &        t_min,t_max,nbin1,nbin2,i_md,i_s1,i_s2,narg,w_cv
      LOGICAL pmf,inpgrid
      CHARACTER*120 :: arg
      REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
      REAL*8, PARAMETER :: au_to_kcal = 627.51
      REAL*8, PARAMETER :: kj_to_kcal = 0.239006

!      character(len=50) f1

      OPEN(11,FILE='COLVAR',STATUS='unknown')
      OPEN(12,FILE='cv.dat',STATUS='unknown')
!      OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')
!
      CALL get_steps(11,md_steps)
!

      kt0=300.D0
      kt=300.D0
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
        ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_min
        ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_max
           IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
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
        END IF
      END DO

!      md_steps=t_max
      WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
      WRITE(*,'(A,I10)')'No: of max MD  steps        =',t_max
      WRITE(*,'(A,I10)')'No: of min MD  steps        =',t_min

      ALLOCATE(cv1(md_steps),cv2(md_steps))

      DO i_md=1,md_steps
        READ(11,*)dummy11,cv1(i_md),cv2(i_md)





        IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
        IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
        IF( cv2(i_md) .gt.  3.14d0)  cv2(i_md) = cv2(i_md) - 6.28d0
        IF( cv2(i_md) .lt. -3.14d0 ) cv2(i_md) = cv2(i_md) + 6.28d0

        WRITE(12,*)dummy11,cv1(i_md),cv2(i_md)
      END DO
      WRITE(*,'(A)')'CV values written in cv.dat'

      nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
      nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1
      WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
      WRITE(*,'(A10,3F8.4,I10)')'US  COORD:',
     &          gridmin1,gridmax1,griddif1,nbin1
      WRITE(*,'(A10,3F8.4,I10)')'RE',
     &          gridmin2,gridmax2,griddif2,nbin2
      ALLOCATE(prob(nbin1,nbin2))
      ALLOCATE(fes(nbin1,nbin2))
!calculate prob        
      den=0.d0
      prob=0.d0
      DO i_md=1,md_steps
      if ((i_md.ge.t_min).and.(i_md.le.t_max)) then
          index1 = nint((cv1(i_md)-gridmin1)/griddif1) +1
          index2 = nint((cv2(i_md)-gridmin2)/griddif2) +1
          prob(index1,index2)=prob(index1,index2)+ 1.d0
      end if
      END DO

      DO index1=1,nbin1
        DO index2=1,nbin2
         den=den+prob(index1,index2)
        end do
      end do

!      dum=dfloat(t_max-t_min+1)!*griddif1
!      write(*,*) den, dum
!      dum=dfloat(den)*griddif1

      OPEN(2,FILE='Pu.dat',STATUS='unknown')
        DO i_s1=1,nbin1
           s1=DFLOAT(i_s1-1)*griddif1+gridmin1
           DO i_s2=1,nbin2
            s2=DFLOAT(i_s2-1)*griddif2+gridmin2
!           WRITE(2,'(2E16.8)')s1,s2,prob(i_s1,i_s2)
           WRITE(2,*)s1,s2,prob(i_s1,i_s2)/(den*griddif1*griddif2)
        END DO
        WRITE(2,*)
      END DO
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
#--------------------------------------------------#
# Function 
function Make_Probabilities() {
mkdir US/ANALYSIS
cd US/ANALYSIS
mkdir PROB

Write_histogram_Code
gfortran make_histograms.f -o histo.x
cat >> whaminput <<  EOF1
0.00001
EOF1

for i in `seq -3.2 0.2 3.2 `; do
cp ../US_$i/COLVAR .
sed -i "/^#/d" COLVAR
a=`wc -l COLVAR |awk '{print $1}'`
b=`echo $(($a - 1))`
execute_histo $a
cat >>whaminput<<EOF
PROB/PROB_$i
$i  500  $b
EOF
rm COLVAR cv.dat
cp Pu.dat PROB/PROB_$i
done
mv whaminput PROB
rm make_histograms.f 
cd ../..
}
#--------------------------------------------------#
# Run WHAM
function Run_WHAM() {
mkdir US/ANALYSIS/WHAM
cd US/ANALYSIS/WHAM
cp -rp ../PROB .
cp ../PROB/whaminput .
Write_WHAM_inputs
Write_WHAM_Code
mpif90  wham_us.f90 -o wham.x
rm wham_us.f90
./wham.x

gnuplot_2d_script
gnuplot plot_fes.gnp

rm -rf PROB
cd ../../..

}

#-----------------------------------------------------#
function Write_WHAM_inputs() {
cat >>input<<EOF
2 33 300
-3.14 3.14 0.10
-3.14 3.14 0.10
EOF
}
#------------------------------------------------------#
function execute_histo() {
./histo.x -T 300 -tmin 1 -tmax $1 -grid -3.2 3.2 0.05 -3.2 3.2 0.05 -pfrqMD 10
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
Perform_Minimisation
Perform_Equilibration
#----------------
Perform_US 
Make_Probabilities
Run_WHAM

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
