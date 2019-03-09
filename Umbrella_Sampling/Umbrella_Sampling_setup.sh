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

Write_WHAM_inputs





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
#Analyse_runs
Make_Probabilities
#Run_WHAM

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
