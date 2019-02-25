#  Tutorial on Replica Exchange Molecular Dynamics (t-REMD) on Alanine dipeptide in vacuum.
#
#            Authour:  Anji Babu Kapakayala
#                      C/O Prof. Nisanth N. Nair
#                      Dept. of Chemistry
#                      IIT Kanpur, India.
#                      
#             USAGE:  sh remd_setup.sh     [ default runs with plumed ]
#                     sh remd_setup.sh --without-plumed   [ Runs without plumed i.e Uses only Gromacs]
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
cat >> nvt_$2.mdp <<EOF
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
ref_t                    = $3
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
cp ../INITIAL_STRUCTURES/* .
write_mdp min
gmx_mpi grompp -f min.mdp -c conf.gro -p topol.top -o min.tpr
gmx_mpi mdrun -deffnm min -v
cd ..
}
#--------------------------------------------------------------#
# Prepare the input files for REMD simulation (NVT files for defferent Temp)
function Prepare_inputs_for_remd() {
mkdir REMD
cd REMD
cp ../MIN/min.gro .
cp ../INITIAL_STRUCTURES/topol.top .
cp ../INITIAL_STRUCTURES/*.itp .
#------------------------------------------#
awk 'BEGIN { T0=300.00; c=0.2; N=4; for (i=1; i<=N; ++i) { print T0; T0=T0*exp(i*c) }}' > temp.dat
TEMP_FILE="temp.dat"
#
#===== Writing multiple MDP files====>>>
#
j=0
while read temp;do
#echo "$temp"
write_mdp nvt $j $temp
j=$(($j+1))
done < $TEMP_FILE
}
#------------------------------------------#
function Perform_REMD() {
#CREATE TPR files
cd REMD
for i in 0 1 2 3;do
gmx_mpi grompp -f nvt_${i}.mdp -c min.gro -p topol.top -o remd_${i}.tpr -maxwarn 10
done
rm \#*
# SUBMIT REMD with PLUMED
# Submit with Plumed
case "$1" in
plumed|with_plumed)
#write plumed.dat
write_plumed_input
mpirun -np 4 gmx_mpi mdrun -v -deffnm remd_ -plumed plumed.dat -multi 4 -replex 100
#gmx_mpi mdrun -v -deffnm md_ -plumed plumed.dat -multi 2 -replex 50
;;
*)
mpirun -np 4 gmx_mpi mdrun -v -deffnm remd_ -multi 4 -replex 100
;;
esac
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
gnuplot plot.gnp
cp All_histograms.png ../ANALYSIS/
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
}

#-------------------------------------------------------#
# Function to produce fortran code to calculate histogram for given data
function  write_histogram_code() {
cat >>histogram.f90 <<EOF
!FOTRAN90 program to calculate Distribution of given data.
!
!Authour: ANJI BABU KAPAKAYALA
!         IIT KANPUR, INDIA.
!
!Writing Histogram 
PROGRAM histogram
  implicit none
  real*8              ::x,xmin,xmax,width,u,t
  real*8, allocatable ::prob(:)
  integer*8           ::i,j,N,bin,nbin  !N=md_steps
  character (len=125) :: input_filename,output_filename
  open(12,file="inputfile.dat",STATUS="OLD")
  open(13,file="histogram.dat",STATUS="NEW")
  call get_steps(12,N)
  print*,"Total steps =",N
  call get_xmin_xmax(12,xmin,xmax,width)
  print*, "xmin= ",xmin,"xmax=",xmax,"Width=",width
  nbin=nint(((xmax-xmin)/width)) + 1
  print *, 'nbin =',nbin
  allocate(prob(nbin+1))
  prob = 0.0
     do i=1,N
        read(12,*)x
        bin =int((x-xmin)/width)+1
           if(bin < nbin+1)then
             prob(bin) = prob(bin) +1.00
           end if
     end do
     do j=1,nbin
        write(13,*) real(j-1)*width+xmin, prob(j)/dfloat(N) !printing x vs P(x) vaules
     end do
  close(12)
  close(13)
  deallocate(prob)
END PROGRAM histogram
!===========================================================================!
SUBROUTINE get_steps(iunt,nsteps)
IMPLICIT NONE
INTEGER,INTENT(IN)  ::iunt
INTEGER*8,INTENT(OUT)  ::nsteps
INTEGER ::ios
nsteps=0
REWIND(iunt)
DO
READ(iunt,*,IOSTAT=ios)
IF(ios /= 0) EXIT
nsteps=nsteps+1
ENDDO
REWIND(iunt)
END SUBROUTINE get_steps
!============================================================================!
SUBROUTINE get_xmin_xmax(iunt,xmin,xmax,xdif)
 IMPLICIT NONE
 INTEGER,INTENT(IN)  :: iunt
 REAL*8,INTENT(INOUT)  :: xmin,xmax,xdif
 INTEGER :: ios
 REAL*8  :: x
 INTEGER, PARAMETER :: Def_No_bins=101
 REWIND(iunt)
 READ(iunt,*,IOSTAT=ios)x
 if(ios /= 0)stop 'ERROR reading INPUT'
 xmin=x
 xmax=x
 RLoop: DO
   READ(iunt,*,IOSTAT=ios)x
   if(ios /= 0)EXIT RLoop
   xmin=MIN(xmin,x)
   xmax=MAX(xmax,x)
 END DO RLoop
 xdif=(xmax-xmin)/DFLOAT(Def_No_bins)
 REWIND(iunt)
END SUBROUTINE get_xmin_xmax
!============================================================================!
!                        ANJI BABU KAPAKAYALA                                !
!============================================================================!                                                         
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
    echo "Exiting..."
    exit
  fi
rm temp
else
echo "Checking plumed: Plumed NOT found"
echo "Exiting...!"
exit
fi
}
#-----------------------------------------------#
function check_gromacs() {
gmx_mpi mdrun -h &> /dev/null
if [ $? -eq 0 ];then
echo "Checking Gromacs: Found"
else
echo "Checking Gromacs: Gromacs NOT found"
echo "Exiting..!"
exit
fi
}
#----------------------------------------------#
#================================================================#
#     MAIN CODE STARTS FROM HERE  
#================================================================#
# Check for gromacs & gmx_mpi
check_gromacs
check_plumed
Prepare_initial_structure
Perform_Minimisation
Prepare_inputs_for_remd
#----------------
case "$1" in
--without-plumed)
Perform_REMD ;;
*)
Perform_REMD plumed
;;
esac
Analyse_runs

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#



