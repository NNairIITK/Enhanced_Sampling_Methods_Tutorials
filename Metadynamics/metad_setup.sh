#  Bash script to perform Well Tempered metadynamics Simulation on alanine dipeptide in vacuum.
#
#            Authour:  Anji Babu Kapakayala
#                      C/O Prof. Nisanth N. Nair
#                      Dept. of Chemistry
#                      IIT Kanpur, India.
#                      
#             USAGE:  sh metad_setup.sh
#             
#             Requirements:  Plumed and Gromacs should be installed with mpi
#                            gfortran, gnuplot
#
#	      INPUTS: ala_dipep.top & ala_dipep.gro
#                          
#             
#             Description:  On the fly, this script will generate required input files and runs the energy minimization, equilibration (NVT), Metadynamics 
#             for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns. Then, It will also does the basic analysis and plots the phi and ps with time,
#	      free energy along phi , along psi and also along phi and psi (2D) and stores the plots in ANALYSIS directory.
#             
#
#             It will create following directories and stores related files in those directories:
#             
#             MIN : write necessary input files and does the energy minimization.
#             NVT : Does the Equilibration for 2ns 
#	      METAD: Performs Metadynamics simulation and store all the related files.
#             ANALYSIS: Contains the analysed plots in this directory.
#             
#             
#  
# 		Cheers !!!
#		Anji Babu
#		March, 2019.
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
nsteps		= 500000	  	; Maximum number of (minimization) steps to perform

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
ref_t                    = 300.0
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

remd)
cat >>plumed.dat<<EOF
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17

# monitor the two variables
PRINT STRIDE=10 ARG=phi,psi FILE=COLVAR
EOF
;;
metad)
cat >>plumed.dat<<EOF
# This input file has taken from below reference
#https://plumed.github.io/doc-v2.5/user-doc/html/belfast-6.html
#
# set up two variables for Phi and Psi dihedral angles 
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
#
# Activate metadynamics in phi and psi
# depositing a Gaussian every 500 time steps,
# with height equal to 1.2 kJoule/mol,
# and width 0.35 rad for both CVs.
# Well-tempered metadynamics is activated,
# and the bias factor is set to 6.0
#
metad: METAD ARG=phi,psi PACE=500 HEIGHT=2.4 SIGMA=0.5,0.5 FILE=HILLS BIASFACTOR=6.0 TEMP=300.0
#metad: METAD ARG=phi,psi PACE=500 HEIGHT=2.4 SIGMA=0.5,0.5 FILE=HILLS 
#
# monitor the two variables and the metadynamics bias potential
PRINT STRIDE=10 ARG=phi,psi,metad.bias FILE=COLVAR
EOF
;;
*) echo "Invalid argument with write_plumed_input function"
;;
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
#cp ../INITIAL_STRUCTURES/* .
cp ../*.top .
cp ../*.gro .
write_mdp min
gmx_mpi grompp -f min.mdp -c *.gro -p *.top -o min.tpr
gmx_mpi mdrun -deffnm min -v
cd ..
}
#--------------------------------------------------------------#
function Perform_NVT() {
mkdir NVT
cd NVT
cp ../MIN/min.gro .
cp ../MIN/*.top .
#cp ../INITIAL_STRUCTURES/topol.top .
#cp ../INITIAL_STRUCTURES/*.itp .
write_mdp nvt
gmx_mpi grompp -f nvt.mdp -c min.gro -p *.top -o nvt.tpr
gmx_mpi mdrun -deffnm nvt -v
cd ..
}
#--------------------------------------------------------------#
# Prepare the input files for METAD 
function Prepare_inputs_for_metad() {
mkdir METAD
cd METAD
cp ../NVT/nvt.gro .
cp ../NVT/*.top .
#cp ../INITIAL_STRUCTURES/topol.top .
#cp ../INITIAL_STRUCTURES/*.itp .
write_mdp nvt
write_plumed_input metad
cd ..
}
#------------------------------------------#

#------------------------------------------#
function Perform_METAD() {
#CREATE TPR files
cd METAD
gmx_mpi grompp -f nvt.mdp -c nvt.gro -p *.top -o metad.tpr -maxwarn 10
rm \#*
# SUBMIT METAD with PLUMED
mpirun -np 1 gmx_mpi mdrun -v -deffnm metad -plumed plumed.dat 
cd ..
}
#-------------------------------------------------#
#  Function for Analysing basic results
function Analyse_runs() {
mkdir ANALYSIS
cd METAD
#
write_gnuplot_script
#
gnuplot plot_cv1.gnp
gnuplot plot_cv2.gnp
gnuplot plot_hills.gnp
#FES using sum_hilles
plumed sum_hills --hills HILLS --outfile fes_2D.dat
#plumed sum_hills --histo COLVAR --sigma 0.2,0.2 --kt 2.5 --outhisto fes_2D.dat
plumed sum_hills --hills HILLS --idw phi --kt 2.5 --outfile fes_phi.dat
#plumed sum_hills --histo COLVAR --idw phi --kt 2.5 --outfile fes_phi.dat
#plumed sum_hills --histo COLVAR --idw psi --kt 2.5 --outfile fes_psi.dat
plumed sum_hills --hills HILLS --idw psi --kt 2.5 --outfile fes_psi.dat
gnuplot plot_fes1.gnp
gnuplot plot_fes2.gnp
gnuplot plot_fes.gnp

mv *.png *.dat *.gnp ../ANALYSIS/

cd ..
}
#-------------------------------------------------#
function write_gnuplot_script() {
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
splot 'fes_2D.dat' us 1:2:3 with pm3d
EOF
#
cat >>plot_fes1.gnp<<EOF
reset
set term png #output terminal and file
set output "fes_phi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Phi" 
set title "MEATD Simulation"
set ylabel " Free Energy" 
#Plot
plot "fes_phi.dat" u 1:2 w l lw 2 notitle
EOF
#
cat >>plot_fes2.gnp<<EOF
reset
set term png #output terminal and file
set output "fes_psi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Psi" 
set title "MEATD Simulation"
set ylabel " Free Energy" 
#Plot
plot "fes_psi.dat" u 1:2 w l lw 2 notitle
EOF
#
# Plotting COLVARS with Time
cat >>plot_cv1.gnp<<EOF
reset
set term png #output terminal and file
set output "phi.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Time (ps)" 
set title "MEATD Simulation"
set ylabel " Phi" 
#Plot
plot "COLVAR" u 1:2 w l lw 2 title'Phi'
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
set title "METAD Simulation" 
set ylabel " Psi" 
#Plot
plot "COLVAR" u 1:3 w l lw 2 title'Psi'
EOF
#
cat >>plot_hills.gnp<<EOF
reset
set term png #output terminal and file
set output "Hill_hieght.png"
set border lw 2
set grid
#set tics font' ,15'
set xlabel " Time(ps)" 
set title "METAD Simulation" 
set ylabel " Gaussian Height" 
#Plot
plot "HILLS" u 1:6 w l lw 2 title'Height of Gaussian'
EOF
}
#-------------------------------------------------------#
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
#================================================================#
#     MAIN CODE STARTS FROM HERE  
#================================================================#
# Check for PDB file, gromacs & gmx_mpi and plumed
check_for_pdb
check_gromacs
check_plumed
#Prepare_initial_structure
Perform_Minimisation
Perform_NVT
Prepare_inputs_for_metad
Perform_METAD
Analyse_runs
#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
