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
cat >> nvt.mdp <<EOF
itle           = nvt simulation for rest2
integrator        = md           
nsteps            =20000     
dt                = 0.001 
; Output control
;nstxout          = 5000     
;nstvout          = 5000      
nstenergy         = 5000       
nstlog            = 5000  
nstxout-compressed  = 5000      
compressed-x-grps   = System  
; Bond parameters
continuation            = no          
constraint_algorithm    = lincs     ; holonomic constraints
constraints                 = all-bonds ; all bonds (even heavy atom-H bonds) constrained
lincs_iter                  = 1             ; accuracy of LINCS
lincs_order                 = 4             ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type             = grid              ; search neighboring grid cells
nstlist             = 10            ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb            = 1.0               ; short-range electrostatic cutoff (in nm)
rvdw                = 1.0               ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype         = PME               ; Particle Mesh Ewald for long-range electrostatics
pme_order           = 4             ; cubic interpolation
fourierspacing  = 0.16          ; grid spacing for FFT
; Temperature coupling is on
tcoupl          = V-rescale                 ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1     0.1           ; time constant, in ps
ref_t           = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is OFF
pcoupl = no
; Periodic boundary conditions
pbc             = xyz           ; 3-D PBC
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
; Velocity generation
gen_vel         = no            ; Velocity generation is off
EOF
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
#--------------------------------------------------------------#
# Preparing the input files for REST2 simulation (Topo files for defferent lamba)
function Prepare_Inputs() {
mkdir SCALED_TOPO
cd SCALED_TOPO
cp ../ala_wat.gro .
cp ../ala_wat.top .
write_mdp
Scale_Topology
Prepare_tpr_files
cd ..

}
#------------------------------------------#
# Function to make TPR files
function Prepare_tpr_files() {
j=0
list=`find *scaled*`
for i in $list;do
gmx_mpi  grompp  -f nvt.mdp -c ala_wat.gro -p ${i} -maxwarn 10 -o rest2_${i}.tpr
mv rest2_${i}.tpr rest2_${j}.tpr
j=$(($j+1))
done
rm \#*
}
#---------------------------------------------------------------------#
# SUBMIT REST2
function Perform_REST2() {
mkdir REST2
cd REST2
mv ../SCALED_TOPO/*.tpr .
write_plumed_input
# Submit
mpirun -np 5 gmx_mpi mdrun -v -deffnm rest2_ -plumed plumed.dat -multi 5 -replex 100 -hrex
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
   grep "hrex " temp &>/dev/null
        if [ $? -eq 0 ];then
        echo "Checking plumed: -hrex found"
        else
        echo "Checking plumed: -hrex NOT found"
	echo "Checking Plumed: Use latest version of Plumed for running HREX"
	fi

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
#-----------------------------------------------#
# Function to check the given topology file
function Check_Topology() {
if [  -f ./ala_wat.gro ]; then
 echo "Checking Gro: Found"
else    
echo "Checking Gro: gro file NOT found! Sorry."
exit
fi
#
if [  -f ./ala_wat.top ]; then
 echo "Checking Topology: Found"
grep -A10 " atoms " ala_wat.top |awk '{print $2}'|grep "_" &>/dev/null
#-----
	if [ $? -eq 0 ];then
	echo "Checking Topology: HOT atoms found"
	else
	echo "Checking Topology: Not found HOT atoms in Topology"
	echo "Checking Topology: Choose HOT atoms by adding (_) to the atomtype (i.e C_, H_. HC_ etc.) in the [ atoms ] section"
	echo "Thank you !"
	exit
	fi
#-------
else    
echo "Checking Topology: gromacs topology file NOT found! Sorry."
exit
fi
#
}
#-----------------------------------------------------------------#
#  Function to scale the Topologies
function Scale_Topology() {
# five replicas
nrep=5
# "effective" temperature range
tmin=300
tmax=1000
#
# build geometric progression
list=$(
awk -v n=$nrep \
    -v tmin=$tmin \
    -v tmax=$tmax \
  'BEGIN{for(i=0;i<n;i++){
    t=tmin*exp(i*log(tmax/tmin)/(n-1));
    printf(t); if(i<n-1)printf(",");
  }
}'
)
echo "$list" >> TEMP.dat
# clean directory
rm -fr \#*
rm -fr topol*

for((i=0;i<nrep;i++))
do

# choose lambda as T[0]/T[i]
# remember that high temperature is equivalent to low lambda
  lambda=$(echo $list | awk 'BEGIN{FS=",";}{print $1/$'$((i+1))';}')
echo "$lambda" >> LAMBDA_VALUES.dat
# process topology
plumed partial_tempering $lambda < ala_wat.top > ala_wat_scaled_$i.top
# clean the scaled topologyies ( anji)
sed -i '1,12 d' ala_wat_scaled_$i.top

done
}

#================================================================#
#     MAIN CODE STARTS FROM HERE  
#================================================================#
# Check for PDB file, gromacs & gmx_mpi and plumed

#check_for_pdb
check_gromacs
check_plumed
Check_Topology
Prepare_Inputs
Perform_REST2
#Analyse

#----------------

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
