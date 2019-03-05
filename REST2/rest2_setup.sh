#
#**Bash script to perform Replica Exchange Solute Scaling (REST2) simulation on alanine dipeptide in explicit solvet.**
#
#* **Authour:**
#   
#      Anji Babu Kapakayala
#     C/O Prof. Nisanth N. Nair
#      Dept. of Chemistry
#      IIT Kanpur, India.
#       
#                      
#* **USAGE :**    
#                         
#      sh rest2_setup.sh                               
#       
#       
# **Requirements**:     
#   
#      * Latest versions of Plumed and Gromacs should be installed with mpi
#      * Well equilibrated ala_wat.gro & ala_wat.top 
#      * Gfortran, Gnuplot
#                              
#             
#* **Description** :   
#    
#      On fly this script will generate required input files and performs the REST2 simulation for alanine 
#      di-peptide in explicit solvent (Water in this case) using gromacs patched with plumed for 2 ns with
#      the 5 replica. And the lambda values will be used in this tutorial are corresponds to the effective
#      temperature range between 300 and 1000 K.
#             
#* **The following directories will be generated:**
#             
#      SCALED_TOPO       :   Contains the scaled topology of each lambda value
#      REST2             :   Contains all the simulation files which run for 2ns each replica.
#      ANALYSIS          :   Contains post processed files and plots.
#             
#           
#* **Note:**
#           
#      This script uses the geometric progression to generate the effective temperature range for obtaining 
#      the lamda values.
#      
#      
# * **REST2 TUTORIAL:**
# 
# 
# [Click me for the Tutorial](https://github.com/NNairIITK/Enhanced_Sampling_Methods_Tutorials/blob/master/Replica_Exchange_MD/REMD_Tutorial.pdf)
#
#
#!/bin/bash
function write_mdp() {
cat >> nvt.mdp <<EOF
itle           = nvt simulation for rest2
integrator        = md           
nsteps            =2000000     
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
cd REST2
# Prints Avg. exchange probabilities into file
grep -A9 "average probabilities" *.log > ../ANALYSIS/Avg_Exchanges.dat
#
# Produce the plots of  potential energy overlap
# Uses gmx energy, gnuplot, and inhouse fortran code for calculating distribution
write_gnuplot_script
write_histogram_code
gfortran histogram.f90 -o hist.x
#-----PE
for i in '0 1 2 3 4';do
gmx_mpi energy -f rest2_${i}.edr -s rest2_${i}.tpr -o pe_${i}.xvg <<EOF
10

EOF
#----TEMP
gmx_mpi energy -f rest2_${i}.edr -s rest2_${i}.tpr -o temp_${i}.xvg <<EOF
14

EOF
sed -i "/^@/d" *.xvg
sed -i "/^#/d" *.xvg
mv pe_${i}.xvg pe_${i}.dat
mv temp_${i}.xvg temp_${i}.dat
awk '{print $2}' pe_${i}.dat > inputfile.dat
./hist.x
mv histogram.dat pe_hist_${i}.dat
awk '{print $2}' temp_${i}.dat > inputfile.dat
./hist.x
mv histogram.dat temp_hist_${i}.dat
done


# FES
plumed sum_hills --histo COLVAR.0 --idw psi --sigma 0.2 --kt 2.5 --outhisto fes_psi.dat
plumed sum_hills --histo COLVAR.0 --idw phi --sigma 0.2 --kt 2.5 --outhisto fes_phi.dat
gnuplot plot_fes_phi.gnp
gnuplot plot_fes_psi.gnp

gnuplot plot_psis.gnp
gnuplot plot_phis.gnp
gnuplot plot_phi-psi.gnp

cp *.png ../ANALYSIS/
cd ..
}
#-------------------------------------------------#
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
#-----------------------------------------------------------------#
function write_gnuplot_script() {
cat >>plot_PE.gnp<<EOF
reset
set term png #output terminal and file
set output "PE_overlap.png"
set border lw 2
#set tics font' ,15'
set ylabel " Probability Distribution" 
set xlabel " Potential Energy" 
set title "P.E Overlap REST2" font' ,15'
#Plot
plot "pe_hist_0.dat" u 1:2 w l lw 2 title'Lambda=1.0', "pe_hist_1.dat" u 1:2 w l lw 2 title'Lambda=0.7',"pe_hist_2.dat" u 1:2 w l lw 2 title'Lambda=0.5',"pe_hist_3.dat" u 1:2 w l lw 2 title'Lambda=0.4',"pe_hist_4.dat" u 1:2 w l lw 2 title'Lambda=0.3'
EOF
#--------------------------------------------#
cat >>plot_Temp.gnp<<EOF
reset
set term png #output terminal and file
set output "Temp_Overlap.png"
set border lw 2
#set tics font' ,15'
set ylabel " Probability Distribution" 
set xlabel " Temperature (K)" 
set title " Temp. Overlap REST2" font' ,15'
#Plot
plot "temp_hist_0.dat" u 1:2 w l lw 2 title'Lambda=1.0', "temp_hist_1.dat" u 1:2 w l lw 2 title'Lambda=0.7',"temp_hist_2.dat" u 1:2 w l lw 2 title'Lambda=0.5',"temp_hist_3.dat" u 1:2 w l lw 2 title'Lambda=0.4',"temp_hist_4.dat" u 1:2 w l lw 2 title'Lambda=0.3'
EOF
#------------------------------------------#
cat >>plot_phis.gnp<<EOF
set term png
set output 'ALL_PHI.png'
set border lw 2
set xrange [0:2000]
#set yrange [-3.14:3.4]
set xlabel "Time (ps)"
set ylabel "Phi (rad)"
#set title "REST2 Phi " font' ,20'
### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set label 1 'a' at graph 0.92,0.9 font ',8'
plot "COLVAR.0" u 1:2 w l title'Lamda=1.0'
#--- GRAPH b
set label 1 'b' at graph 0.92,0.9 font ',8'
plot "COLVAR.1" u 1:2 w l  title'Lamda=0.7'
# --- GRAPH c
set label 1 'c' at graph 0.92,0.9 font ',8'
plot "COLVAR.2" u 1:2 w l title'Lamda=0.5'
# --- GRAPH d
set label 1 'd' at graph 0.92,0.9 font ',8'
plot "COLVAR.4" u 1:2 w l title'Lamda=0.3'
unset multiplot
### End multiplot
EOF
#------------------------------------------------#
cat >>plot_psis.gnp<<EOF
set term png
set output 'ALL_PSI.png'
set border lw 2
set xrange [0:2000]
#set yrange [-3.14:3.4]
set xlabel "Time (ps)"
set ylabel "Psi (rad)"
#set title "REST2 Psi " font' ,20'
### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set label 1 'a' at graph 0.92,0.9 font ',8'
plot "COLVAR.0" u 1:3 w l title'Lamda=1.0'
#--- GRAPH b
set label 1 'b' at graph 0.92,0.9 font ',8'
plot "COLVAR.1" u 1:3 w l  title'Lamda=0.7'
# --- GRAPH c
set label 1 'c' at graph 0.92,0.9 font ',8'
plot "COLVAR.2" u 1:3 w l title'Lamda=0.5'
# --- GRAPH d
set label 1 'd' at graph 0.92,0.9 font ',8'
plot "COLVAR.4" u 1:3 w l title'Lamda=0.3'
unset multiplot
### End multiplot
EOF
#------------------------------------------------#
cat >>plot_phi-psi.gnp<<EOF
set term png
set output 'ALL_PHI_PSI.png'
set border lw 2
set xrange [0:2000]
#set yrange [-3.14:3.4]
set xlabel "Phi (rad)"
set ylabel "Psi (rad)"
#set title "REST2 Phi vs Psi " font' ,20'
### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
set label 1 'a' at graph 0.92,0.9 font ',8'
plot "COLVAR.0" u 2:3 w p title'Lamda=1.0'
#--- GRAPH b
set label 1 'b' at graph 0.92,0.9 font ',8'
plot "COLVAR.1" u 2:3 w p  title'Lamda=0.7'
# --- GRAPH c
set label 1 'c' at graph 0.92,0.9 font ',8'
plot "COLVAR.2" u 2:3 w p title'Lamda=0.5'
# --- GRAPH d
set label 1 'd' at graph 0.92,0.9 font ',8'
plot "COLVAR.4" u 2:3 w p title'Lamda=0.3'
unset multiplot
### End multiplot
EOF
#----------------------------------------#
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
}
#-----------------------------------------------#
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
Analyse_runs
#----------------

#================================================================#
#    Written By Anji Babu Kapakayala			         #
#================================================================#
