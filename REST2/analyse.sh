#!/bin/bash
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
#write_gnuplot_script
#
#cp *.png ../ANALYSIS/
cd ..
}
#-------------------------------------------------#

