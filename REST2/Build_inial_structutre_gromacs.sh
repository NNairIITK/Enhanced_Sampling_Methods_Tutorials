# Bash Script to build the initial structure with explicit solvent
#
#  Authour: Anji Babu Kapakayala
#		IIT Kanpur, India

# Usage: sh Build_inial_structutre_gromacs.sh <pdbfile>
#
#
#

#!/bin/bash

function build_structure() {

gmx pdb2gmx -f $1 -o ala.gro  -water tip3p <<EOF
5
EOF
gmx editconf -f ala.gro -o ala_box.gro -c -d 1.0 -bt cubic
gmx solvate -cp ala_box.gro -cs spc216.gro -o ala_solv.gro -p topol.top
write_mdp
gmx grompp -f ions.mdp -c ala_solv.gro -p topol.top -o ions.tpr -pp ala_processed.top
gmx genion -s ions.tpr -o ala_solv_ions.gro -p ala_processed.top -pname NA -nname CL -nn 8 <<EOF
13
EOF
}
#-------------------------#
function write_mdp() {
cat >>ions.mdp<<EOF
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator      = steep         ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps          = 50000         ; Maximum number of (minimization) steps to perform
; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist             = 1             ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type             = grid              ; Method to determine neighbor list (simple, grid)
coulombtype         = PME               ; Treatment of long range electrostatic interactions
rcoulomb            = 1.0               ; Short-range electrostatic cut-off
rvdw                = 1.0               ; Short-range Van der Waals cut-off
pbc                     = xyz           ; Periodic Boundary Conditions (yes/no)
EOF
}


#------------------
build_structure $1
