**  Bash Script to perform Umbrella Sampling Simulation on Alanine dipeptide in vacuum.**

* **Authour:**
  
        Anji Babu Kapakayala
        C/O Prof. Nisanth N. Nair
        Dept. of Chemistry
        IIT Kanpur, India.
                    
* **USAGE:**

        sh Umbrella_Sampling_setup.sh
             
* **Requirements:**
  
        Plumed and Gromacs should be installed with mpi
        gfortran, gnuplot, mpif90
#
#             INPUTS:  ala_di-pep.gro and ala_di-pep.top
#                          
#             
#             Description:  This script will generate required input files and run the umbrella sampling simulation
#             for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns. And constructs the 2D free energy along Phi and Psi collective variables
#             
#             On the fly, It will create following directories and stores related files in those directories:
#             
#             MIN : write necessary input files and does the energy minimization.
#             NVT : Writes the required input files and performs equilibration here.
#             US  : Performs the Umbrella Sampling Simulations here in respective directory for every CV value
#             US/INPUTS: Stores all required input files here        
#             US/ANALYSIS: Writes all required inputs here to do analysis 
#             US/PROB: Stores all the  probabilities corresponding to all the umbrella windows here
#             US/WHAM: Performs WHAM run here 
#
#
#
#            Cheers !!
#           Anji Babu
#
#!/bin/bash
