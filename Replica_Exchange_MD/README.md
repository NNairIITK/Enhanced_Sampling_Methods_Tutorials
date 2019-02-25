
**Bash script to perform Replica Exchange Molecular Dynamics (t-REMD) simulation on alanine dipeptide in vacuum.**

   ..* **Authour:**
   
       Anji Babu Kapakayala
       C/O Prof. Nisanth N. Nair
       Dept. of Chemistry
       IIT Kanpur, India.
       
                      
   **USAGE :**    
                         
       sh remd_setup.sh                            [ default runs with plumed ]
       sh remd_setup.sh --without-plumed           [ Runs without Plumed i.e Uses only Gromacs ] 
       
       
   **Requirements**:     
   
       Latest versions of Plumed and Gromacs should be installed with mpi
       gfortran, gnuplot, ala_di-peptide.pdb
                          
             
   **Description** :   
   
       On fly this script will generate required input files and performs the remd simulation
       for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns with 4 replica.
             
   **The following directories will be generated:**
             
       INITIAL_STRUCTURES:   Contains the topology(top) & coordinate(gro) files which are made by using gmx2pdb.
       MIN               :   Contains files related to the energy minimization.
       REMD              :   Contains files after running REMD simulation for 2ns each replica.
       ANALYSIS          :   Contains post processed files and plots.
             
             
   **Note:**
           
       This script uses the geometric progression to generate the temperatures.
           
                    
   **More Details will be updated Soon.......!!!!**
                
         Cheers!!!!
         Anji Babu
         23-02-2019
