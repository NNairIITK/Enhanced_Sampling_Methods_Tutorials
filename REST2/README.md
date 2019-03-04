
**Bash script to perform Replica Exchange Solute Scaling (REST2) simulation on alanine dipeptide in explicit solvet.**

* **Authour:**
   
      Anji Babu Kapakayala
      C/O Prof. Nisanth N. Nair
      Dept. of Chemistry
      IIT Kanpur, India.
       
                      
* **USAGE :**    
                         
      sh rest2_setup.sh                               
       
       
* **Requirements**:     
   
      * Latest versions of Plumed and Gromacs should be installed with mpi
      * Well equilibrated ala_wat.gro & ala_wat.top 
      * Gfortran, Gnuplot
                              
             
* **Description** :   
    
      On fly this script will generate required input files and performs the REST2 simulation for alanine 
      di-peptide in explicit solvent (Water in this case) using gromacs patched with plumed for 2 ns with
      the 5 replica. And the lambda values will be used in this tutorial are corresponds to the effective
      temperature range between 300 and 1000 K.
             
* **The following directories will be generated:**
             
      SCALED_TOPO       :   Contains the scaled topology of each lambda value
      REST2             :   Contains all the simulation files which run for 2ns each replica.
      ANALYSIS          :   Contains post processed files and plots.
             
           
* **Note:**
           
      This script uses the geometric progression to generate the effective temperature range for obtaining 
      the lamda values.
      
      
 * **REST2 TUTORIAL:**
 
 
 [Click me for the Tutorial](https://github.com/NNairIITK/Enhanced_Sampling_Methods_Tutorials/blob/master/Replica_Exchange_MD/REMD_Tutorial.pdf)
       
       
       
   **More Details will be updated Soon.......!!!!**
                
         Cheers!!!!
         Anji Babu
         23-02-2019
