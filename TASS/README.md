
# Temperature Accelerated Sliced Sampling (TASS)

**TASS TUTORIAL: (AMBER VERSION)**

[Ckick me for Tutorial](https://sites.google.com/view/the-nnn-group/tutorials/tass)
    
 
**Bash Script to perform TASS Simulation on Alanine dipeptide in vacuum using gromacs & Plumed**

* **Authour:**
            
            Anji Babu Kapakayala
            C/O Prof. Nisanth N. Nair
            Dept. of Chemistry
            IIT Kanpur, India.
                      
* **USAGE:**

            sh tass_setup.sh
             
* **Requirements:**

            Plumed and Gromacs should be installed with mpi
            gfortran, gnuplot, mpif90

* **INPUTS:**

            ala_di-pep.gro and ala_di-pep.top
                          
             
* **Description:**

            This script will generate required input files and runs the umbrella sampling simulation
            for alanine dipeptide in vacuum using gromacs and plumed (optional) for 2 ns to obtain the
		    initial structures for TASS simulations. And performs TASS then constructs the 2D free energy 
   		    along Phi and Psi collective variables
             
* **Directories:**

            On the fly, It will create following directories and stores related files in those directories:
             
                * MIN : write necessary input files and does the energy minimization.
                * NVT : Writes the required input files and performs equilibration here.
	            * US  : Performs the Umbrella Sampling Simulations for every CV value in their 
		            respective directory
	            * TASS:  Runs TASS Simulations in respective sub directories
	            * TASS/INPUTS: Stores all required input files here
	            * TASS/TASS_*: Respective directory for umbrella window        
	            * TASS/ANALYSIS: Writes all required inputs here to do analysis 
	            * TASS/PROB: Stores all the unbiased probabilities after reweighting metadynamics
	            * TASS/WHAM: Performs WHAM run here 


  **More Details will be updated soon..!!!**
 
 	
		Cheers !!!
 		Anji Babu
