**Bash script to perform Well Tempered metadynamics Simulation on alanine dipeptide in vacuum**

* **Authour:**

        Anji Babu Kapakayala
        C/O Prof. Nisanth N. Nair
        Dept. of Chemistry
        IIT Kanpur, India.
                      
* **USAGE:**

        sh metad_setup.sh
        
             
* **Requirements:**

        * Plumed and Gromacs should be installed with mpi
        * gfortran, gnuplot

* **INPUTS:**

        la_dipep.top & ala_dipep.gro
        
                       
* **Description:**

        On the fly, this script will generate required input files and runs the energy minimization,
        equilibration (NVT), Metadynamics for the alanine dipeptide in vacuum using gromacs and plumed 
        for 2 ns. Then, It will also does the basic analysis and plots the phi and ps with time, free energy 
        along phi , along psi and also along phi and psi (2D) and stores the plots in ANALYSIS directory.
             
* **Directories:**

            It will create following directories and stores related files in those directories:
             
         *  MIN : write necessary input files and does the energy minimization.
         *  NVT : Does the Equilibration for 2ns 
         *  METAD: Performs Metadynamics simulation and store all the related files.
         *  ANALYSIS: Contains the analysed plots in this directory.
             
             
  
               Cheers !!!
               Anji Babu
               March, 2019.


