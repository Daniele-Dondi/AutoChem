Linux BASH utilities for AutoChem package
The directory scripts contains BASH scripts for the automatic DFT/ab-initio calculations with Gaussian and ORCA.
It is required a Linux system having a queue manager like tsp or SLURM.
The main script is called 'check'.
The program will automatically submit input files to the processing queue at the first launch.
Input files should have optimization and frequency jobs.
When it is launched subsequently, it will check if the calculations are finished.
If they are finished, it will check if they finished by errors or negative frequency.
In that case, it will automatically resubmit again up to completion.
At the end of all calculations, if REACTIONS.txt is given, it will automatically calculate the reactions deltaG.

For the installation:
1 copy all files from the directory scripts to your working folder
2 render chmodme executable (chmod +x chmodme or sudo chmod +x chmodme, depending on your Linux flavor)
3 launch chmodme (./chmodme or sudo ./chmodme). It will automatically render other scripts executable
4 copy all your input files (together with REACTIONS.txt, if needed) to your working folder
5 indicate which type of queue you are using
  if tsp, type touch TSP, if SLURM type touch SLURM
6 indicate which software you are using
  for Gaussian, type touch GAUSSIAN, if orca type touch ORCA
7 check that file calculation header is corresponding to the job of your input files (open gauheader or orcaheader)
8 launch check by typing ./check
9 after some time repeat step 8 up to the termination of all calculations

NOTES: for SLURM users, please consider to read the paper explaining how to setup properly SLURM scripts.


