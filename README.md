# AutoChem

                             ___         __        ________                 
                            /   | __  __/ /_____  / ____/ /_  ___  ____ ___ 
                           / /| |/ / / / __/ __ \/ /   / __ \/ _ \/ __ `__ \
                          / ___ / /_/ / /_/ /_/ / /___/ / / /  __/ / / / / /
                         /_/  |_\__,_/\__/\____/\____/_/ /_/\___/_/ /_/ /_/ 
                                            



                              Chemical Reaction Network Calculator

 This program generates all the products coming from the combination of a starting pool of compounds 
 and a set of reactions.<br>
  The program can also download modules and pathways from KEGG database.<br>
 Moreover, input files for ORCA, Gaussian, CREST and xTB can be produced for further calculations of free energies.
 For this aim, bash scripts present in the <i>scripts</i> directory help the automation process 
 for the submission of jobs, their checking and automatic restarting in order to have optimized 
 structures with no imaginary frequencies. <br>
 Scripts permit to calculate the reactions deltaG and to prepare a chemical network for 
 subsequent analysis.<br>
 In addition, calculations can be saved and iterated, giving the possibility to mix calculations 
 on KEGG downloaded reactions and generated reactions or to make subsequent calculations with different
 reactions and/or reactants.<br>
 AutoChem was published on the SoftwareX Journal https://doi.org/10.1016/j.softx.2024.101851<br>
 <br>
 IMPORTANT:
 The starting conditions (reactants and reactions) are defined in the plain text files starting_reactants.txt and starting_reactions.txt.<br>
 Anyway, a basic editor is present in the program for inserting, enabling and disabling both reagents 
 and reactions.

 The program creates some supplementary files:<br>
 <b><i>COMPOUNDS.csv</i></b>: A list of the compound number, its molecular formula and MW. Energy at the RDKit level is also added if the optimization flag is enabled<br>
 <b><i>FREQUENCY.csv</i></b>: A frequency distribution of molecular weights (aka: number of compounds having the same<br>
     molecular weight)<br>
 <b><i>REACTIONS.txt</i></b>: A list of all reactions calculated or downloaded from KEGG (useful for the determination of reactions 
     free energies)<br>
 <b><i>SUMMARY.txt</i></b>: A file containing summary informations for the calculations made. This is an important log
 file.<br>
 <b><i>reactions.html</i></b>: A html file containing the visual representation of reactions considered and calculated.<br>
 <b><i>smarts_reactions.txt</i></b>: A text file containing all the reactions inserted by user in SMARTS format together with 
     auxiliary reagents and products
<br> 
 ## Note
 AutoChem is divided in two programs:<br>
 <b>AutoChem</b>, written in Python, is managing the reaction generation<br>
 <b>check</b>, written in BASH, is managing the job submission with the automatic calculation of free energies<br>
<br>
# Installation
<br>
<b>AutoChem</b>: Download and install RdKit library in the proper Conda environment (follow instructions on RdKit website/anaconda installation)<br>
Copy the whole directory of AutoChem in your PC.<br>
Launch Autochem.py from your favourite IDE (we suggest Spyder)<br>
<br>
<br>
<b>check</b>: 
You should have awk, bc available on your Linux environment together with a task spooler like tsp or slurm<br>
If you use Windows for Autochem you should also install dos2unix<br>
To activate <i>check</i> for the first time you should:<br>
1) Copy the whole content of the folder <i>scripts</i> in your Linux working environment for calculations.<br>
2) Render the scripts runnable:<br>
chmod +x chmodme<br>
and then launch <i>chmodme</i><br>
./chmodme<br>
<br>
After doing this, you can put in the same folder the inputs to be calculated.<br>
Then launch check:<br>
./check<br>
and follow instructions.<br>
After some time launch check again and again up to the calculations of all inputs.<br>
<i>Check</i> will take care the resubmission of jobs in the case of negative frequencies.<br>
(if you are expert, you can add launch to chron or to some automatic running of check)<br>


# QUICKSTART
<br>Watch the video: (AutoChem program with SMARTS reactions)<br>
https://raw.githubusercontent.com/Daniele-Dondi/AutoChem/main/videos/AutoChem_example.mp4


https://github.com/user-attachments/assets/3b95d033-b533-40df-97f6-128ca5fede92


<br>
<br>Watch the video: (check program with crest)<br>
https://raw.githubusercontent.com/Daniele-Dondi/AutoChem/main/videos/check.mp4


https://github.com/user-attachments/assets/cb7bd15c-823a-44bd-8b3b-2234a6f599dc


<br>
<br>Watch the video: (check program with xtb and deltaG calculations)<br>
https://raw.githubusercontent.com/Daniele-Dondi/AutoChem/main/videos/xtb_+_deltaG.mp4


https://github.com/user-attachments/assets/003b1c4b-01f4-4b33-b956-cf628a6a7b8d



# SLURM setup<br>
In order to have the slurm queue manager working, you have to edit some files.<br>
<br>
1) slurm_settings. Edit the file inserting you username and accountname:<br>
user=your username<br>
account=your account<br>
<br>
2) Ask system administrators to have an example of slurm script for the program you are going to use<br>
as an example, we are going to use gaussian, so we create the file slurm_gaussian.<br>
In the file, substitute your username with $user and your accountname with $account.<br>
Substitute also the file name with $file, the memory allocation with $mem and the number of processors with $proc.<br>
These variables will be substituted when a file is submitted to the queue.<br>
<br>
3) Check if the slurm script is properly setup.<br>
To do this, edit the file run<br>
and locate the following code (locate the section related to the program you are going to use)<br>
  sbatch $file.tmp #2>/dev/null<br>
  #cat $file.tmp<br>
  rm $file.tmp<br><br>
modify as follows:<br>
  #sbatch $file.tmp #2>/dev/null<br>
  cat $file.tmp<br>
  #rm $file.tmp<br><br>
get ready with a simple example file to launch and write the command:<br>
./run yourtestfile<br>
(note that you should have prepared previously files like SLURM and GAUSSIAN too)<br>
check if the slurm header is well prepared:<br>
cat yourtestfile.tmp<br>
If it is OK, submit to the queue to check:<br>
sbatch yourtestfile.tmp<br>
4) If everything works, revert back to the previous version of run.<br>
5) Start your calculations with ./check and... good luck

<br>
<br>
