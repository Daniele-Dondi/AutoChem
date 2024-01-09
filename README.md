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
 Moreover, input files for ORCA and Gaussian can be produced for further calculations of free energies.
 For this aim, bash scripts present in the scripts directory help the automation process 
 for the submission of jobs, their checking and automatic restarting in order to have optimized 
 structures with no imaginary frequencies. <br>
 Moreover, scripts permits to calculate the reactions deltaG and to prepare a chemical network for 
 subsequent analysis.<br>
 In addition, calculations can be saved and integrated, giving the possibility to mix calculations 
 on KEGG downloaded reactions and generated reactions or to make subsequent calculations with different
 reactions and/or reactants.<br>

 IMPORTANT:
 The starting conditions (reactants and reactions) are defined in the files starting_reactants.txt and starting_reactions.txt.<br>
 Anyway, a basic editor is present in the program for inserting, enabling and disabling both reagents 
 and reactions.

 The program creates some supplementary files:<br>
 <b>COMPOUNDS.csv</b>: A list of the compound number, its molecular formula and MW. Energy at the rdkit level is also added if the optimization flag is enabled<br>
 FREQUENCY.csv: A frequency distribution of molecular weights (aka: number of compounds having the same<br>
     molecular weight)<br>
 REACTIONS.txt: A list of all reactions calculated or downloaded from KEGG (useful for the determination of reactions 
     free energies)<br>
 SUMMARY.txt: A file containing summary informations for the calculations made. This is an important log
 file.<br>
 reactions.html: A html file containing the visual representation of reactions considered and calculated.<br>
 smarts_reactions.txt: A text file containing all the reactions inserted by user in SMARTS format together with 
     auxiliary reagents and products
