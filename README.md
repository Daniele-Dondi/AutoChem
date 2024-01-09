# AutoChem

                             ___         __        ________                 
                            /   | __  __/ /_____  / ____/ /_  ___  ____ ___ 
                           / /| |/ / / / __/ __ \/ /   / __ \/ _ \/ __ `__ \
                          / ___ / /_/ / /_/ /_/ / /___/ / / /  __/ / / / / /
                         /_/  |_\__,_/\__/\____/\____/_/ /_/\___/_/ /_/ /_/ 
                                            



                              Chemical Reaction Network Calculator

 This program generates all the products coming from the combination of a starting pool of compounds 
 and a set of reactions.
 The program can also download modules and pathways from KEGG database.
 Moreover, input files for ORCA and Gaussian can be produced for further calculations of free energies.
 For this aim, bash scripts present in the scripts directory help the automation process 
 for the submission of jobs, their checking and automatic restarting in order to have optimized 
 structures with no imaginary frequencies. 
 Moreover, scripts permits to calculate the reactions deltaG and to prepare a chemical network for 
 subsequent analysis.
 In addition, calculations can be saved and integrated, giving the possibility to mix calculations 
 on KEGG downloaded reactions and generated reactions or to make subsequent calculations with different
 reactions and/or reactants.

 IMPORTANT:
 The starting conditions (reactants and reactions) are defined in the function Set_Initial_Conditions. 
 Anyway, a basic editor is present in the program for inserting/enabling and disabling both reagents 
 and reactions.

 The program creates some supplementary files:
 COMPOUNDS.csv: A list of the compound number, its molecular formula and MW. 
     Energy at the rdkit level is also added if the optimization flag is enabled
 FREQUENCY.csv: A frequency distribution of molecular weights (aka: number of compounds having the same
     molecular weight)
 REACTIONS.txt: A list of all reactions calculated or downloaded from KEGG (useful for the determination of reactions 
     free energies)
 reactions.html: A html file containg the visual representation of reactions considered and calculated.
 SUMMARY.txt: A file containing summary informations for the calculations made. This is an important log
 file.
 smarts_reactions.txt: A text file containing all the reactions inserted by user in smarts format together with 
     auxiliary reagents and products
