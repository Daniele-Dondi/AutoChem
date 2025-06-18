# -*- coding: utf-8 -*-
#
#
#                             ___         __        ________                 
#                            /   | __  __/ /_____  / ____/ /_  ___  ____ ___ 
#                           / /| |/ / / / __/ __ \/ /   / __ \/ _ \/ __ `__ \
#                          / ___ / /_/ / /_/ /_/ / /___/ / / /  __/ / / / / /
#                         /_/  |_\__,_/\__/\____/\____/_/ /_/\___/_/ /_/ /_/ 
#                                                                                              
# Version 2.0 2025
#
# Big text made by using https://www.fancytextpro.com/BigTextGenerator/
#
#                              Chemical Reaction Network Calculator
#
# This program can operate in two different modes:
# SMARTS mode: AutoChem generates all the products coming from the combination of a starting pool of compounds 
# and a set of reactions acting as a virtual chemical reactor. (See the section IMPORTANT below)
# KEGG mode: AutoChem downloads modules or pathways from KEGG database
# There is also the possibility to read data created by previous jobs and continue from that point. 
# This allows the creation of mixed jobs.
# Moreover, input files for ORCA, Gaussian, CREST and xTB can be produced for further calculations of free energies.
# 
# AutoChem possesses a second module called check, an intelligent script to launch, control and
# retrieve data from Gaussian and ORCA calculations (now extended to XTB and CREST).
# For this aim, bash scripts present in the scripts directory help the automation process 
# for the submission of jobs, their checking and automatic restarting in order to have optimized 
# structures with no imaginary frequencies. 
# Moreover, scripts permit to calculate the reactions deltaG and to prepare a chemical network for 
# subsequent analysis.
# The second module check is not dependent from AutoChem and it is able to manage multiple input files 
# created by other methods too.
# Refer to README.txt in scripts folder to have furher information about this module.
#
# IMPORTANT:
# The starting conditions (reactants and reactions) for the calculations of SMARTS reactions are read
# from files starting_reactants.txt and smarts_reactions.txt
# A basic editor is present in the program for inserting/enabling and disabling both reagents 
# and reactions.
# We suggest to edit directly these files for customisation.
#
# The program creates some supplementary files:
# COMPOUNDS.csv: A list of the compound number, its molecular formula and MW. 
#     Energy at the RDKit level is also added if the optimization flag is enabled
# FREQUENCY.csv: A frequency distribution of molecular weights (aka: number of compounds having the same
#     molecular weight)
# REACTIONS.txt: A list of all reactions calculated or downloaded from KEGG (useful for the determination of reactions 
#     free energies)
# SUMMARY.txt: A file containing summary informations for the calculations made. This is an important log
# file.
# reactions.html: A html file containing the visual representation of reactions considered and calculated.
#
# The program reads some files (in the case of SMARTS reactions):
# starting_ractants.txt: A text file containing all the starting reactants in SMILES format
# smarts_reactions.txt: A text file containing all the reactions inserted by user in SMARTS format together with 
#     auxiliary reagents and products

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import os
import glob
import time
import datetime
import shutil
import urllib.request #for KEGG
import ssl

ssl._create_default_https_context = ssl._create_stdlib_context


#global vars
c_names=[] #used if rdkit energies are calculated (to locate the correct energy from the molecule name)
pool=[] #database of all reagents in SMILES format
reactant_is_enabled=[] #reactants can be present but not enabled for the calculations
manually_added_reactants=[]
Formulas=[] #database for all molecular formulas in human-readable format
MWs=[] #database of all reagents molecular weights
energies=[] #databases of all rdkit energies (if calculated)
reactions=[] #database of all reactions for rdkit
reactions_smarts=[] #database of all reactions in smarts format
reaction_is_enabled=[] #reactions can be present but not enabled for the calculations
num_reagents=[] #database indicating the number of reagents for each reaction
aux_reagents=[] #database containing the number of auxiliary reagents for each reaction 
aux_products=[] #database containing the number of auxiliary products for each reaction
num_aux=0 #total number of auxiliary compounds
computed_reactions=[] #database of rdkit generated reactions (for an eventual calculation of free energy)
path=os.getcwd() #path for the current dataset (projectname)
opt_mol=True #if true make the MM optimization of the molecules saving them as .mol
gau_pre="%nprocshared=4"
gau_job="# opt freq b3lyp 6-31G(d)" #default gaussian job
orca_header=["! XTB2 opt numfreq","* XYZ "]
warm_start=False #if False starts from the beginning, if True it means we load a reaction already made and we start from that
debug=False #if debug print all the reaction products (see do_reaction)
printnewreaction=True
printreaction=False
printreactionnumber=False
limit=0
limitC=0
KEGG=False
KEGG_Title=''
warnings=[] #list of warning messages
EnergiesArePresent=False

def Load_SMARTS(filename): #open a text file containing SMART reactions and load them into memory
    filein = open(filename,'r') 
    lines = filein.readlines()
    for line in lines:
      try:
       line=line.split('#')[0]
      except:
       pass
      line=line.strip()
      try:
       field=line.split('/')
      except:
       pass
      else:
       if len(field)==3:   #add reaction only if all three fields are present
        Add_Reaction(field[0].strip(),field[1].strip(),field[2].strip())  
    filein.close()   

  
#Initial conditions defined here ----------------------------------------
def Set_Initial_Conditions():
 #for each reaction to be considered, insert the reaction in SMARTS and then call function Add_Reaction
 #when a reaction is added, auxiliary reagents and products could be indicated, respectively, as second and third parameter
    
 global pool,warm_start
 if (warm_start==False)or(warm_start and not(os.path.isfile(path+os.sep+"smarts_reactions.txt.old"))): 
  # if it is a cold start, or the previous job has no SMARTS reactions, load the default reactions from starting_reactions.txt
  Load_SMARTS('starting_reactions.txt')

 else: #warm start
  Load_SMARTS(path+os.sep+"smarts_reactions.txt.old") 

 if (warm_start==False): 
  #For the cold start, reagents are loaded from starting_reactants.txt. 
  #For the warm start, the already calculated reagents are loaded from the previous calculation
  filein = open('starting_reactants.txt','r') 
  lines = filein.readlines()
  for line in lines:
    try:
     line=line.split('#')[0].strip()
    except:
     pass
    if len(line)>=1:   #add reactant
      Write_Formula("",line,"","")   
  filein.close()     

#end of initial conditions ----------------------------------------------  

def Write_Formula(fname,smil,mol,m2):
 global opt_mol,compounds_file,path,energies,c_names,MWs,Formulas,reactant_is_enabled,warnings
 if smil=='': 
        warning="ERROR from Write_Formula. Smiles cannot be empty"
        print(warning)
        warnings.append(warning)
        1/0
 if mol=="":
  mol=Chem.MolFromSmiles(smil)
 if m2=="":
  m2=Chem.AddHs(mol)
 if fname=="": fname=len(pool)
 fname=str(fname)
 AllChem.Compute2DCoords(mol)
 Draw.MolToFile(mol,path+os.sep+fname+".png") # save reactant image
 MW=Descriptors.MolWt(mol)
 Formula=Chem.rdMolDescriptors.CalcMolFormula(mol)
 compounds_file.write(fname+","+Formula+","+smil+","+str(MW))
 c_names.append(fname)             # update databases
 Formulas.append(Formula)          # update databases
 pool.append(smil)                 # update databases
 MWs.append(MW)                    # update databases
 reactant_is_enabled.append(True)  # update databases  
 if (opt_mol):   
  try:  
   AllChem.EmbedMolecule(m2)
   AllChem.MMFFOptimizeMolecule(m2)
   props = AllChem.MMFFGetMoleculeProperties(m2) 
   mp = AllChem.MMFFGetMoleculeForceField(m2, props)
  except:
   warning="WARNING: Cannot optimize molecule "+fname
   print(warning) # Chem.AllChem.EmbedMultipleConfs(ligand, numConfs=10)
   warnings.append(warning)
  try:  # for homonuclear molecules CalcEnergy returns an error
   en=-mp.CalcEnergy() #for strange reasons we have to put the negative sign
  except:
   en=0     
  energies.append(en)
  compounds_file.write(","+str(en))
  try:
   molfile=Chem.MolToMolBlock(m2)
   file = open(path+os.sep+fname+".mol", "w")
   file.writelines(molfile)
   file.close()
  except: 
   warning="Error writing "+fname+".mol"
   print(warning)  
   warnings.append(warning)
 compounds_file.write("\n")
    
def Write_Reaction(SMARTS): #save images of the reaction
  global path,reactions_smarts
  try:
      reactions_smarts.append(SMARTS)
      num_reactions=len(reactions_smarts)-1
      rpath=path+os.sep+"REACTION"+str(num_reactions)
      if not os.path.exists(rpath):
        os.makedirs(rpath)  
      RP=SMARTS.split(">>")  
      Reactants=RP[0].split(".")
      Products=RP[1].split(".")  
      image_files=[]
      i=0
      for Reactant in Reactants:
       i+=1
       mol=Chem.MolFromSmiles(Reactant) 
       if mol is None:
        print("Unable to convert the SMILES string:"+Reactant)
       else: 
        filename=rpath+os.sep+"REACT"+str(i)+".png"
        image_files.append(filename)
        Draw.MolToFile(mol,filename)  
      i=0
      for Product in Products:
       i+=1
       mol=Chem.MolFromSmiles(Product)  
       if mol is None:
        print("Unable to convert the SMILES string:"+Product)
       else:
        filename=rpath+os.sep+"PROD"+str(i)+".png"
        image_files.append(filename)
        Draw.MolToFile(mol,filename) 
      images = [Image.open(im) for im in image_files]
      widths, heights = zip(*(i.size for i in images))
      total_width = sum(widths)+(len(image_files)-1)*20
      max_height = max(heights)
      new_im = Image.new('RGB', (total_width, max_height))
      draw = ImageDraw.Draw(new_im)
      draw.rectangle([(0,0),(total_width, max_height)],fill='white', outline='white')  #clear image
      new_font = ImageFont.truetype('arial.ttf', 40)    
      x_offset = 0
      c=0
      old=''  
      for im in images:
       new_im.paste(im, (x_offset,0))
       text=os.path.basename(image_files[c]).split(".",2)[0]
       if 'P' in text: new='P' 
       else: new='R'
       #draw.text((x_offset+im.size[0]/2, im.size[1]*3/4),text,fill='black')   
       conn='='   
       if old==new in text: conn='+'   
       old=new 
       if c>0: draw.text((x_offset-21, max_height/2-12),conn,fill='black',font=new_font)   
       x_offset += im.size[0]+20
       c+=1
      new_im.save(rpath+os.sep+"REACTION"+str(num_reactions)+".png")
  except:
      print("Error in write reaction")

def Add_Reaction(smarts,aux_reag,aux_prod):
 global reactions,num_reagents,aux_reagents,aux_products,num_aux,pool,reactions_file
 reactions_file.write(smarts+'/'+aux_reag+'/'+aux_prod+'\n')   
 Write_Reaction(smarts)
 t=smarts.split(">>") #try to understand the number of reactants involved
 num_reactants=t[0].count('.')+1 #num reactants = number of dots + 1 present in the reagents (before >>)
 r=rdChemReactions.ReactionFromSmarts(smarts)
 reactions.append(r)
 num_reagents.append(num_reactants)
 aux_name=""   
 if aux_reag:
     if not aux_reag in pool:
        aux_name="A"+str(num_aux)   
        Write_Formula(aux_name,aux_reag,"","")
        num_aux+=1
     else:
        aux_name=c_names[pool.index(aux_reag)] 
 aux_reagents.append(aux_name)       
 aux_name=""   
 if aux_prod:
     if not aux_prod in pool:
        aux_name="A"+str(num_aux)   
        Write_Formula(aux_name,aux_prod,"","")
        num_aux+=1       
     else:
        aux_name=c_names[pool.index(aux_prod)] 
 aux_products.append(aux_name)
 reaction_is_enabled.append(True)    

def Do_Reaction(reacts,rnum):
 global pool,reactions,computed_reactions,aux_reagents,aux_products,limitC
 reactants=[]
 rname=reactions[rnum]   
 for reagentnumber in reacts:
  reagent=Chem.MolFromSmiles(pool[reagentnumber])
  reactants.append(reagent)
 products = rname.RunReactants(reactants)
 reag="" #reaction reagent string
 for r in reacts:
  reag+=c_names[r]+"_" 
 if aux_reagents[rnum]:
  reag+=aux_reagents[rnum]+"_"
 reag+="R"+str(rnum)
 predicted_products_flat = [Chem.MolToSmiles(mol) for product_set in products for mol in product_set]
 output=[]
 for product in predicted_products_flat:
     mol=Chem.MolFromSmiles(product)
     try:
         m2=Chem.AddHs(mol) #gives an error for carbon extra valence
         if limitC>0:
          patt=Chem.MolFromSmarts("[C]") # get the number of carbon atoms
          if len(mol.GetSubstructMatches(patt))>limitC:
           if (debug): print("Exceeds the number of carbon atoms")
           continue
         output.append(product)
     except:
         if (debug): print("extra valence")
 for P in output: 
   rname=reag #creating the reaction string
   smil=P
   mol=Chem.MolFromSmiles(smil)
   m2=Chem.AddHs(mol) #gives an error for carbon extra valence
   if not(smil in pool): #add the compound if not already present in the database
     prodnum=str(len(pool))
     if (debug): print(len(prodnum),smil)   
     Write_Formula(prodnum,smil,mol,m2)
   else:
     prodnum=c_names[pool.index(smil)]
     if (debug): print(prodnum,smil,pool[int(prodnum)])
   rname+="_"+prodnum
   if aux_products[rnum]:
    rname+="_"+aux_products[rnum]
   if (debug): print(rname)   
   if not(rname in computed_reactions): #add the reaction if not already present in the database
     computed_reactions.append(rname)
     if(printnewreaction): print(rname)
    
def Calc_Frequencies():
 global MWs,path
 tolerance=0.01
 MW_list=[[]]
 del MW_list[0]   #unknowingly the first element is empty when the list is created
 for MW in MWs:
  position=0
  found=False  
  for x in MW_list:
   if len(x)>0: 
    if (float(x[0])+tolerance>=MW)and(float(x[0])-tolerance<=MW):
     MW_list[position][1]+=1
     found=True   
     break
   position+=1  
  if (found==False):
   MW_list.append([MW,1])
 frequency_file=open(path+os.sep+'FREQUENCY.csv','w') 
 frequency_file.write("MW,number of occurrences\n")
 MW_list = sorted(MW_list, key = lambda x: x[0])   
 for x in MW_list:
         frequency_file.write(str(x[0])+","+str(x[1])+"\n")
 frequency_file.close()    

def Conv_mol2input(type):
    global gau_job,gau_pre,orca_header,path
    if type=="gaussian":
        fileext="gjf"
    elif type=="orca":
        fileext="inp"
    elif type=="xyz":
        fileext="xyz"        
    else:
        print("Error Conv_mol2input: type must be gaussian, orca or xyz")
    mol_list = glob.glob(path+os.sep+'*.mol')
    for file in mol_list:
        filein = open(file, 'r')
        Lines = filein.readlines()
        Charge=0
        for s in Lines:
            if 'M  CHG' in s:
              sp=s.split()
              for i in range(int(sp[2])):
                    Charge+=int(sp[4+i])
        numatoms=int(Lines[3][:3])
        outfile = open(file[:-3]+fileext, "w")        
        if type=="gaussian":
         outfile.write(gau_pre+"\n"+gau_job+"\n \n"+"CheReNetw\n"+" \n"+str(Charge)+" 1\n")
        if type=="orca":
         outfile.write(orca_header[0]+"\n"+orca_header[1]+str(Charge)+" 1\n")
        if type=="xyz":
         outfile.write(str(numatoms)+"\nCreated by AutoChem\n")         
        for i in range(numatoms):
            tmp=Lines[4+i].split()
            outfile.write(tmp[3]+" "+tmp[0]+" "+tmp[1]+" "+tmp[2]+"\n")
        if type=="gaussian":
         outfile.write("\n\n")
        if type=="orca":
         outfile.write("*\n")
        if type=="xyz":
         outfile.write("\n")         
        outfile.close()

def getURL(url):
 fp = urllib.request.urlopen(url)
 mybytes = fp.read()
 mystr = mybytes.decode("utf8")
 fp.close()
 return(mystr)

def getKEGGpathway(pathwayname):
 mystr=getURL("https://www.kegg.jp/pathway/"+pathwayname)
 print (len(mystr),' bytes read')
 modules=[]
 pos=0
 while not(mystr.find('/module/M',pos)==-1):
    pos=mystr.find('/module/M',pos)
    module=mystr[pos+8:mystr.find('"',pos)]
    modules.append(module)
    pos+=5
 print('Found',len(modules),'modules')
 for m in modules:
     getKEGGmodule(m)

def getKEGGmodule(modulename):
 global computed_reactions,c_names,KEGG_Title   
 mystr=getURL('https://www.kegg.jp/module/'+modulename)
 print (len(mystr),' bytes read')
 titlestart=mystr.find('<td>Name</td>')
 titlestart=mystr.find('<td>',titlestart+4)
 titlestop=mystr.find('</td>',titlestart)
 title=mystr[titlestart+4:titlestop]
 title=title.replace('&gt;','>')
 print(title)
 KEGG_Title+=title+'\n'
 start=mystr.find('<td>Reaction</td>')
 stop=mystr.find('</tr>',start)
 reactions=mystr[start:stop]
 ReactionsURLToSearch=[]
 pos=0
 while not(reactions.find('/entry/R',pos)==-1):
    pos=reactions.find('/entry/R',pos)
    reaction=reactions[pos:reactions.find('">',pos)]
    ReactionsURLToSearch.append('https://www.kegg.jp'+reaction)
    pos+=5
 print('\n',len(ReactionsURLToSearch),' reactions found')
 for URL in ReactionsURLToSearch:
    string=getURL(URL)
    reactionName=URL[URL.find('/R')+1:]
    start=string.find('Equation')
    stop=string.find('</tr>',start)
    string=string[start:stop]
    pos=0
    reactionstechio=''
    while not(string.find('/entry/C',pos)==-1):
     pos=string.find('/entry/C',pos)
     pos2=string.find('">',pos)
     reactant=string[pos+7:pos2]    
     if  not(reactant in c_names):
         c_names.append(reactant)
     pos3=string.find('</a>',pos)
     reactionstechio+=string[pos2+2:pos3]
     pos4=string.find('<a',pos)
     if pos4>0: reactionstechio+=string[pos3+4:pos4]
     pos+=5
    reactionstechio=reactionstechio.replace('&lt;=&gt;','_'+reactionName+'_')
    reactionstechio=reactionstechio.replace('+','_')
    pos=reactionstechio.find(' 2 ')
    while (pos>0):
     pos2=reactionstechio.find(' ',pos+5)
     if pos2==-1: pos2=len(reactionstechio)
     reac=reactionstechio[pos+3:pos2]
     reactionstechio=reactionstechio.replace(' 2 ',reac+'_',1)
     pos=reactionstechio.find(' 2 ',pos+5)
    reactionstechio=reactionstechio.replace(' ','') 
    if reactionstechio not in computed_reactions:
     computed_reactions.append(reactionstechio)
     print(reactionstechio) # print only reactions not present in database

def ReactionEditor():
    global reactant_is_enabled, reaction_is_enabled, pool, reactions_smarts,c_names
    reactant_is_enabled=[True]*len(pool)
    #reaction_is_enabled=[True]*len(reactions_smarts)
    MaxItems=20
    Editing=True
    while Editing:
     print('')
     print(str(len(pool))+' Reactants loaded: ')   
     if len(pool)>MaxItems:
      print(' (only the first '+str(MaxItems)+' will be printed)')   
      max=MaxItems
     else:
      max=len(pool)   
     for i in range(0,max):
           aux='Auxiliary compound' if 'A' in c_names[i] else ''
           disabled = "" if reactant_is_enabled[i] else " ---- DISABLED ----" 
           print(i,pool[i],Formulas[i],aux,disabled)
     print('')     
     print(str(len(reactions_smarts))+' Reactions loaded: ')   
     if len(reactions_smarts)>MaxItems:
      print(' (only the first '+str(reactions_smarts)+' will be printed)')   
      max=MaxItems
     else:
      max=len(reactions_smarts)   
     for i in range(0,max):
           disabled = "" if reaction_is_enabled[i] else " ---- DISABLED ----" 
           print(i,reactions_smarts[i],disabled)
     ask=input("\nEdit Reactants or rEactions? [r/e], ENTER to accept ").upper()
     if ask=="R":
         print("   editing reactants")
         ask=input("Action: Disable/enable or Add? [d/a] ").upper()
         if ask=="A":
             molecule=input("Insert a molecule in SMILES format: ")
             #add reactant
             Write_Formula("",molecule,"","")
             manually_added_reactants.append(molecule)
         if ask=="D":    
           try:   
            num=int(input("Reactant number? [0-"+str(len(pool)-1)+"] "))
           except:
            num=-1  
           if num<0 or num>len(pool)-1:
              print("  ------- insert a number in the range")
           else:
              #disable/enable reactant
              reactant_is_enabled[num]=not(reactant_is_enabled[num])               
     elif ask=="E":
         print("   editing reactions")
         ask=input("Action: Disable/enable or Add? [d/a] ").upper()
         if ask=="A":
             smarts=input("Insert a valid reaction in SMARTS format: ")
             aux_r=input("Insert (if there is) an aux reagent in SMILES format: ")
             aux_p=input("Insert (if there is) an aux product in SMILES format: ")
             #add reaction
             Add_Reaction(smarts,aux_r,aux_p)
         else:    
           try:   
            num=int(input("Reaction number? [0-"+str(len(reactions_smarts)-1)+"] "))
           except:
            num=-1   
           if num<0 or num>len(reactions_smarts)-1:
              print("  ------- insert a number in the range")
           else:
              #disable/enable reaction
              reaction_is_enabled[num]=not(reaction_is_enabled[num])
     elif ask=="":
         Editing=False
     else:
         print(" ------- please type R,E or enter")

    

# --------------------------------------------------------------------------------------------------#
#                                                                                                   #
#                                    M A I N       P R O G R A M                                    #
#                                                                                                   #
# --------------------------------------------------------------------------------------------------#
print("    ___         __        ________                      ___    ____  ")
print("   /   | __  __/ /_____  / ____/ /_  ___  ____ ___     |__ \  / __ \ ")
print("  / /| |/ / / / __/ __ \/ /   / __ \/ _ \/ __ `__ \    __/ / / / / / ")
print(" / ___ / /_/ / /_/ /_/ / /___/ / / /  __/ / / / / /   / __/_/ /_/ /  ")
print("/_/  |_\__,_/\__/\____/\____/_/ /_/\___/_/ /_/ /_/   /____(_)____/   ")
print()
warnings=[] #clear all warning messages
projectname=input("Insert project name (files will be stored in a directory with the project name) ")
if projectname=="":
    print("Error. project name cannot be blank")
    1/0
forbidden=["?",'"',":","/","*","|","<",">"]
if any([char in projectname for char in forbidden]):
    print("Error. project contains illegal characters")
    1/0    
if not os.path.exists(projectname):
    os.makedirs(projectname)  
path=projectname
to_be_removed=['*.mol','*.png','*.gjf','*.inp','REACTIONS.txt','smarts_reactions.txt','COMPOUNDS.csv','FREQUENCY.csv']
file_lists=[]
dir_list=[]
items = os.listdir(path)
for item in items: #search directories to be removed (aka: graphical representation of reactions)
    if os.path.isdir(os.path.join(path,item)) and 'REACTION' in item:
        dir_list.append(item)
for file_extension in to_be_removed:
    filelist=glob.glob(path+os.sep+file_extension)
    if len(filelist)>0:
     file_lists.append(filelist)
if len(file_lists)+len(dir_list)>0:
    ask=input("Output files already present. \nPress Y to delete, any other keys to keep (common files will still be overwritten): ").upper()
    if ask=="Y":
     totfiles=0
     for l in file_lists:
        for f in l:
            totfiles+=1
     tmp=""       
     if len(dir_list)>0: tmp=" AND "+str(len(dir_list))+" DIRECTORIES"       
     ask=input("ARE YOU SURE DO YOU WANT TO DELETE "+str(totfiles)+" FILES"+tmp+"?\nPress Y to confirm ").upper()
     if ask=="Y":
      for file_list in file_lists:
        for file in file_list:
            try:
                os.remove(file)
            except:
                print("Error deleting the file",file)
      for directory in dir_list:
            try:
                shutil.rmtree(os.path.join(path,directory))
            except:
                print("Error deleting directory ",directory)
                
    else:
      if (os.path.exists(path+os.sep+"COMPOUNDS.csv")) and (os.path.exists(path+os.sep+"REACTIONS.txt")):
        ask2=input("Load the compounds already calculated? [Y/n]\nIf reply y or hit ENTER, the calculation will start with the already present molecules ").upper()
        if ask2=='Y' or ask2=='':  #                load the compounds already calculated
            with open(path+os.sep+"COMPOUNDS.csv", "r") as fp:   # read compounds
               lines = fp.readlines()  
            fp.close()   
            if len(lines[1].split(','))>4: EnergiesArePresent=True  # if there are more than 4 fields we have energies too
            else: EnergiesArePresent=False   
            lines=lines[1:] # remove header
            for line in lines:
                sp=line.split(',')
                c_names.append(str(sp[0]))
                Formulas.append(str(sp[1]))
                pool.append(str(sp[2]))
                MWs.append(float(sp[3]))
                if EnergiesArePresent: energies.append(float(sp[4]))
            with open(path+os.sep+"REACTIONS.txt", "r") as fp:   # read reactions
               computed_reactions=[] #initialize array
               lines = fp.readlines()
               for line in lines:
                 try:
                   line=line.split(',')[0]
                 except:
                   pass                   
                 line=line.rstrip()  #remove CR LF
                 computed_reactions.append(line)  
            fp.close()     
            warm_start=True

ask=input("Do you want to download KEGG reactions or start a calculation? [k/C] ").upper()
if ask=="K":
    KEGG=True
    KEGG_Type=input("KEGG Module or Pathway? [m/p] ").upper()
    if KEGG_Type=="M": 
        KEGG_Type="module"
        example="M00004"
    elif KEGG_Type=="P": 
        KEGG_Type="pathway"
        example="map00010"
    else:
      print("Please choose M or P")  
      1/0
    KEGG_Name=input("Please insert the "+KEGG_Type+" name without spaces [ex.: "+example+"] ")    
    if KEGG_Name=='':
        print("Insert a valid name")
        1/0

if (warm_start==True):
    mode='a'
else:
    mode='w'
compounds_file=open(path+os.sep+"COMPOUNDS.csv",mode) #a file created with the corrispondence between the compound number and its molecular formula and MW
if (EnergiesArePresent==False):
 ask_mol = input("Optimize molecules? [Y/n]: ").upper()
 if ask_mol=="N":
    opt_mol=False
 else:
    opt_mol=True
else: opt_mol=True # if energies were present previously, force the optimization    
if opt_mol:
    if (warm_start) and not(EnergiesArePresent):
     print(" We have to calculate energies before proceeding\n it might take some minutes")
     energies=[]
     for i in range(0,len(pool)):
        smil=pool[i]
        mol=Chem.MolFromSmiles(smil)
        m2=Chem.AddHs(mol)
        AllChem.EmbedMolecule(m2)
        AllChem.MMFFOptimizeMolecule(m2)
        props = AllChem.MMFFGetMoleculeProperties(m2) 
        mp = AllChem.MMFFGetMoleculeForceField(m2, props)
        try:  # for homonuclear molecules CalcEnergy returns an error
         en=-mp.CalcEnergy() #for strange reasons we have to put the negative sign
        except:
         en=0     
        energies.append(en)
        molfile=Chem.MolToMolBlock(m2) # write mol file
        file = open(path+os.sep+c_names[i]+".mol", "w")
        file.writelines(molfile)
        file.close()        
     compounds_file.close()    
     compounds_file=open(path+os.sep+"COMPOUNDS.csv","w") # we have to rewrite the whole file because some reagents haven't energies
     compounds_file.write("name,formula,smiles,MW,Energy\n")
     for i in range(0,len(pool)):
      compounds_file.write(c_names[i]+','+Formulas[i]+','+pool[i]+','+str(MWs[i])+','+str(energies[i])+"\n")  
    
    ask_gau=input('Do you want to create Gaussian input file? [Y/n] ').upper()
    if ask_gau=='N':
     gau=False
    else:
     gau=True
     ask=str(input('Gaussian precommand: ('+gau_pre+'): '))
     if ask!="":
        gau_pre=ask    
     ask=str(input('Gaussian command: ('+gau_job+'): '))
     if ask!="":
        gau_job=ask
    ask_orca=input('Do you want to create ORCA input file? [Y/n] ').upper()
    if ask_orca=='N':
     orca=False
    else:
     orca=True
     ask=str(input('Orca command: ('+orca_header[0]+'): '))
     if ask!="":
        orca_header[0]=ask    
    ask_xyz=input('Do you want to create xyz (CREST, xTB) input file? [Y/n] ').upper()
    if ask_xyz=='N':
     xyz=False
    else:
     xyz=True
        
else:
    gau=False
    orca=False
    xyz=False
    
if (warm_start==False): 
    compounds_file.write("name,formula,smiles,MW")
    if (opt_mol):
     compounds_file.write(",Energy")
    compounds_file.write("\n")       

if KEGG==False:    # if not KEGG we are going to use the reactions written in smarts format
 rfilename=path+os.sep+"smarts_reactions.txt"   
 if os.path.isfile(rfilename): #if the file exists rename in .old
        if os.path.isfile(rfilename+'.old'): #if the db.old file exists, remove it
            os.remove(rfilename+'.old')
        if os.stat(rfilename).st_size > 0:    
         os.rename(rfilename, rfilename+'.old')
 reactions_file=open(rfilename,'w') #a file containing the reactions to be computed      
 Set_Initial_Conditions()    
 print()   
 ReactionEditor()
 reactions_file.close()   
 while True:   
  try:   
   cycles = int(input('How many cycles to calculate: '))
   if cycles<=0: cycles=1/0        
   break   
  except:
   print('Insert a valid value')   
 ask_l = input("Limit the number of reactants for bimolecular reactions? (insert the number or leave blank for unlimited): ")
 try:
  limit=int(ask_l)
 except:    
  limit=0
 ask_l = input("Limit the number of Carbon atoms for products? (insert the number or leave blank for unlimited): ")
 try:
  limitC=int(ask_l)
 except:    
  limitC=0
     
if not(limit==0): print("\nReactions are limited to the first ",limit," reagents in pool\n")
if not(limitC==0): print("\nReactions are limited to the products having max ",limitC," carbon atoms\n")
if KEGG: print("\nReactions are taken from KEGG database "+KEGG_Type+" "+KEGG_Name+"\n")
start = time.time()
if KEGG==False:
 # perform reactions    perform reactions    perform reactions    perform reactions    perform reactions    
 for c in range(cycles):
  print("Cycle ",c+1," of ",cycles," cycles")
  l=len(pool)
  for reaction in range(len(reactions)):
   if reaction_is_enabled[reaction]:
    if (printreactionnumber): print('Reaction ',reaction+1,':')      
    if num_reagents[reaction]==1:  #reaction having one reagent
     for x in range(0,l,1):
      if reactant_is_enabled[x]:   
       if (printreaction): print(pool[x]," =")
       Do_Reaction([x],reaction)
    elif num_reagents[reaction]==2: #reaction having two reagents
     last=l    
     if (not(limit)==0)and(l>limit): 
          last=limit
     for x in range(0,last,1):
      for y in range(x,l,1):
       if reactant_is_enabled[x] and reactant_is_enabled[y]:   
        if (printreaction): print(pool[x]," + ",pool[y]," =")
        Do_Reaction([x,y],reaction)  #swap reagents and call reaction 2 times
        Do_Reaction([y,x],reaction)
    elif num_reagents[reaction]==3: #reaction having three reagents, never tested
     last=l    
     if (not(limit)==0)and(l>limit): 
          last=limit
     for x in range(0,last,1):
      for y in range(x,l,1):
       for z in range(y,l,1):          
        if reactant_is_enabled[x] and reactant_is_enabled[y] and reactant_is_enabled[z]:   
         if (printreaction): print(pool[x]," + ",pool[y]," + ",pool[z]," =")
         Do_Reaction([x,y,z],reaction)  
         Do_Reaction([x,z,y],reaction)  
         Do_Reaction([y,x,z],reaction)  
         Do_Reaction([y,z,x],reaction)  
         Do_Reaction([z,x,y],reaction)  
         Do_Reaction([z,y,x],reaction)  
    else:
     print('ERROR: number of reagents not covered. Current version deals up to trimolecular reactions')
else:
 # KEGG database   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG   KEGG
 len_names=len(c_names)  #get the length of c_names array before KEGG (in the case of previous jobs)
 if KEGG_Type=='module': getKEGGmodule(KEGG_Name)
 else: getKEGGpathway(KEGG_Name)
 newR=len(c_names)-len_names   #we will add only newly found compounds
 print('\n',newR,' reactants found')
 does_we_have_protons='' #this variable will be set to the compound name if we have protons (incomputable), that will be substituted automatically with hydronium ion and water
 newR_names=c_names[len_names:]
 c_names=c_names[:len_names] # we removed the newly added compounds from c_names: some compounds might not be real molecules in KEGG
 for compound in newR_names: # download new reactants
    URL='https://www.genome.jp/entry/-f+m+'+compound
    string=getURL(URL)
    string=string.replace('R ','H ') # substitute general R with H or C (CH3)
    fname=URL[URL.find('m+')+2:]+'.mol'
    if len(string)>0:
      print(fname)
      print(string,  file=open(path+os.sep+fname, 'w'))
      m=Chem.rdmolfiles.MolFromMolFile(path+os.sep+fname)
      smil=Chem.rdmolfiles.MolToSmiles(m) 
      if smil=='[H+]':
          if not(does_we_have_protons==''):
              tmp='Gosh, we have more protons with different file names. Program will generate messy output'
              print(tmp)
              warnings.append(tmp)
          warnings.append('H+ was converted automatically to H3O+')              
          does_we_have_protons=fname[:fname.find('.')] # the fname is a proton, we have to change into hydronium ion
          smil='[OH3+]'
      Write_Formula(fname[:fname.find('.')],smil,"","")        
    else:
        warningtext='WARNING   Some problem encountered with '+fname+' file size is 0.\n Please check on KEGG website the reason\n https://www.genome.jp/entry/'+fname[:fname.find('.')]
        print(warningtext)
        warnings.append(warningtext)

 if not(does_we_have_protons==''): #if we have protons, that were converted to hydronium ions, we must add water from the other side of the reaction to mainain the balance
    try:
        water=c_names[pool.index('O')]
    except:
        #we haven't water molecules, so we that to add it
        print('Added water to pool')
        Write_Formula('', 'O','','')
        water=c_names[len(pool)]
    for reaction in range(len(computed_reactions)):
        rp=computed_reactions[reaction].split('_R')# _R is the delimeter between reactants and products
        r=rp[0]
        p=rp[1]
        n_H_in_r=r.count(does_we_have_protons)
        n_H_in_p=p.count(does_we_have_protons)
        if (n_H_in_p+n_H_in_r)>0:        
         for i in range(n_H_in_r):
            p+='_'+water # for each H3O+ we add water from the other side
         for i in range(n_H_in_p):
            r+='_'+water
         warnings.append('water was added to reaction R'+computed_reactions[reaction].split('_R')[1].split('_')[0]+' to allow the formation of H3O+')   
         computed_reactions[reaction]=r+'_R'+p

reaction_file=open(path+os.sep+'REACTIONS.txt','w') #create a reaction files containing all the reactions to be calculated
for element in computed_reactions:
 reaction_file.write(element)
 if (opt_mol): #calculates deltaG with rdkit energy
    ele=element.split("_")
    sign=-1
    deltaG=0
    for rp in ele:
     if "R" in rp: sign=-sign
     else:
      try:  
       deltaG+=energies[c_names.index(rp)]*sign
      except:
       deltaG='some compounds are missing'     
       break
    reaction_file.write(","+str(deltaG))
 reaction_file.write('\n')   
reaction_file.close()
compounds_file.close()    

Calc_Frequencies() #calculate a distribution curve of compounds according to their molecular weight
if gau: 
    Conv_mol2input("gaussian")
if orca: 
    Conv_mol2input("orca")  
if xyz: 
    Conv_mol2input("xyz")      
end = time.time()
time_string='\nTOTAL NUMBER OF COMPOUNDS: '+str(len(pool))+' involved in '+str(len(computed_reactions))+' REACTIONS\n calculated in '+str(round(end-start,2))+' seconds\n\n'
print(time_string)
#                                        save a summary html file with all the reactions
reaction_html=open(path+os.sep+'reactions.html','w')
reaction_html.write('<html><body>\n')
if len(reactions_smarts)>0:
  for r in range(len(reactions_smarts)):
    r_image="REACTION"+str(r)+os.sep+"REACTION"+str(r)+".png"    
    reaction_html.write('<h1>Reaction type '+str(r)+'</h1><br>')
    reaction_html.write('<img src="'+r_image+'"><br>')    
    reaction_html.write('<br>\n')
  reaction_html.write('<hr>\n')      
reaction_html.write('<table>')
count=0
for element in computed_reactions:
  reaction_html.write('<tr><td><h1>'+str(count)+')</h1></td>')  
  count+=1  
  ele=element.split("_")
  first=True  
  for rp in ele:
    if "R" in rp:
      reaction_html.write('<td><h1>=</h1></td>')
      first=True  
      reaction_type=rp  
    else:    
      if not(first): reaction_html.write('<td><h1>+</h1></td>')  
      reaction_html.write('<td><img src=\"'+rp+'.png\"></td>')  
      first=False    
  reaction_html.write('<td><h1>('+reaction_type+')</h1></td></tr>\n')      
reaction_html.write('</table>\n</body></html>')
reaction_html.close()
#save summary to a file
summary_file=open(path+os.sep+'SUMMARY.txt',mode)
summary_file.write('Starting time: '+str(datetime.datetime.now())+'\n')
if KEGG:
 summary_file.write("Type of job: KEGG download "+KEGG_Type+" "+KEGG_Name+'\n') 
 summary_file.write('KEGG names of module(s) explored: \n'+KEGG_Title)
else:    
 summary_file.write("Type of job: rdkit calculation\n")   
 summary_file.write("Reactions considered in this calculation (SMARTS format):\n")
 count=0   
 for reaction in reactions_smarts:
    summary_file.write(str(count)+') ') 
    count+=1 
    summary_file.write(reaction+"\n")
 if not(limit==0): summary_file.write("\nReactions are limited to the first "+str(limit)+" reagents in pool\n")
 if not(limitC==0): summary_file.write("\nReactions are limited to the products having max "+str(limitC)+" carbon atoms\n")  
 for r in range(len(reactions)):    
     if not(reaction_is_enabled[r]): summary_file.write('Reaction '+str(r)+' DISABLED\n')
 for r in range(len(reactant_is_enabled)):    
     if not(reactant_is_enabled[r]): summary_file.write('Reactant '+str(r)+' ('+Formulas[r]+') DISABLED\n')
 if len(manually_added_reactants)>0:
     summary_file.write('Manually added reactants:\n')
     for r in manually_added_reactants:
         summary_file.write(r+'\n')

 summary_file.write("\nNumber of cycles calculated: "+str(cycles)+"\n")
if len(warnings)>0:
     summary_file.write('\n    Warning messages:\n')     
     for w in warnings:
         summary_file.write(w+'\n') 
     summary_file.write('\n')     
summary_file.write(time_string)
summary_file.close()
