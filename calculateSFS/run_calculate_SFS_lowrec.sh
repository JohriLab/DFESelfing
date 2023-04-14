#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=0-00:10
#SBATCH -o /home/pjohri1/LOGFILES/sfs_selfing_%j.out
#SBATCH -e /home/pjohri1/LOGFILES/sfs_selfing_%j.err

module load perl/5.22.1

######################################
#To be run on AGAVE!
#1 job per submitted job!
#######################################



folder="/scratch/pjohri1/DFESelfing/gammaDFE_lowrec"
subfolder="eqm_lowrec001_DFE3"

#unzip folder:                                                                  
cd ${folder}                                                                    
unzip ${subfolder}.zip                                                          
                                                                                
cd /home/pjohri1/SlimStats                                                      
echo "getting SFS"

#neutral
python get_sfs_from_full_output_general_reps.py -inputFolder ${folder}/${subfolder} -inputFileExt ".txt" -outputFolder ${folder}/SFS -outputPrefix ${subfolder}_neu -mutnTypes m1 -numGenomes 100

#selected
python get_sfs_from_full_output_general_reps.py -inputFolder ${folder}/${subfolder} -inputFileExt ".txt" -outputFolder ${folder}/SFS -outputPrefix ${subfolder}_sel -mutnTypes m2 -numGenomes 100

cd /home/pjohri1/DFESelfing/calculateSFS                                        
echo "getting the fixed class"

#neutral fixed class 
python get_fixed_class.py -inputFolder ${folder}/${subfolder} -inputFileExt ".fixed" -outputFolder ${folder}/SFS -outputPrefix ${subfolder}_neu -mutnTypes m1 -numGenomes 100 -NAnc 5000
                                                                                
#selected fixed class                                                           
python get_fixed_class.py -inputFolder ${folder}/${subfolder} -inputFileExt ".fixed" -outputFolder ${folder}/SFS -outputPrefix ${subfolder}_sel -mutnTypes m2 -numGenomes 100 -NAnc 5000

#remove folder                                                                  
cd ${folder}                                                                    
rm -r ${subfolder}

echo "Finished getting the SFS"




