#!/bin/bash
#SBATCH --mail-user=adaigle@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=11-00:00
#SBATCH -a 1-5%5
#SBATCH -o /nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/LOGFILES/dfeselfing_%A_rep%a.out
#SBATCH -e /nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/LOGFILES/dfeselfing_%A_rep%a.err
#SBATCH --mem=16g 

module load perl/5.36.0
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
 
declare -i repID=0+$SLURM_ARRAY_TASK_ID

module load slim/4.0.1

cd /nas/longleaf/home/adaigle/DFESelfing/simulations/scaling_scripts
echo "starting simulation " $repID
folder="/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/scaling/"
                                                        
slim -d d_Ncur=5000 -d d_selfing_prob=0.95 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/eqm_selfing95_DFE2'" demo_DFE_gamma_selfing_Celegans.slim
 
echo "Finished simulation " $repID





