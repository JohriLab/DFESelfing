#!/bin/bash
#SBATCH --mail-user=adaigle@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -n 2 #number of tasks
#SBATCH --time=1-06:00
#SBATCH -a 1-10000%20
#SBATCH -o /nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/LOGFILES/dfeselfing_%A_rep%a.out
#SBATCH -e /nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/LOGFILES/dfeselfing_%A_rep%a.err

module load perl/5.36.0
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_IDmodule load perl/5.36.0
module load slim/4.0.1
cd /nas/longleaf/home/adaigle/simulate_DFEs/nolinkage_scripts
folder="/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/nolinkage"

declare -i repID=0+$SLURM_ARRAY_TASK_ID
#slim -d d_Ncur=5000 -d d_selfing_prob=0 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/eqm_selfing0_DFE1_nolinkageneu'" demo_DFE_gamma_selfing_Celegans_unlinked_neu.slim
#slim -d d_Ncur=5000 -d d_selfing_prob=0.5 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/DFE1_50/eqm_selfing50_DFE1_nolinkageneu'" demo_DFE_gamma_selfing_Celegans_unlinked_neu.slim
slim -d d_Ncur=5000 -d d_selfing_prob=0.99 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/DFE1_99/eqm_selfing99_DFE1_nolinkageneu'" demo_DFE_gamma_selfing_Celegans_unlinked_neu.slim
