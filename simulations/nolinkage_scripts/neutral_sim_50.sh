#!/bin/bash
#SBATCH --mail-user=adaigle@unc.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1 #number of tasks
#SBATCH --time=1-06:00

module load perl/5.36.0
echo "SLURM_JOBID: " $SLURM_JOBID

######################################
#To be run on AGAVE!
#10 simulations per submitted job!
#######################################


declare -i repID=0+$SLURM_JOBID

module load slim/4.0.1

cd /nas/longleaf/home/adaigle/simulate_DFEs/nolinkage_scripts
echo "starting simulation " $repID
folder="/nas/longleaf/home/adaigle/work/johri_elegans/sim_outputs/nolinkage/neutral_outputs"

#eqm, no selfing
#slim -d d_Ncur=5000 -d d_selfing_prob=0 -d d_beta=0.5 -d d_gamma_del=0 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/eqm_selfing0'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 99%
#slim -d d_Ncur=5000 -d d_selfing_prob=0.99 -d d_beta=0.5 -d d_gamma_del=0 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/eqm_selfing99'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 50%, DFE1
slim -d d_Ncur=5000 -d d_selfing_prob=0.5 -d d_beta=0.5 -d d_gamma_del=0 -d d_f_pos=0.0 -d d_gamma_pos=0 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='${folder}/eqm_selfing50'" demo_DFE_gamma_selfing_Celegans.slim



echo "Finished simulation " $repID





