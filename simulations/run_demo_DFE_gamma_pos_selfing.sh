#!/bin/bash
#SBATCH --mail-user=pjohri1@asu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 2 #number of tasks
#SBATCH --time=2-06:00
#SBATCH -a 1-1%1
#SBATCH -o /home/pjohri1/LOGFILES/dfeselfingpos_%A_rep%a.out
#SBATCH -e /home/pjohri1/LOGFILES/dfeselfingpos_%A_rep%a.err

module load perl/5.22.1
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
######################################
#To be run on AGAVE!
#10 simulations per submitted job!
#######################################


declare -i repID=0+$SLURM_ARRAY_TASK_ID

module load slim/3.7

cd /home/pjohri1/DFESelfing/simulations
echo "starting simulation " $repID

#eqm, no selfing, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing0_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, no selfing, DFE2                                                          
#slim -d d_Ncur=5000 -d d_selfing_prob=0 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing0_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, no selfing, DFE3                                                          
slim -d d_Ncur=5000 -d d_selfing_prob=0 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing0_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 99%, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0.99 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing99_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 99%, DFE2
#slim -d d_Ncur=5000 -d d_selfing_prob=0.99 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing99_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 99%, DFE3
#slim -d d_Ncur=5000 -d d_selfing_prob=0.99 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing99_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 50%, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0.5 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing50_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 50%, DFE2
#slim -d d_Ncur=5000 -d d_selfing_prob=0.5 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing50_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 50%, DFE3
#slim -d d_Ncur=5000 -d d_selfing_prob=0.5 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing50_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 80%, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0.8 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d d_dom_del=0.5 -d d_dom_pos=0.5 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing80_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 80%, DFE2
#slim -d d_Ncur=5000 -d d_selfing_prob=0.8 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing80_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 80%, DFE3
#slim -d d_Ncur=5000 -d d_selfing_prob=0.8 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing80_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 90%, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0.9 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing90_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 90%, DFE2
#slim -d d_Ncur=5000 -d d_selfing_prob=0.9 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing90_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 90%, DFE3
#slim -d d_Ncur=5000 -d d_selfing_prob=0.9 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing90_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 95%, DFE1
#slim -d d_Ncur=5000 -d d_selfing_prob=0.95 -d d_beta=0.9 -d d_gamma_del=5 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing95_DFE1_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 95%, DFE2
#slim -d d_Ncur=5000 -d d_selfing_prob=0.95 -d d_beta=0.5 -d d_gamma_del=50 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing95_DFE2_pos'" demo_DFE_gamma_selfing_Celegans.slim

#eqm, selfing rate 95%, DFE3
#slim -d d_Ncur=5000 -d d_selfing_prob=0.95 -d d_beta=0.3 -d d_gamma_del=1000 -d d_f_pos=0.01 -d d_gamma_pos=200 -d "d_repID='${repID}'" -d "d_folder='/scratch/pjohri1/DFESelfing/gammaDFE/eqm_selfing95_DFE3_pos'" demo_DFE_gamma_selfing_Celegans.slim


echo "Finished simulation " $repID





