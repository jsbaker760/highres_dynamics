#!/bin/bash
#SBATCH -p sched_mem1TB_centos7,sched_mem1TB,defq,sched_mem4TB
#SBATCH --time=1:00:00
#SBATCH -o err1.txt
#SBATCH -e err1.txt
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:3
#SBATCH --cpus-per-task=3

mkdir -p ./$SLURM_JOB_ID
mkdir -p ./MATLAB_PREFDIR

module add mit/matlab/2016b
export MATLAB_PREFDIR='/scratch/mit_lieberman/projects/jsb_cuti/all_samples/dm/MATLAB_PREFDIR'
ulimit -u 63536 
matlab -logfile distance_matrix.log -nodisplay -nodesktop -nosplash -softwareopengl  < get_distance_matrices_GPU.m

rm -rf ./$SLURM_JOB_ID
