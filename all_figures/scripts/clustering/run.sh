#!/bin/bash

mkdir -p ./MATLAB_PREFDIR

export MATLAB_PREFDIR='/scratch/mit_lieberman/projects/jsb_cuti/all_samples/clustering/sweep/MATLAB_PREFDIR'

module add mit/matlab/2018a

sbatch -p sched_mem1TB_centos7,sched_mem1TB,defq,sched_mem4TB \
--time=5-00:00:00 \
-o sout.txt \
-e srr1.txt \
--mem=300G \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=12 \
--wrap=" mkdir -p ./$SLURM_JOB_ID 

ulimit -u 63536 
matlab -nodisplay -nodesktop -nosplash -softwareopengl  < \
SepidermidisATCC12228_ParameterSweep.m"

sbatch -p sched_mem1TB_centos7,sched_mem1TB,defq,sched_mem4TB \
--time=5-00:00:00 \
-o cout.txt \
-e crr.txt \
--mem=300G \
--nodes=1 \
--ntasks=1 \
--cpus-per-task=12 \
--wrap=" mkdir -p ./$SLURM_JOB_ID 

ulimit -u 63536 
matlab -nodisplay -nodesktop -nosplash -softwareopengl  < \
Pacnes_C1_ParameterSweep.m"