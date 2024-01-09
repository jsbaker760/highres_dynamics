  #!/bin/bash

# ntasks
SM_PARAMS="job-name ntasks partition time mail-user mail-type error output"

SM_ARGS="  -x node083 --no-requeue --parsable --cpus-per-task {cluster.cpus-per-task} --mem-per-cpu {cluster.mem-per-cpu-mb}"

for P in ${SM_PARAMS}; do SM_ARGS="$SM_ARGS --$P {cluster.$P}"; done
echo "SM_ARGS: ${SM_ARGS}"

# our SLURM error/output paths expect a logs/ subdir in PWD
mkdir -p logs
#       
### run snakemake
# -j defines total number of jobs executed in parallel on cluster
# -n dryrun
# -p print command lines

snakemake -p \
    $* \
     --latency-wait 10 \
    -j 3000 \
    --cluster-config $(dirname $0)/cluster.slurm.json \
    --cluster "sbatch $SM_ARGS" \
    --cluster-status /scratch/mit_lieberman/scripts/slurm_status.py \
    --reason \
    --use-conda \
    --restart-times 2\
    -s Snakefile.py \
    --max-jobs-per-second 1000\
    --max-status-checks-per-second 1000\
    --keep-going\
    --rerun-incomplete\