#!/bin/sh
# requesting the number of nodes needed
#SBATCH -A lu2022-2-24
#SBATCH -N 1
#SBATCH --tasks-per-node=10
#SBATCH --mem-per-cpu=6200
#
# job time, change for what your job farm requires
#SBATCH -t 06:00:00
#
# job name and output file names
#SBATCH -J jobFarm
#SBATCH -o res_jobFarm_%j.out
#SBATCH -e res_jobFarm_%j.out
cat $0

# set the number of jobs - change for your requirements
export NB_of_jobs=3

# Loop over the job number

for ((i=0; i<$NB_of_jobs; i++))
do
    srun -Q --exclusive --overlap -n 1 -N 1 worker_script.sh $i &> out_worker_${SLURM_JOB_ID}_${i} &
    sleep 1
done

# keep the wait statement, it is important!

wait
