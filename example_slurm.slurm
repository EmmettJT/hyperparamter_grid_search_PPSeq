#!/bin/bash

#SBATCH --job-name=p1			#name of the job to find when calling >>>sacct or >>>squeue
#SBATCH --nodes=1				#number of nodes, i.e. computers, to request off the cluster (nodes typically have ~20 singled threaded cores)
#SBATCH --ntasks=1				#how many independent script you are hoping to run 
#SBATCH --cpus-per-task=20			#how many threads to multithread across (no point more than number of cores available. also, you cannot thread across nodes) 
#SBATCH --time=180:00:00				#compute time
#SBATCH --mem=8gb				#memory to request
#SBATCH --output=run_out.log		#where to save output log files (julia script prints here)  
#SBATCH --error=run_out.err		#where to save output error files
#SBATCH --array=0-44            #task ID array for array scripting (can be passed to script below with command line argument --slurm-array-task-id $SLURM_ARRAY_TASK_ID)

pwd; hostname; date

export JULIA_NUM_THREADS=10

julia SimpleSweep_P1.jl --data-directory "/ceph/sjones/projects/sequence_squad/revision_data/organised_data/ppseq_data/prepared_data/hyper_search/short_data_simple_sweeps/" --data-filename seq006_1_3 --num-threads 10 --results-directory "/ceph/sjones/projects/sequence_squad/revision_data/organised_data/ppseq_data/output_data/HyperSearch/simple_sweeps/p1/" --results-name grid_search_$SLURM_ARRAY_TASK_ID --slurm-array-task-id $SLURM_ARRAY_TASK_ID


sstat -j $SLURM_JOB_ID.batch --format=JobID,MaxVMSize,
