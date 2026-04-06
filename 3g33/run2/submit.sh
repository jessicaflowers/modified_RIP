#!/bin/bash
#$ -S /bin/bash         #-- the shell for the job
#$ -q gpu.q             #-- use the gpu queue
##$ -pe smp 2            #-- when using the gpu queue with more than 1 gpu, this tells the scheduler how many gpus are needed 
#$ -o sge_logs          #-- output directory (fill in)
#$ -j y                 #-- tell the system that the STDERR and STDOUT should be joined
#$ -cwd                 #-- tell the job that it should start in your working directory
#$ -l mem_free=4G       #-- submits on nodes with enough free memory
#$ -l h_rt=2:00:00    #-- runtime limit - max 2 weeks == 336 hours
#$ -R yes               #-- SGE host reservation
##$ -l hostname=!(‘qb3-atgpu*‘|'qb3-atgpu**‘|'qb3-iogpu4'|'qb3-idgpu11'|'qb3-idgpu15')
#$ -l hostname=!('qb3-idgpu11'|'qb3-idgpu6'|'qb3-idgpu10'|'qb3-iogpu1')
#$ -N rip_3g33_run2
#$ -t 1-150
##$ -t 1-1  
#$ -tc 1

# Load modules
module load mpi
module load Sali
module load cuda/12.5 # cuda12.X for Gromacs2023/4/5


# Environment Variables
export OMP_NUM_THREADS=8
export CUDA_VISIBLE_DEVICES=$SGE_GPU
echo "Starting Task ID: $SGE_TASK_ID"
echo "Running on Host: $HOSTNAME"
echo "Assigned GPU ID: $SGE_GPU"

# Run job

# starting at SGE_TASK_ID, search forward for the first existing
# pulse_res_<n> directory and run the dynamics there
MAX_SEARCH=310
START_ID=${SGE_TASK_ID}
FOUND_DIR=""
for (( offset=0; offset<MAX_SEARCH; offset++ )); do
	IDX=$((START_ID + offset))
	PULSE_DIR="pulse_res_${IDX}"
	if [ -d "${PULSE_DIR}" ]; then
		FOUND_DIR="${PULSE_DIR}"
		break
	fi
done

if [ -z "${FOUND_DIR}" ]; then
	echo "ERROR: no pulse_res_<n> directory found in range ${START_ID}..$((START_ID+MAX_SEARCH-1)). Exiting."
	exit 1
fi

echo "Found pulse directory: ${FOUND_DIR} (from starting ID ${START_ID})"
cd "${FOUND_DIR}"
# copy run_dynamics.sh into the pulse dir if not already present so the script
if [ ! -f run_dynamics.sh ]; then
	cp -n ../run_dynamics.sh . || true
fi
echo "Running in: $(pwd)"
bash run_dynamics.sh

# Job summary
printenv
qstat -j $JOB_ID
