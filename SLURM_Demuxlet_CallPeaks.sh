#!/bin/bash
#SBATCH --job-name=CallPeaks_array   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=25600           # Memory per processor
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=./logs/CallPeaks_array_%A-%a.out    # Standard output and error log
#SBATCH --partition=howchang,sfgf   #array jobs are only submitted to howchang and sfgf to avoid job limit quotas on normal
##############################################################################################
#Ryan Corces 11/27/18
#---SLURM_SampleGenotypeMixup_CallPeaks.sh
#	Script expects argument #1 of input to be a full path to a file containing the manifest of BAM files to be processed. In format <BAMFilePath>
#	Script expects argument #2 of input to be a full path to the output directory.
#
#THIS SCRIPT IS NOT MEANT TO BE CALLED FROM COMMAND LINE. THIS IS PART OF A PIPELINE.
#
#Usage: bash /share/PI/howchang/users/mcorces/scripts/ATAC/SLURM_SampleGenotypeMixup_CallPeaks.sh <Full_Path_To_Input_Manifest> <Full_Path_To_Output_Directory>
#
#Code applies to the following software versions:
#macs2 2.1.1.20160309
##############################################################################################
#load python 2.7.13 since MACS2 isnt compatible with python3 at the moment
PATH=/share/PI/howchang/users/mcorces/tools/python/virtualenv/python-2.7.13/bin:$PATH
##############################################################################################


MANIFEST=$1
OUTDIR=$2

#Get the current sample name using the SLURM_ARRAY_TASK_ID variable
SAMPLE=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $1}'`
BAM=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $2}'`

# Print the task and run range
echo -e "This is task ${SLURM_ARRAY_TASK_ID}, which will perform MACS2 peak calling for sample ${SAMPLE}"

if [ -e $OUTDIR/${SAMPLE}_peaks.narrowPeak ]
then
	echo -e "Peak file already exists! Skipping ${SAMPLE}"
else
	macs2 callpeak --nomodel -t $BAM -n ${OUTDIR}/${SAMPLE} --nolambda --keep-dup all --call-summits
fi

