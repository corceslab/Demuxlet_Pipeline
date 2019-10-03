#!/bin/bash
#SBATCH --job-name=mpileup_array   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=6400           # Memory per processor
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00             # Time limit hrs:min:sec
#SBATCH --output=./logs/mpileup_array_%A-%a.out    # Standard output and error log
#SBATCH --partition=howchang,sfgf   #array jobs are only submitted to howchang and sfgf to avoid job limit quotas on normal
##############################################################################################
#Ryan Corces 12/2/17
#---SLURM_SampleGenotypeMixup_mpileup.sh
#	Script expects argument #1 of input to be a full path to a file containing the manifest of BAM files to be genotyped. In format <FilePath BAM>
#	Script expects argument #2 of input to be a BED file. mpileup will only be generated in single base locations denoted in this bed file.
#	Script expects argument #3 of input to be the full path to an output directory. All files will be placed there, including log files.
#
#Usage: bash /share/PI/howchang/users/mcorces/scripts/ATAC/SLURM_SampleGenotypeMixup_mpileup.sh <Input_Manifest> <BED_File> <OUTDIR>
#
#THIS SCRIPT IS NOT MEANT TO BE CALLED FROM COMMAND LINE. THIS IS PART OF A PIPELINE.
#
#Code applies to the following software versions:
#bedtools v2.26.0
#samtools v1.5
##############################################################################################

MANIFEST=$1
REGIONS=$2
OUTDIR=$3
GENOME=$4

##############################################################################################
FASTA=${GENOMES}/${GENOME}/${GENOME}.fa
##############################################################################################

#Get the current sample name using the SLURM_ARRAY_TASK_ID variable
SAMPLE=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $1}'`
BAM=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $2}'`

# Print the task and run range
echo -e "This is task ${SLURM_ARRAY_TASK_ID}, which will perform mpileup genotyping on specific regions for sample ${SAMPLE}"

MPILEUP=${SAMPLE}.mpileup.vcf

if [ -e $OUTDIR/$MPILEUP ]
then
	echo -e "VCF mpileup genotyping file already exists! Skipping ${SAMPLE}"
else
	echo -e "No VCF mpileup genotyping file found for ${SAMPLE}. Submitting!"
	module --ignore-cache load bedtools/S2/2.26.0
	module --ignore-cache load samtools/S2/1.5
	samtools mpileup --VCF --skip-indels --uncompressed --output-tags AD --positions $REGIONS --fasta-ref $FASTA $BAM | grep -vE "#" >$OUTDIR/$MPILEUP
fi

