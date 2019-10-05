#!/bin/bash
#SBATCH --job-name=varscanGenotype_array   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=25600           # Memory per processor
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=./logs/varscanGenotype_array_%A-%a.out    # Standard output and error log
#SBATCH --partition=howchang,sfgf   #array jobs are only submitted to howchang and sfgf to avoid job limit quotas on normal
##############################################################################################
#Ryan Corces 11/27/18
#---SLURM_Demuxlet_varscanGenotype.sh
#	Script expects argument #1 of input to be a full path to a file containing the manifest of BAM files to be genotyped. In format <FilePath BAM>
#	Script expects argument #2 of input to be a BED file. Genotyping will only be performed in regions annotated in the input BED file to save time.
#	Script expects argument #3 of input to be the full path to an output directory. All files will be placed there, including log files.
#
#Usage: bash SLURM_Demuxlet_varscanGenotype.sh <Input_Manifest> <BED_File> <OUTDIR>
#
#THIS SCRIPT IS NOT MEANT TO BE CALLED FROM COMMAND LINE. THIS IS PART OF A PIPELINE.
#
#Code applies to the following software versions:
#bedtools v2.26.0
#samtools v1.5
##############################################################################################

MANIFEST=$1
PEAKS=$2
OUTDIR=$3
GENOME=$4

##############################################################################################
FASTA=${GENOMES}/${GENOME}/${GENOME}.fa
##############################################################################################



#Get the current sample name using the SLURM_ARRAY_TASK_ID variable
SAMPLE=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $1}'`
BAM=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $2}'`

# Print the task and run range
echo -e "This is task ${SLURM_ARRAY_TASK_ID}, which will perform Varscan genotyping for sample ${SAMPLE}"

VCF=${SAMPLE}.vcf

if [ -e ${OUTDIR}/${VCF} ]
then
	echo -e "VCF genotyping file already exists! Skipping ${SAMPLE}"
else
	echo -e "No VCF genotyping file found for ${SAMPLE}. Submitting!"
	module --ignore-cache load bedtools/S2/2.26.0
	module --ignore-cache load samtools/S2/1.5
	samtools mpileup -l $PEAKS -B -f ${FASTA} $BAM | java -jar /share/PI/howchang/users/mcorces/tools/varscan/VarScan.v2.4.3.jar mpileup2snp --min-coverage 5 --min-reads2 2 --min-var-freq 0.1 --strand-filter 1 --output-vcf 1 | grep -vE \"^#\" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >${OUTDIR}/${VCF} 
fi
