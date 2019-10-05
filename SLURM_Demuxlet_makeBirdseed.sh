#!/bin/bash
#SBATCH --job-name=makeBirdseed_array   # Job name
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=6400           # Memory per processor
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00             # Time limit hrs:min:sec
#SBATCH --output=./logs/makeBirdseed_array_%A-%a.out    # Standard output and error log
#SBATCH --partition=howchang,sfgf   #array jobs are only submitted to howchang and sfgf to avoid job limit quotas on normal
##############################################################################################
#Ryan Corces 7/5/18
#---SLURM_Demuxlet_makeBirdseed.sh
#This script is meant to be used when you are only comparing ATAC data to itself rather than to SNP array data
#
#	Script expects argument #1 of input to be a full path to a file containing the manifest of vcf files to be converted to birdseed format. In format <FilePath VCF>
#	Script expects argument #2 of input to be the full path to the mpileup file directory.
#	Script expects argument #3 of input to be an Affymetrix library type file with information about all of the probes. File should be 180122_AFFX_GenomeWideSNPs_filteredFinal.txt.
#			AFFX FILE FORMAT: <Probe_Set_ID> <dbSNP_RS_ID> <hg19_chr> <hg19_start> <hg19_stop> <hg19_strand> <Allele_A> <Allele_B> <hg38_chr> <hg38_start> <hg38_stop> <hg38_strand> <type> <freqCNV> <par>
#								First row should have a header!
#	Script expects argument #4 of input to be the full path to an output directory. All files will be placed there, including log files.
#
#Usage: bash SLURM_Demuxlet_makeBirdseed.sh <Input_Manifest> <affymetrix file> <OUTDIR>
#
#THIS SCRIPT IS NOT MEANT TO BE CALLED FROM COMMAND LINE. THIS IS PART OF A PIPELINE.
#
#Code applies to the following software versions:
#
##############################################################################################

##############################################################################################

MANIFEST=$1
MPILEUP_DIR=$2
AFFX=$3
OUTDIR=$4
SOURCE_DIR=$5

#Get the current sample name using the SLURM_ARRAY_TASK_ID variable
SAMPLE=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $1}'`
BAM=`cat ${MANIFEST} | sed -n ${SLURM_ARRAY_TASK_ID}p | awk -F"\t" '{print $2}'`

# Print the task and run range
echo -e "This is task ${SLURM_ARRAY_TASK_ID}, which will perform vcf to birdseed conversion for sample ${SAMPLE}"

BIRD=${SAMPLE}.affx
MPILEUP=${MPILEUP_DIR}/${SAMPLE}.mpileup.vcf

if [ -e $OUTDIR/$BIRD ]
then
	echo -e "Birdseed-SNP file already exists! Skipping ${SAMPLE}"
else
	echo -e "No Birdseed-SNP file found for ${SAMPLE}. Submitting!"
	Rscript ${SOURCE_DIR}/SLURM_Demuxlet_makeBirdseed.R --input $MPILEUP --affx $AFFX --outdir $OUTDIR
fi

