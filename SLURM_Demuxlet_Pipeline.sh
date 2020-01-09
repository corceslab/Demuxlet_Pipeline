#!/bin/bash
#Ryan Corces 11/27/18
#---SLURM_Demuxlet_Pipeline.sh
#
#Usage: bash <githubFolder>/SLURM_Demuxlet_Pipeline.sh -m <Input_Manifest> ... <other options>
#	where Input_Manifest is a tab-delimited file with each line representing a single BAM file in the format <Sample Name> \t <File_path_BAM>
#
#Other options include:
#	-v  Full path to the varscan .jar file to be used for genotyping (for example /share/PI/howchang/users/mcorces/tools/varscan/VarScan.v2.4.3.jar).
#	-p 	Full directory path to a directory containing MACS2 peak calls (.narrowPeak format) for each sample in format <Sample Name>_peaks.narrowPeak
#		This option, when used, will skip the peak calling step.
#	-o 	Full directory path to the desired output directory
#	-x	Skip to this numbered step (must be an integer)
#	-g <GENOME> = This option tells the pipeline which genome to use. Currently, only supported value is "hg38"!
#
#Code applies to the following software versions:
#bedtools v2.26.0
#samtools v1.5
#macs2 2.1.1.20160309
#R 3.5.1
#
#NOTES:
#	1) You may brush up against the per-user job number submission limit depending on how many samples you are trying to run.
#		I'm not sure what the actual limit is but this certainly happens when running more than 750 samples. In this case,
#		the pipeline wont finish properly.
#	2) This script expects all of the auxilary scripts to be present in the same directory
#
#REQUIREMENTS:
#macs2 needs to runable when called as "macs2"
#utilizes my GENOMES directory / directory structure
#
##############################################################################################
#GLOBAL VARIABLES:
#Default is to clean directories and NOT call loops. Can be over-ridden by command line options
PEAKS=1 #variable that decides whether to call peaks during pipeline. 1 = TRUE, 0 = FALSE
DEPENDS="" #variable that holds the dependency for the next job. If NULL, then the job depends on nothing
SKIP=1 #variable that decides which steps to skip. 2 means you will skip step 1 and go directly to step 2
NUMCORES=20 #the number of cores used for parallel correlation calculations
##############################################################################################
# --- DO NOT CHANGE ---
#PIPELINE VARIABLES:
#retrieve the source directory of this pipeline script
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
CALL_PEAKS_PATH=${SOURCE_DIR}/SLURM_Demuxlet_CallPeaks.sh
PROCESS_PEAKS_PATH=${SOURCE_DIR}/SLURM_Demuxlet_MergePeaks.sh
GENOTYPE_PEAKS_PATH=${SOURCE_DIR}/SLURM_Demuxlet_varscanGenotype.sh
MPILEUP_REGIONS_PATH=${SOURCE_DIR}/SLURM_Demuxlet_mpileup.sh
MAKE_BIRDSEED_PATH=${SOURCE_DIR}/SLURM_Demuxlet_makeBirdseed.sh
CONCAT_BIRDSEED_PATH=${SOURCE_DIR}/SLURM_Demuxlet_concatBirdseed.R
DEMUXLET_VCF_PATH=${SOURCE_DIR}/SLURM_Demuxlet_makeDemuxletVCF.R
VARSCAN_PATH=/share/PI/howchang/users/mcorces/tools/varscan/VarScan.v2.4.3.jar
GENOME="hg38" #This can be changed at runtime using -g
##############################################################################################
#MAIN:
#Handle command line inputs using getopts
while getopts ":m:v:o:p:x:g:" opt; do
	case $opt in
	m)
		if [[ -f ${OPTARG} ]];
			then
				MANIFEST=${OPTARG}
				echo "MANIFEST: ${MANIFEST}" >&2
			else
				echo "ERROR --- -m flag observed but suggested MANIFEST does not exist: ${OPTARG}" >&2
				exit 1
		fi
		;;
	v)
		if [[ -f ${OPTARG} ]];
			then
				echo "-v flag observed. VARSCAN_PATH set to ${OPTARG}." >&2
				VARSCAN_PATH=${OPTARG}
			else
				echo "ERROR --- -v flag observedbut VARSCAN_PATH jar file not does not exist: ${OPTARG}" >&2
				exit 1
		fi
		OUTPUT_DIR=${OPTARG}
		;;
	o)
		if [[ -d ${OPTARG} ]];
			then
				echo "-o flag observed. OUTPUT_DIR exists: ${OPTARG}." >&2
			else
				echo "-o flag observed. Creating OUTPUT_DIR: ${OPTARG}" >&2
				mkdir -p $OPTARG
		fi
		OUTPUT_DIR=${OPTARG}
		;;
	p)
		if [[ -d ${OPTARG} ]];
			then
				PEAK_DIR=${OPTARG}
				PEAKS=0
				echo "-p flag observed. Skipping peak calling and using files in ${OPTARG} instead." >&2
			else
				echo "ERROR --- -p flag observed but suggested PEAK_DIR does not exist: ${OPTARG}" >&2
				exit 1
		fi
		;;
	x)
		re='^[0-9]+$'
		if [[ ${OPTARG} =~ $re ]];
			then
				SKIP=${OPTARG}
				echo "-x flag observed. Skipping to STEP# ${OPTARG}" >&2
			else
				echo "ERROR --- -x flag observed but provided argument is not a valid integer" >&2
				exit 1
		fi
		;;
	g)
		case ${OPTARG} in
			hg38 | hg19) GENOME=${OPTARG}; echo "-g flag observed. GENOME has been set to ${GENOME}." >&2;;
			*) echo "ERROR --- -g flag observed but GENOME does not match one of the supported genomes: hg38, hg19" >&2; exit 1;;
		esac
		;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done

#---------------------------------------------------------------------------------------------

#Check that all BAM files have been properly indexed before proceeding
echo Checking for BAM index files...
FAILED=0
while read LINE
do
	SAMPLE=`echo "$LINE" | awk -F"\t" '{print $1}'`
	BAM=`echo "$LINE" | awk -F"\t" '{print $2}'`
	INDEX=`echo "$BAM" | sed 's/\.bam/\.bam\.bai/g'`
	
	if [ ! -f $INDEX ];
	then
		echo -e "Index file does not exist for ${SAMPLE}"
		FAILED=1
	fi
done < $MANIFEST

if [ "$FAILED" -eq "0" ]
	then
		echo -e "All BAM files have been indexed..."
	else
		echo -e "Index files before proceeding!"
		exit 1
fi
#---------------------------------------------------------------------------------------------
#move into the output directory
cd $OUTPUT_DIR
mkdir -p ./logs

#count the number of samples
NUM_SAMPLES=`wc -l ${MANIFEST} | awk -F" " '{print $1}'`
#---------------------------------------------------------------------------------------------
#STEP1 --- Call peaks
STEP=1
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	if [ $PEAKS -eq 1 ]
	then
		#make PEAK_DIR
		PEAK_DIR=${OUTPUT_DIR}/peakCalls
		mkdir -p $PEAK_DIR
		#dependency currently set to NULL
		JOB_STRING_CALLPEAKS=$(sbatch --dependency=${DEPENDS} --array=1-${NUM_SAMPLES} ${CALL_PEAKS_PATH} ${MANIFEST} ${PEAK_DIR})
		JOB_ID_CALLPEAKS=`echo $JOB_STRING_CALLPEAKS | awk '{print $4}'`
		DEPENDS="afterok:${JOB_ID_CALLPEAKS}"
	else
		#check that narrow peak file exists for each sample
		echo Peak calling was skipped with command line input. Checking for narrowPeak files for all samples...
		FAILED=0
		while read LINE
		do
			SAMPLE=`echo "$LINE" | awk -F"\t" '{print $1}'`
			PEAK_FILE=${PEAK_DIR}/${SAMPLE}_peaks.narrowPeak
			
			if [ ! -f ${PEAK_FILE} ];
			then
				echo -e "Peak file does not exist for ${SAMPLE}"
				FAILED=1
			fi
		done < $MANIFEST
		#if not all peak files were found, then throw an error
		if [ "$FAILED" -eq "0" ]
		then
			echo -e "All peak files exist. Proceeding with pipeline..."
		else
			echo -e "Peak files missing from $PEAK_DIR. Check command line input. Consider re-calling peaks!"
			exit 1
		fi
		DEPENDS=""
	fi
else
	echo "Step ${STEP} [Peak Calling] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP2 --- Merge Peaks
STEP=2
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	if [ $PEAKS -eq 1 ]
	then
		JOB_STRING_MERGEPEAKS=$(sbatch --dependency=${DEPENDS} --output=${OUTPUT_DIR}/logs/peakMerge.log --mem=128000 --cpus-per-task=20 --time=01:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=peakMerge --wrap="bash ${PROCESS_PEAKS_PATH} ${PEAK_DIR} ${GENOME}")
		JOB_ID_MERGEPEAKS=`echo $JOB_STRING_MERGEPEAKS | awk '{print $4}'`
		DEPENDS="afterok:${JOB_ID_MERGEPEAKS}"
	else
		if [ ! -f ${PEAK_DIR}/MergedPeaks.srt.mrg.flt.narrowPeak ];
		then
			JOB_STRING_MERGEPEAKS=$(sbatch --dependency=${DEPENDS} --output=${OUTPUT_DIR}/logs/peakMerge.log --mem=128000 --cpus-per-task=20 --time=01:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=peakMerge --wrap="bash ${PROCESS_PEAKS_PATH} ${PEAK_DIR} ${GENOME}")
			JOB_ID_MERGEPEAKS=`echo $JOB_STRING_MERGEPEAKS | awk '{print $4}'`
			DEPENDS="afterok:${JOB_ID_MERGEPEAKS}"
		else
			echo -e "Merged peak file found in $PEAK_DIR. Skipping peak merge step as well."
			DEPENDS=""
		fi
	fi
else
	echo "Step ${STEP} [Peak Merging] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP3 --- Genotype Samples at Defined Peak Regions
STEP=3
VCF_DIR=${OUTPUT_DIR}/vcf
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	mkdir -p ${VCF_DIR}
	JOB_STRING_VARSCAN=$(sbatch --dependency=${DEPENDS} --array=1-${NUM_SAMPLES} ${GENOTYPE_PEAKS_PATH} ${MANIFEST} ${PEAK_DIR}/MergedPeaks.srt.mrg.flt.narrowPeak ${VCF_DIR} ${GENOME} ${VARSCAN_PATH})
	JOB_ID_VARSCAN=`echo $JOB_STRING_VARSCAN | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_VARSCAN}"
else
	echo "Step ${STEP} [Peak Genotyping] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP4 --- Merge VCF Files
STEP=4
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	JOB_STRING_MERGEVCF=$(sbatch --dependency=${DEPENDS} --output=${OUTPUT_DIR}/logs/vcfMerge.log --mem=128000 --cpus-per-task=20 --time=01:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=vcfMerge --wrap="cat ${VCF_DIR}/*.vcf | sort -k1,1 -k2,2n >${VCF_DIR}/All_SNPs_Merged.vcf")
	JOB_ID_MERGEVCF=`echo $JOB_STRING_MERGEVCF | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_MERGEVCF}"
else
	echo "Step ${STEP} [VCF Merging] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP5 --- Bedtools Merge without bookending and remove bad bed entries with ambiguous bases
STEP=5
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	JOB_STRING_BEDTOOLSMERGE=$(sbatch --dependency=${DEPENDS} --output=${OUTPUT_DIR}/logs/bedtoolsMerge.log --mem=128000 --cpus-per-task=20 --time=01:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=bedtoolsMerge --wrap="module --ignore-cache load bedtools/S2/2.26.0; bedtools merge -d -1 -c 5,6,5,6 -o distinct,distinct,count_distinct,count_distinct -i ${VCF_DIR}/All_SNPs_Merged.vcf | awk -v OFS="\'"\t"\'" "\''$7 == 1 {print $1, $2, $3, $4, $5}'\'" | awk "\''{ if ($4 !~ /[UMRWSYKVHDBN]/ && $5 !~ /[UMRWSYKVHDBN]/) print $0 }'\'">${VCF_DIR}/All_SNPs_Merged_mergedPositions.bed")
	JOB_ID_BEDTOOLSMERGE=`echo $JOB_STRING_BEDTOOLSMERGE | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_BEDTOOLSMERGE}"
else
	echo "Step ${STEP} [Bedtools Merging] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP6 --- Genotype all BAM files at these merged positions using mpileup
STEP=6
MPILEUP_DIR=${OUTPUT_DIR}/mpileup
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	mkdir -p ${MPILEUP_DIR}
	JOB_STRING_MPILEUP=$(sbatch --dependency=${DEPENDS} --array=1-${NUM_SAMPLES} ${MPILEUP_REGIONS_PATH} ${MANIFEST} ${VCF_DIR}/All_SNPs_Merged_mergedPositions.bed ${MPILEUP_DIR} ${GENOME})
	JOB_ID_MPILEUP=`echo $JOB_STRING_MPILEUP | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_MPILEUP}"
else
	echo "Step ${STEP} [Mpileup Genotyping] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP7 --- Convert mpileup VCF files to birdseed style files
STEP=7
BIRDSEED_DIR=${OUTPUT_DIR}/birdseed
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	mkdir -p ${BIRDSEED_DIR}
	JOB_STRING_BIRDSEED=$(sbatch --dependency=${DEPENDS} --array=1-${NUM_SAMPLES} ${MAKE_BIRDSEED_PATH} ${MANIFEST} ${MPILEUP_DIR} ${VCF_DIR}/All_SNPs_Merged_mergedPositions.bed ${BIRDSEED_DIR} ${SOURCE_DIR})
	JOB_ID_BIRDSEED=`echo $JOB_STRING_BIRDSEED | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_BIRDSEED}"
else
	echo "Step ${STEP} [Birdseed Conversion] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------
#STEP8 --- Create Demuxlet Input VCF
STEP=8
#Check if SKIP is less than STEP, if so run this step
if [ "$SKIP" -le "$STEP" ];
then
	JOB_STRING_DEMUXLETVCF=$(sbatch --dependency=${DEPENDS} --output=${OUTPUT_DIR}/logs/createDemuxletVCF.log --mem=25600 --cpus-per-task=4 --time=04:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=createDemuxletVCF --wrap="Rscript ${DEMUXLET_VCF_PATH} --input ${OUTPUT_DIR}/birdseed --source ${SOURCE_DIR} --outdir ${OUTPUT_DIR}")
	JOB_ID_DEMUXLETVCF=`echo $JOB_STRING_DEMUXLETVCF | awk '{print $4}'`
	DEPENDS="afterok:${JOB_ID_DEMUXLETVCF}"
else
	echo "Step ${STEP} [Make Demuxlet VCF] was skipped due to command line input (-x ${SKIP})."
	DEPENDS=""
fi
#---------------------------------------------------------------------------------------------


