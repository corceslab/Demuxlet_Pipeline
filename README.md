# SampleGenotypeMixup
Pipeline used to identify sample mixups based on genotyping of NGS reads.

This pipeline contains a controller script [SLURM_SampleGenotypeMixup_Pipeline.sh] and a series of auxilary bash and R scripts that are called at run time and queued as jobs with SLURM dependencies.

## PIPELINE COMPONENTS:

	SLURM_SampleGenotypeMixup_Pipeline.sh
		STEP1: Call peaks in data. This pipeline was designed for ATAC-seq. Calling peaks lets you quickly identify regions of high read count which could be used for genotyping purposes.
			SLURM_SampleGenotypeMixup_CallPeaks.sh
		STEP2: Merge peaks called in all samples. Merging the peak regions creates a union set of regions to be genotyped across all samples.
			SLURM_SampleGenotypeMixup_MergePeaks.sh
		STEP3: Genotyping the merged peak regions in every sample using Varscan.
			SLURM_SampleGenotypeMixup_varscanGenotype.sh
		STEP4: Merge all resultant VCF files from Varscan. This doesnt use a script, just a multipart bash one-liner
		STEP5: Use bedtools to merge (without bookending) all of the unique single-base positions that should be re-genotyped across all samples. This represents the union set of putative positions containing SNPs.
		STEP6: Use samtools mpileup to genotype at these specific positions across all samples.
			SLURM_SampleGenotypeMixup_mpileup.sh
		STEP7: Convert the allelic depth at each of these positions into a "birdseed"-style genotyping call using an R script. Bash script serves as a wrapper.
			SLURM_SampleGenotypeMixup_makeBirdseed.sh
			SLURM_SampleGenotypeMixup_makeBirdseed.R
		STEP8: Concatenate all of the birdseed-style calls into a single quaternary matrix (0, 1, or 2 for the genotyping calls or -1 for no call).
			SLURM_SampleGenotypeMixup_concatBirdseed.R
		STEP9: Correlate all by all and output plots and summary files to analyze sample mixups.
			SLURM_SampleGenotypeMixup_SelfMatrixCorrelation.R

## REQUIREMENTS:

The auxilary scripts utilize my GENOMES directory / directory structure which is shared on Sherlock. This relies on an environment variable called GENOMES that points to a directory that looks like this and has the following files:

	${GENOMES}/hg38/hg38.fa
	${GENOMES}/hg38/hg38_Common_SNPs.txt
	${GENOMES}/hg38/hg38.blacklist.bed
	${GENOMES}/hg38/hg38_NumtS_Regions.bed

## DEPENDENCIES:
This code applies to the following software versions which must be loadable with the given commands:

	bedtools v2.26.0	- module load bedtools/S2/2.26.0
	samtools v1.5		- module load samtools/S2/1.5
	macs2 2.1.1.20160309	- Not loaded. Instead, macs2 needs to runable when called as "macs2".
	R 3.5.1			- module load R/3.5.1
	libpng v1.2.57		- module load libpng/S2/1.2.57 (required for certain tools on Sherlock2)
	UCSC Tools v 6.21.17	- module load ucscTools/6.21.17


This code requires the following R packages to be installed:

	optparse, foreach, doParallel, ggplot2, seqinr, matrixStats

## NOTES:
1) You may brush up against the per-user job number submission limit depending on how many samples you are trying to run. I'm not sure what the actual limit is but this certainly happens when running more than 750 samples. In this case, the pipeline wont finish properly.
2) This script expects all of the auxilary scripts to be present in the same directory
3) Only hg38 is currently supported

## USAGE:
	bash /path/to/SLURM_SampleGenotypeMixup_Pipeline.sh -m <Input_Manifest> ... <other options>
	where Input_Manifest is a tab-delimited file with each line representing a single BAM file in the format <Sample Name> \t <File_path_BAM> \t <condition> \t <bioRepID>
	where condition is the sample type (for example, cancer type from TCGA, or brain region from AD, or brain cohort). Just some way of grouping samples for plotting
	where bioRepID is an identifier that is shared across all samples from the same biological donor. These IDs are used to make the "in groups" for correlation
	
	Other options include:
	-p 	Full directory path to a directory containing MACS2 peak calls (.narrowPeak format) for each sample in format <Sample Name>_peaks.narrowPeak
		This option, when used, will skip the peak calling step.
	-o 	Full directory path to the desired output directory
	-x	Skip to this numbered step (must be an integer)

## TEST USAGE:
To download and test the pipeline, follow these steps:

0) Make sure all of the required R libraries are installed and you have satisfied the requirements and dependencies.
1) Clone the github repository using: `git clone https://github.com/rcorces/SampleGenotypeMixup.git`
2) Make a folder to perform the analysis.
3) Move into this folder and download and unzip the test data using `wget https://s3.amazonaws.com/changseq/RyanCorces/Shared/DIAN-CAUD_test_data.zip` and `unzip DIAN-CAUD_test_data.zip`. You should end up with 2 folders and 1 test manifest text file. `DIAN-CAUD_test_input` should contain 8 bam files and 8 bam index files. `DIAN-CAUD_test_output` should contain 2 text files and 1 PDF. The files in this `DIAN-CAUD_test_output` folder contain examples of what the output should look like after running the pipeline for these test samples.
4) Edit the second column of the test manifest file to provide the correct full path to each bam file.
5) Run the test data through the pipeline using `bash /full/path/to/SLURM_SampleGenotypeMixup_Pipeline.sh -m /full/path/to/DIAN-CAUD_test_manifest.txt -o /full/path/to/output/directory/`.
6) Compare your output files to those in the `DIAN-CAUD_test_output` folder.
