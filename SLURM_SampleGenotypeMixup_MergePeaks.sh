#!/bin/bash
#Ryan Corces 11/27/18
#---SLURM_SampleGenotypeMixup_MergePeaks.sh
#	Script expects argument #1 of input to be a path to the directory containing narrowPeak files
#	SCRIPT EXPECTS GENOME TO BE hg38
#Usage: bash /share/PI/howchang/users/mcorces/scripts/ATAC/ATAC_processPeaks.sh <PeakDIR>
#
#THIS SCRIPT IS NOT MEANT TO BE CALLED FROM COMMAND LINE. THIS IS PART OF A PIPELINE.
#
#Code applies to the following software versions:
#macs2 2.1.1.20160309
#bedtools_2.26.0

echo Loading Modules...
module load bedtools/S2/2.26.0
module load libpng/S2/1.2.57
module load ucscTools/6.21.17

PEAKS_DIR=$1
GENOME=$2

BLACKLIST=${GENOMES}/${GENOME}/${GENOME}.blacklist.bed
NUMTS=${GENOMES}/${GENOME}/${GENOME}_NumtS_Regions.bed

cd $PEAKS_DIR

echo Concatenating peak files...
cat *_peaks.narrowPeak >MergedPeaks.narrowPeak

echo Sorting merged peak file...
bedSort MergedPeaks.narrowPeak MergedPeaks.srt.narrowPeak

echo Merging sorted peak file...
bedtools merge -i MergedPeaks.srt.narrowPeak -c 6,7,8,9 -o count,max,max,max >MergedPeaks.srt.mrg.narrowPeak

echo Filtering merged peak file...
bedtools intersect -v -a MergedPeaks.srt.mrg.narrowPeak -b $BLACKLIST $NUMTS | grep -v chrM >MergedPeaks.srt.mrg.flt.narrowPeak
