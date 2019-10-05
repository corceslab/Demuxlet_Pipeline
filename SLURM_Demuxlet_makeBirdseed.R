################################################################################################################
#1/22/18
#Last Updated 11/27/18
#Ryan Corces
#SLURM_SampleGenotypeMixup_makeBirdseed.R
#
#This program is meant to convert mpileup files to be compatible with birdseed data files. It takes a VCF file from
#samtools mpileup generated using SLURM_ATAC_mpileupFromBed.sh and makes "confident" birdseed-style genotype calls
#where possible. This means that it filters SNPs based off a minimum depth of <depthCutoff> and then performs a genotype call
#for "AA" "AB" "BB" or "no call". These types of calls are only made on a subset of the SNPs and the remainder are given
#"no call". 
#
#NOTES:
#usage - 
#example - 
#Rscript /path/to/SLURM_SampleGenotypeMixup_makeBirdseed.R --input </path/to/inputMpileupMatrix.txt> --affx </path/to/affxSNP6_metadata.txt> --outdir </path/to/outputDirectory>
#where:
#--input </path/to/inputMpileupMatrix.txt>
#--affx </path/to/affxSNP6_metadata.txt>
#--outdir </path/to/outputDirectory>
#
#output files are automatically named by changing the filename suffix from ".vcf" to ".snp" or ".affx" depending
#
#The idea is to create birdseed-style output files and then concatenate them into a single matrix which can
#then be correlated on a sample by sample basis with birdseed-style SNP6 from TCGA. To concatenate the individual
#birdseed-style output files, use script TCGA_concatenateBirdseed.R
#
################################################################################################
## Load dependencies
library(optparse)
library(seqinr)
################################################################################################
#Global Variables and functions
depthCutoff <- 20 #cutoff for depth in ATAC-seq data to make a genotyping call

################################################################################################
## Input variables

print(paste("Start",date(),sep=" - "))

## uses optparse package to read in command line arguments

option_list = list(
  make_option(c("--input"), action="store", default=NULL, type='character',
              help="/path/to/inputMpileupMatrix.txt"),
  make_option(c("--affx"), action="store", default=NULL, type='character',
              help="/path/to/affxSNP6_metadata.txt"),
  make_option(c("--outdir"), action="store", default=NULL, type='character',
              help="/path/to/outputDirectory")
)

opt = parse_args(OptionParser(option_list=option_list))

if(is.null(opt$input) || is.null(opt$affx) || is.null(opt$outdir)){
  # print usage information if an argument is not provided
  stop("\n Usage: Rscript /path/to/SLURM_SampleGenotypeMixup_makeBirdseed.R --input </path/to/inputMpileupMatrix.txt> --affx </path/to/affxSNP6_metadata.txt> --outdir </path/to/outputDirectory> \n
       All arguments are required! \n\n")
}

setwd(opt$outdir)
AFFX_file <- opt$affx
mpileup_file <- opt$input
snp_outfile <- gsub(".vcf",".snp",paste(opt$outdir,"/",basename(mpileup_file),sep=""))
affx_outfile <- gsub(".vcf",".affx",paste(opt$outdir,"/",basename(mpileup_file),sep=""))

print(paste("input = ",opt$input),sep="")
print(paste("affx data = ",opt$affx),sep="")
print(paste("output dir = ",opt$outdir),sep="")
print(paste("snp outfile = ",snp_outfile),sep="")
print(paste("affx outfile = ",affx_outfile),sep="")


################################################################################################

#---- This Code is meant for NON-TCGA data where correlations are being made just from ATAC-seq. If using SNP arrays, use the above snippet.
#Read in AFFX library file and remove unneccessary columns
affx <- read.table(AFFX_file, header = FALSE, row.names=NULL, check.names = TRUE, colClasses= c("character","integer","integer","character","character"), sep="\t")
colnames(affx) <- c("chr","start","stop","Allele_A","Allele_B")
affx$Probe_Set_ID <- rownames(affx)

#make rownames indexable by joining chr and stop position
rownames(affx) <- paste(affx$chr,affx$stop,sep="-")

#read in mpileup file
vcf <- read.table(mpileup_file, header = FALSE, row.names = NULL, colClasses= c("character","integer","character","character","character","character","character","character","character","character"))
colnames(vcf) <- c("chr","pos","1","ref","varString","2","3","4","5","PL_AD")

AD <- matrix(unlist(strsplit(vcf$PL_AD, ":")), nrow=nrow(vcf), byrow=TRUE)[,2]

#make matrix that is more usefully formatted for snp interpretation
snps <- data.frame(matrix(0,nrow = nrow(vcf), ncol = 6))
colnames(snps) <- c("chr","pos","A","C","G","T")

snps$chr <- vcf$chr
snps$pos <- vcf$pos

for (j in 1:nrow(vcf))
{
  if(j %% 10000 == 0)
  {
    print(j)
  }
  alleles <- c(vcf$ref[j],unlist(strsplit(vcf$varString[j],",")))
  counts <- unlist(strsplit(AD[j],","))
  for(k in 1:length(alleles))
  {
    if (alleles[k] == "A")
    {
      snps$A[j] <- counts[k]
    }else if (alleles[k] == "C") {
      snps$C[j] <- counts[k]
    }else if (alleles[k] == "G") {
      snps$G[j] <- counts[k]
    }else if (alleles[k] == "T") {
      snps$T[j] <- counts[k]
    }
  }
}
rownames(snps) <- paste(snps$chr,snps$pos,sep="-")

write.table(snps, snp_outfile,sep="\t", col.names = NA, quote = FALSE)

#truncate snp table based off absolute depth of seq at position
snpsTruncated <- snps[which(rowSums(data.matrix(snps[,c("A","C","G","T")])) > depthCutoff),]

#make new affy-based matrix with one entry per affx snp probe
affxSNPs <- data.frame(matrix(NA,nrow = nrow(affx), ncol = 8))
colnames(affxSNPs) <- c("probeID","chr","pos","allele_A","allele_A_counts","allele_B","allele_B_counts","call")
affxSNPs$probeID <- affx$Probe_Set_ID
affxSNPs$chr <- affx$chr
affxSNPs$pos <- affx$stop
affxSNPs$allele_A <- affx$Allele_A
affxSNPs$allele_B <- affx$Allele_B
affxSNPs$call <- rep("-1",nrow(affxSNPs))
rownames(affxSNPs) <- rownames(affx)

#for each SNP in snps
for (y in 1:nrow(snpsTruncated))
{
  if(y %% 1000 == 0)
  {
    print(y)
  }
  #retrieve affx index position and check that it is a unique match bewtween ATAC snp and AFFX snp
  affxIndex <- which(rownames(affx) == rownames(snpsTruncated)[y])
  try(if(length(affxIndex) != 1) stop(paste("ERROR - Duplicate entries identified in AFFX for snp ",rownames(snpsTruncated)[y],"!",sep="")))
  
  #make call 0=AA, 1=AB, 2=BB, -1=no call
  #if either A or B counts is zero, then the call is easy
  #if the difference in counts between A and B is lower than 50% of the depth, then its heterzygous
  alleleA <- affx$Allele_A[affxIndex]
  alleleB <- affx$Allele_B[affxIndex]
  alleleA_counts <- as.numeric(snpsTruncated[y,alleleA])
  alleleB_counts <- as.numeric(snpsTruncated[y,alleleB])
  
  #make call 0=AA, 1=AB, 2=BB, -1=no call
  #minimum depth = 6 reads
  if(alleleA_counts == 0) 
  {
    currentCall <- 2 #if A = 0 counts, then it must be BB given minimum depth
  } else if (alleleB_counts == 0) {
    currentCall <- 0 #if B = 0 counts, then it must be AA given minimum depth
  } else if ((abs(alleleA_counts - alleleB_counts)/(alleleA_counts + alleleB_counts)) < 0.5) {
    currentCall <- 1 #if the counts for A and B are close to each other, it is heterozygous AB
  } else if (alleleA_counts > alleleB_counts) {
    currentCall <- 0 #otherwise, if A is greater than B, then it is AA
  } else if (alleleA_counts < alleleB_counts) {
    currentCall <- 2 #otherwise, if B is greater than A, then it is BB
  }#if no call was made by now, the call remains -1
  
  affxSNPs$call[affxIndex] <- currentCall
  affxSNPs$allele_A_counts[affxIndex] <- alleleA_counts
  affxSNPs$allele_B_counts[affxIndex] <- alleleB_counts
}

write.table(affxSNPs, affx_outfile,sep="\t", col.names = NA, quote = FALSE)
print(paste("Start",date(),sep=" - "))