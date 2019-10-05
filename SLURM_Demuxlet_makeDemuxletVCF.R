################################################################################################################
#10/05/19
#Last Updated 10/05/19
#Ryan Corces
#SLURM_Demuxlet_makeDemuxletVCF.R
#
#This program is meant to create a sample VCF where each row is a different SNP and each column is a different
#sample and the values contain mock genotype probabilities based on birdseed-style genotyping calls
#
#NOTES:
#usage - 
#Rscript SLURM_Demuxlet_makeDemuxletVCF.R --input </path/to/birdseed/directory/> --source </path/to/directory/containing/this/Rscript/> --outdir </path/to/outputDirectory>
#where:
# --input </path/to/birdseed/directory/> -- directory contains ".affx" files
# --source </path/to/directory/containing/this/Rscript/>
# --outdir </path/to/outputDirectory>
#
#Optional Arguments:
# --techRep </path/to/manifest/file/containing/techRepInfo.txt> -- Column 3 should contain the tech rep info for each sample 
#
#output files are automatically named by changing the filename suffix from ".vcf" to ".snp" or ".affx" depending
#
#The idea is to create birdseed-style output files and then concatenate them into a single matrix which can
#then be correlated on a sample by sample basis with birdseed-style SNP6 from TCGA. To concatenate the individual
#birdseed-style output files, use script TCGA_concatenateBirdseed.R
#
################################################################################################
## Load dependencies
library(matrixStats)
library(ArchRx)
################################################################################################
#Global Variables and functions
################################################################################################
## Input variables

print(paste("Start",date(),sep=" - "))

## uses optparse package to read in command line arguments

option_list = list(
  make_option(c("--input"), action="store", default=NULL, type='character',
              help="/path/to/birdseed/directory/"),
  make_option(c("--source"), action="store", default=NULL, type='character',
              help="/path/to/directory/containing/this/Rscript/"),
  make_option(c("--techRep"), action="store", default=NULL, type='character',
              help="/path/to/manifest/file/containing/techRepInfo.txt"),
  make_option(c("--outdir"), action="store", default=NULL, type='character',
              help="/path/to/outputDirectory")
)

opt = parse_args(OptionParser(option_list=option_list))

if(is.null(opt$input) || is.null(opt$source) || is.null(opt$outdir)){
  # print usage information if an argument is not provided
  stop("\n Usage: Rscript /path/to/SLURM_Demuxlet_makeDemuxletVCF.R --input </path/to/birdseed/directory/> --source </path/to/directory/containing/this/Rscript/> --outdir </path/to/outputDirectory> \n
       Input, source, and output directories are required! \n\n")
}

in_dir <- opt$input
source_dir <- opt$source
out_dir <- opt$outdir

print(paste("input dir = ",opt$input),sep="")
print(paste("output dir = ",opt$outdir),sep="")

dir.create(path = out_dir, showWarnings = FALSE)

################################################################################################
setwd(in_dir)

fileList <- list.files(pattern = ".*affx")

#read in a template file to get the formating of the depth and call matrices correct up front
template <- read.table(fileList[1], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
numRows <- nrow(template)
numCols <- length(fileList)

#innitiate the depth and call matrices to store the total depth and genotype call at each position
depth <- as.data.frame(matrix(nrow = numRows, ncol = numCols, NA))
colnames(depth) <- gsub(pattern = ".mpileup.affx", replacement = "", x = fileList)
rownames(depth) <- rownames(template)
call <- depth

#read in each affx file and add to the depth and call matrices
for (i in 1:length(fileList)) {
  print(i)
  current <- read.table(fileList[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  current$allele_A_counts[which(is.na(current$allele_A_counts))] <- 0
  current$allele_B_counts[which(is.na(current$allele_B_counts))] <- 0
  
  depth[,i] <- current$allele_A_counts + current$allele_B_counts
  call[,i] <- current$call
}

saveRDS(object = depth, file = paste0(out_dir,"/demuxlet_SNP-depth.rds"))
saveRDS(object = call, file = paste0(out_dir,"/demuxlet_SNP-call.rds"))

nonZero_depth <- depth[which(rowSums(depth) > 0),]
nonZero_call <- call[which(rowSums(depth) > 0),]
saveRDS(object = nonZero_depth, file = paste0(out_dir,"/demuxlet_SNP-depth_nonZeroDepth.rds"))
saveRDS(object = nonZero_call, file = paste0(out_dir,"/demuxlet_SNP-call_nonZeroDepth.rds"))

#--------------------------------------------------------------------
#Merge techReps if necessary
#THIS SECTION HAS NOT BEEN FIXED YET
# if(!is.null(techRep)) {

#   #Merge data from tech reps to only keep positions that are able to be called in both reps
#   #make matrices to hold data corresponding to tech rep comparisons
#   techRep_depth <- as.data.frame(matrix(nrow = nrow(nonZero_depth), ncol = ncol(nonZero_depth)/2, NA))
#   rownames(techRep_depth) <- rownames(nonZero_depth)
#   colnames(techRep_depth) <- unique(gsub(pattern = "_Rep.", replacement = "", x = colnames(nonZero_depth), perl = TRUE))
  
#   techRep_call <- as.data.frame(matrix(nrow = nrow(nonZero_depth), ncol = ncol(nonZero_depth)/2, NA))
#   rownames(techRep_call) <- rownames(nonZero_depth)
#   colnames(techRep_call) <- unique(gsub(pattern = "_Rep.", replacement = "", x = colnames(nonZero_depth), perl = TRUE))
  
#   #make dummy vectors to store depth of matches and mismatches
#   match_depth <- as.data.frame(matrix(nrow=0, ncol = 4, NA))
#   colnames(match_depth) <- c("id","depth","call","type")
  
#   #go through all tech rep pairs and check matches vs mismatches
#   for (i in 0:(length(techRep_depth) - 1)) {
#     x <- (i * 2) + 1 #techRep 1 index
#     y <- (x + 1) #techRep 2 index
#     #subset the calls to those where both tech reps can make a call
#     subset <- intersect(which(nonZero_call[,x] != -1), which(nonZero_call[,y] != -1))
#     subset_depth <- nonZero_depth[subset,]
#     subset_call <- nonZero_call[subset,]
#     #identify indicies of matches and mismatches
#     matches <- which(subset_call[,x] == subset_call[,y])
#     mismatches <- which(subset_call[,x] != subset_call[,y])
#     #update the corresponding cells in the techRep matrices with the values identified in the subsets for matches and mismatches
#     #for matches, the call is (obviously) the matched call
#     techRep_call[rownames(subset_depth)[matches],(i+1)] <- subset_call[matches,x]
#     #for matches, the depth stored is the minimum of the depths of the two tech reps
#     techRep_depth[rownames(subset_depth)[matches],(i+1)] <- rowMins(as.matrix(data.frame(A = subset_depth[matches,x], B = subset_depth[matches,y])))
#     #for mismatches, the value of the call is made -1
#     techRep_call[rownames(subset_depth)[mismatches],(i+1)] <- -1
#     #for mismatches, the depth stored is the minimum of the depths of the two tech reps
#     techRep_depth[rownames(subset_depth)[mismatches],(i+1)] <- rowMins(as.matrix(data.frame(A = subset_depth[mismatches,x], B = subset_depth[mismatches,y])))
    
#     #for each techRep pair, add the match and mismatch data to the depth matrix
    
#     match_bind <- data.frame(id = rownames(subset_depth)[matches],
#                              depth = rowMins(as.matrix(data.frame(A = subset_depth[matches,x], B = subset_depth[matches,y]))),
#                              call = subset_call[matches,x],
#                              type = rep("match",length(matches)))
#     mismatch_bind <- data.frame(id = rownames(subset_depth)[mismatches],
#                                 depth = rowMins(as.matrix(data.frame(A = subset_depth[mismatches,x], B = subset_depth[mismatches,y]))),
#                                 call = subset_call[mismatches,x],
#                                 type = rep("mismatch",length(mismatches)))
    
#     match_depth <- rbind(match_depth, match_bind, mismatch_bind)
    
#     print(paste0(x, ": ",as.character(round(length(matches) / (length(matches) + length(mismatches)) * 100, digits = 2))," matched."))
#   }
  
#   saveRDS(techRep_call, file = paste0(out_dir,"/demuxlet_techRepsMerged_GenotypeCall.rds"))
#   saveRDS(techRep_depth, file = paste0(out_dir,"/demuxlet_techRepsMerged_GenotypeDepth.rds"))
  
# }

#--------------------------------------------------------------------

#preallocate a matrix to hold the number of differences between a given sample and all other samples
distance <- as.data.frame(matrix(nrow = ncol(nonZero_call), ncol = ncol(nonZero_call), 0))
rownames(distance) <- colnames(nonZero_call)
colnames(distance) <- colnames(nonZero_call)
#distance is a matrix where each cell represents the number of genotype call differences between the column and row
for (a in 1:nrow(distance)) {
  print(a)
  for(b in a:ncol(distance)) {
    #only one of the two following distance calculations should be uncommented
    #distance could be calculated based only on positions where BOTH samples can make a call
    # call_made <- intersect(which(nonZero_call[,a] >= 0), which(nonZero_call[,b] >= 0))
    #or distance could be calculated based on positions where EITHER sample can make a call
    call_made <- union(which(nonZero_call[,a] >= 0), which(nonZero_call[,b] >= 0))
    
    distance[a,b] <- length(which(nonZero_call[call_made,a] != nonZero_call[call_made,b]))
    distance[b,a] <- distance[a,b]
  }
}
saveRDS(distance, file = paste0(out_dir,"/demuxlet_GenotypeDistance_anyCall.rds"))
write.table(x = distance, file = paste0(out_dir,"/demuxlet_GenotypeDistance_anyCall.txt"), quote = FALSE, sep = "\t", col.names = NA)

#--------------------------------------------------------------------
#Generate "genotype probabilities" for use with demuxlet

#find all of the rows that have at least one genotyping call in one of the used samples
list <-c()
for (i in 1:nrow(nonZero_call)) {
  if ((i %% 1000) == 0) {
    print(i)
  }
  if(length(which(is.na(nonZero_call[i,]))) != ncol(nonZero_call)) {
    list <- c(list,i)
  }
}
#subset nonZero_call to only contain rows that had at least one sample with a genotyping call
nonZero_call <- nonZero_call[list,]

#read in the SNPs and subset to the ones still left in techRepCall
snp_file <- "../vcf/All_SNPs_Merged_mergedPositions.bed"
snps <- read.table(file = snp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-2]
colnames(snps) <- c("chr","pos","allele_A","allele_B")
rownames(snps) <- paste0(snps$chr,"-",snps$pos)
snps <- snps[rownames(nonZero_call),]

#convert NA to -1 for lookup purposes (cant lookup NA as a named vector)
nonZero_call[is.na(nonZero_call)] <- "-1"

#create data frame to hold genotype probabilites
genotype_prob <- as.data.frame(matrix(nrow = nrow(nonZero_call), ncol = ncol(nonZero_call), NA))
colnames(genotype_prob) <- colnames(nonZero_call)
rownames(genotype_prob) <- rownames(nonZero_call)

#create a named vector to allow facile conversion from birdseed scores to genotype probabilities
lookup <- c("-1" = "0|0:1:0.333,0.333,0.333", "0" = "0|0:0:1,0,0", "1" = "0|1:0.333:0,1,0", "2" = "1|1:2:0,0,1")

for (x in 1:ncol(genotype_prob)) {
  print(x)
  for (y in 1:nrow(genotype_prob)) {
    genotype_prob[y,x] <- lookup[as.character(nonZero_call[y,x])]
  }
}

#create final dataframe
demuxlet_df <- data.frame(CHROM = snps$chr,
                          POS = snps$pos,
                          ID = paste0(snps$chr,"-",snps$pos),
                          REF = snps$allele_A,
                          ALT = snps$allele_B,
                          QUAL = rep(".",nrow(genotype_prob)),
                          FILTER = rep("PASS",nrow(genotype_prob)),
                          INFO = rep(".",nrow(genotype_prob)),
                          FORMAT = rep("GT:DS:GP",nrow(genotype_prob))
)
demuxlet_df <- cbind(demuxlet_df, genotype_prob)

#remove lines that correspond to non standard chromosomes
demuxlet_df <- demuxlet_df[grep(pattern = "_", x = demuxlet_df$CHROM, invert = TRUE),]

#reorder snps to be in chromosome numeric order
chr_numbers <- gsub(pattern = "chr", replacement = "", x = demuxlet_df$CHROM)
chr_numbers <- gsub(pattern = "X", replacement = "99", x = chr_numbers)
demuxlet_df <- demuxlet_df[order(as.numeric(chr_numbers)),]

colnames(demuxlet_df)[1] <- paste0("#",colnames(demuxlet_df)[1])
write.table(x = demuxlet_df, file = paste0(out_dir,"/demuxlet_input.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

system2("cat", args = c(paste0(source_dir,"/","config_VCF_header.txt"), paste0(out_dir,"/demuxlet_input.txt")), stdout = paste0(out_dir,"/demuxlet_input.vcf"))


