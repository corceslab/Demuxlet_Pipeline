################################################################################################################
#10/05/19
#Last Updated 10/05/19
#Ryan Corces
#SLURM_Demuxlet_makeDemuxletVCF.R
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
#Rscript SLURM_Demuxlet_makeBirdseed.R --input </path/to/inputMpileupMatrix.txt> --affx </path/to/affxSNP6_metadata.txt> --outdir </path/to/outputDirectory>
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
library(matrixStats)
library(ArchRx)
################################################################################################
#Global Variables and functions
################################################################################################




# setwd("K:/Team Drives/TCGA-scATAC/Analysis/CellHashing/")
setwd("K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/")

# fileList <- list.files(pattern = ".*affx", path = "K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/affx/", full.names = TRUE)
fileList <- list.files(pattern = ".*affx", path = "K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/191003_affx_trimmed/", full.names = TRUE)


template <- read.table(fileList[1], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
numRows <- nrow(template)
numCols <- length(fileList)

depth <- as.data.frame(matrix(nrow = numRows, ncol = numCols, NA))
colnames(depth) <- gsub(pattern = ".pe.q10.sort.rmdup.mpileup.affx", replacement = "", x = basename(fileList))
rownames(depth) <- rownames(template)
call <- depth

for (i in 1:length(fileList)) {
  print(i)
  current <- read.table(fileList[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  current$allele_A_counts[which(is.na(current$allele_A_counts))] <- 0
  current$allele_B_counts[which(is.na(current$allele_B_counts))] <- 0
  
  depth[,i] <- current$allele_A_counts + current$allele_B_counts
  call[,i] <- current$call
}

saveRDS(object = depth, file = paste0(format(Sys.Date(), "%y%m%d"),"_SNP-depth.rds"))
saveRDS(object = call, file = paste0(format(Sys.Date(), "%y%m%d"),"_SNP-call.rds"))

nonZero_depth <- depth[which(rowSums(depth) > 0),]
nonZero_call <- call[which(rowSums(depth) > 0),]
saveRDS(object = nonZero_depth, file = paste0(format(Sys.Date(), "%y%m%d"),"_SNP-depth_nonZeroDepth.rds"))
saveRDS(object = nonZero_call, file = paste0(format(Sys.Date(), "%y%m%d"),"_SNP-call_nonZeroDepth.rds"))

#--------------------------------------------------------------------
# setwd("K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/")
#read output from Sherlock back in
# nonZero_depth <- readRDS("K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/190925_SNP-depth_nonZeroDepth.rds")
# nonZero_call <- readRDS("K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/190925_SNP-call_nonZeroDepth.rds")

#we dont have a techRep for PBMC_CellHash_Bulk_Donor6539 so lets just make a fake one that is the same as the single rep
temp <- data.frame( PBMC_CellHash_Bulk_Donor6539_RepA= nonZero_depth[,"PBMC_CellHash_Bulk_Donor6539_RepB"])
# nonZero_depth <- cbind(nonZero_depth[,1:30],temp,nonZero_depth[,31:39])
nonZero_depth <- cbind(nonZero_depth[,1:24],temp,nonZero_depth[,25:31])
temp2 <- data.frame( PBMC_CellHash_Bulk_Donor6539_RepA= nonZero_call[,"PBMC_CellHash_Bulk_Donor6539_RepB"])
# nonZero_call <- cbind(nonZero_call[,1:30],temp2,nonZero_call[,31:39])
nonZero_call <- cbind(nonZero_call[,1:24],temp2,nonZero_call[,25:31])
#this creates problems downstream because there are presumed to be at least one mismatch so we're just going to change a single value to avoid problems
nonZero_call[1,"PBMC_CellHash_Bulk_Donor6539_RepA"] <- 1

#Merge data from tech reps to only keep positions that are able to be called in both reps
#make matrices to hold data corresponding to tech rep comparisons
techRep_depth <- as.data.frame(matrix(nrow = nrow(nonZero_depth), ncol = ncol(nonZero_depth)/2, NA))
rownames(techRep_depth) <- rownames(nonZero_depth)
colnames(techRep_depth) <- unique(gsub(pattern = "_Rep.", replacement = "", x = colnames(nonZero_depth), perl = TRUE))

techRep_call <- as.data.frame(matrix(nrow = nrow(nonZero_depth), ncol = ncol(nonZero_depth)/2, NA))
rownames(techRep_call) <- rownames(nonZero_depth)
colnames(techRep_call) <- unique(gsub(pattern = "_Rep.", replacement = "", x = colnames(nonZero_depth), perl = TRUE))

#make dummy vectors to store depth of matches and mismatches
match_depth <- as.data.frame(matrix(nrow=0, ncol = 4, NA))
colnames(match_depth) <- c("id","depth","call","type")

#go through all tech rep pairs and check matches vs mismatches
for (i in 0:(length(techRep_depth) - 1)) {
  x <- (i * 2) + 1 #techRep 1 index
  y <- (x + 1) #techRep 2 index
  #subset the calls to those where both tech reps can make a call
  subset <- intersect(which(nonZero_call[,x] != -1), which(nonZero_call[,y] != -1))
  subset_depth <- nonZero_depth[subset,]
  subset_call <- nonZero_call[subset,]
  #identify indicies of matches and mismatches
  matches <- which(subset_call[,x] == subset_call[,y])
  mismatches <- which(subset_call[,x] != subset_call[,y])
  #update the corresponding cells in the techRep matrices with the values identified in the subsets for matches and mismatches
  #for matches, the call is (obviously) the matched call
  techRep_call[rownames(subset_depth)[matches],(i+1)] <- subset_call[matches,x]
  #for matches, the depth stored is the minimum of the depths of the two tech reps
  techRep_depth[rownames(subset_depth)[matches],(i+1)] <- rowMins(as.matrix(data.frame(A = subset_depth[matches,x], B = subset_depth[matches,y])))
  #for mismatches, the value of the call is made -1
  techRep_call[rownames(subset_depth)[mismatches],(i+1)] <- -1
  #for mismatches, the depth stored is the minimum of the depths of the two tech reps
  techRep_depth[rownames(subset_depth)[mismatches],(i+1)] <- rowMins(as.matrix(data.frame(A = subset_depth[mismatches,x], B = subset_depth[mismatches,y])))
  
  #for each techRep pair, add the match and mismatch data to the depth matrix
  
  match_bind <- data.frame(id = rownames(subset_depth)[matches],
                           depth = rowMins(as.matrix(data.frame(A = subset_depth[matches,x], B = subset_depth[matches,y]))),
                           call = subset_call[matches,x],
                           type = rep("match",length(matches)))
  mismatch_bind <- data.frame(id = rownames(subset_depth)[mismatches],
                              depth = rowMins(as.matrix(data.frame(A = subset_depth[mismatches,x], B = subset_depth[mismatches,y]))),
                              call = subset_call[mismatches,x],
                              type = rep("mismatch",length(mismatches)))
  
  match_depth <- rbind(match_depth, match_bind, mismatch_bind)
  
  print(paste0(x, ": ",as.character(round(length(matches) / (length(matches) + length(mismatches)) * 100, digits = 2))," matched."))
}

saveRDS(techRep_call, file = paste0(format(Sys.Date(), "%y%m%d"),"_techRepsMerged_GenotypeCall.rds"))
saveRDS(techRep_depth, file = paste0(format(Sys.Date(), "%y%m%d"),"_techRepsMerged_GenotypeDepth.rds"))

#preallocate a matrix to hold the number of differences between a given sample and all other samples
distance <- as.data.frame(matrix(nrow = ncol(techRep_call), ncol = ncol(techRep_call), 0))
rownames(distance) <- colnames(techRep_call)
colnames(distance) <- colnames(techRep_call)
#distance is a matrix where each cell represents the number of genotype call differences between the column and row
for (a in 1:nrow(distance)) {
  print(a)
  for(b in a:ncol(distance)) {
    #only one of the two following distance calculations should be uncommented
    #distance could be calculated based only on positions where BOTH samples can make a call
    # call_made <- intersect(which(techRep_call[,a] >= 0), which(techRep_call[,b] >= 0))
    #or distance could be calculated based on positions where EITHER sample can make a call
    call_made <- union(which(techRep_call[,a] >= 0), which(techRep_call[,b] >= 0))
    
    distance[a,b] <- length(which(techRep_call[call_made,a] != techRep_call[call_made,b]))
    distance[b,a] <- distance[a,b]
  }
}
saveRDS(distance, file = paste0(format(Sys.Date(), "%y%m%d"),"_techRepsMerged_GenotypeDistance_anyCall.rds"))
write.table(x = distance, file = paste0(format(Sys.Date(), "%y%m%d"),"_techRepsMerged_GenotypeDistance_anyCall.txt"), quote = FALSE, sep = "\t", col.names = NA)

#--------------------------------------------------------------------
#9/27/19
#Generate "genotype probabilities" for use with demuxlet
# setwd("K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/")
# techRep_call <- readRDS("techRepsMerged_GenotypeCall.rds")

#find all of the rows that have at least one genotyping call in one of the used samples
list <-c()
for (i in 1:nrow(techRep_call)) {
  if ((i %% 1000) == 0) {
    print(i)
  }
  if(length(which(is.na(techRep_call[i,]))) != ncol(techRep_call)) {
    list <- c(list,i)
  }
}
#subset techRep_call to only contain rows that had at least one sample with a genotyping call
techRep_call <- techRep_call[list,]

#read in the SNPs and subset to the ones still left in techRepCall
snp_file <- "K:/Shared drives/TCGA-scATAC/Analysis/CellHashing/PBMC_CellLines/All_SNPs_Merged_mergedPositions.bed"
snps <- read.table(file = snp_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,-2]
colnames(snps) <- c("chr","pos","allele_A","allele_B")
rownames(snps) <- paste0(snps$chr,"-",snps$pos)
snps <- snps[rownames(techRep_call),]

#convert NA to -1 for lookup purposes (cant lookup NA as a named vector)
techRep_call[is.na(techRep_call)] <- "-1"

techRep_prob <- as.data.frame(matrix(nrow = nrow(techRep_call), ncol = ncol(techRep_call), NA))
colnames(techRep_prob) <- colnames(techRep_call)
rownames(techRep_prob) <- rownames(techRep_call)




lookup <- c("-1" = "0|0:1:0.333,0.333,0.333", "0" = "0|0:0:1,0,0", "1" = "0|1:0.333:0,1,0", "2" = "1|1:2:0,0,1")

for (x in 1:ncol(techRep_prob)) {
  print(x)
  for (y in 1:nrow(techRep_prob)) {
    techRep_prob[y,x] <- lookup[as.character(techRep_call[y,x])]
  }
}

#create final dataframe
demuxlet_df <- data.frame(CHROM = snps$chr,
                          POS = snps$pos,
                          ID = paste0(snps$chr,"-",snps$pos),
                          REF = snps$allele_A,
                          ALT = snps$allele_B,
                          QUAL = rep(".",nrow(techRep_prob)),
                          FILTER = rep("PASS",nrow(techRep_prob)),
                          INFO = rep(".",nrow(techRep_prob)),
                          FORMAT = rep("GT:DS:GP",nrow(techRep_prob))
)
demuxlet_df <- cbind(demuxlet_df, techRep_prob)

#remove lines that correspond to non standard chromosomes
demuxlet_df <- demuxlet_df[grep(pattern = "_", x = demuxlet_df$CHROM, invert = TRUE),]

#reorder snps to be in chromosome numeric order
chr_numbers <- gsub(pattern = "chr", replacement = "", x = demuxlet_df$CHROM)
chr_numbers <- gsub(pattern = "X", replacement = "99", x = chr_numbers)
demuxlet_df <- demuxlet_df[order(as.numeric(chr_numbers)),]

demuxlet_pbmc <- demuxlet_df[,-grep(pattern = "CellLine", x = colnames(demuxlet_df))]
demuxlet_cellLine <- demuxlet_df[,-grep(pattern = "PBMC", x = colnames(demuxlet_df))]

write.table(x = demuxlet_df, file = paste0(format(Sys.Date(), "%y%m%d"),"_PBMC-CellLine_CellHash_demuxlet.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = demuxlet_cellLine, file = paste0(format(Sys.Date(), "%y%m%d"),"_PBMC-CellLine_CellHash_demuxlet_CellLine.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(x = demuxlet_pbmc, file = paste0(format(Sys.Date(), "%y%m%d"),"_PBMC-CellLine_CellHash_demuxlet_PBMC.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

