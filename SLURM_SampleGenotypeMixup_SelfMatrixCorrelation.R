################################################################################################################
#1/24/18
#Last Updated 11/27/18
#Ryan Corces
#SLURM_SampleGenotypeMixup_SelfMatrixCorrelation.R
#
#This program is meant to take a matrix of birdseed style correlations and perform an all by all correlation of
#each column with all other columns. This correlation is only performed on the positions where both samples being
#compared were able to make a birdseed call above the specified depth in upstream scripts
#
#
#NOTES:
#usage - 
#example - 
#Rscript /share/PI/howchang/users/mcorces/scripts/ATAC/SLURM_SampleGenotypeMixup_SelfMatrixCorrelation.R --inFile --tcgaDir --cores --outDir
#where:
#--inFile </path/to/inputMatrix.txt>
#--cores <numberOfCores>
#--outDir </path/to/outDir/>
#All arguments are required
#
#On Sherlock2 - 
#sbatch --mem=1500000 --cpus-per-task=32 --time=24:00:00 --partition=bigmem --distribution=block --ntasks=1 --job-name=correlateBirdseed  --output=/scratch/users/mcorces/temp/logs/correlateBirdseed.log --wrap="Rscript /share/PI/howchang/users/mcorces/scripts/ATAC/TCGA_correlateBirdseed.R --inFile </path/to/inputMatrix.txt> --cores <numberOfCores> --outDir </path/to/outDir.txt>"
#
#The matrix ouput of this script should be processed with TCGA_plotBirdseedCorrelations.R to plot the results.
#
################################################################################################
## Load dependencies
library(optparse)
library(foreach)
library(doParallel)
library(ggplot2)
################################################################################################
#Global Variables and functions

#Compares the pearson correlation values provided based on the metadata table provided
#expects a symmetrical pearson matrix where the number of rows/columns equals the number of rows provided in the metadata matrix
#expects the metadata matrix to have a column called "condition" which details the biological type of the sample
#expects the metadata matrix to have a column called "bioRepID" which is an identifier that is unique to the biological donor or group to which the samples belongs
compareCorrelations <- function(pearson, metadata)
{
  #check input to funciton
  #check for symmetric pearson matrix (colnames must be identical to rownames)
  try(if(!identical(colnames(pearson),rownames(pearson))) stop("Pearson matrix does not look symmetrical! Rownames do not equal Colnames."))
  #check for existence of "condition" and "bioRepID" as column names in metadata
  try(if(!all(c("condition","bioRepID") %in% colnames(metadata))) stop("Metadata matrix does not contain required columns - 'condition' and 'bioRepID'. Check metadata matrix."))
  
  #reorder the metadata matrix to have the same row order as the pearson matrix
  metadata <- metadata[rownames(pearson),]
  
  #refine pearson matrix to only samples that have a replicate as the columns but keep all rows
  #this allows outgroups to be fully representative (ie you dont exclude samples without a tech rep from contributing to the max out group correlation)
  #but this prevents samples without a tech rep from getting their own point on the plot
  pearsonTrimmed <- pearson[,which((metadata$bioRepID %in% metadata$bioRepID[which(duplicated(metadata$bioRepID))]))]
  #trim the bioRepID vector to match
  bioRepTrimmed <- metadata$bioRepID[which((metadata$bioRepID %in% metadata$bioRepID[which(duplicated(metadata$bioRepID))]))]
  #trim the metadata matrix to match
  metadataTrimmed <- metadata[which((metadata$bioRepID %in% metadata$bioRepID[which(duplicated(metadata$bioRepID))])),]
  
  #for each sample, check that the minimum pearson within the bioRep group is higher than the maximum pearson outside the bioRep group
  #make an empty matrix to hold the tech rep stats
  correlationStats <- as.data.frame(matrix(0,nrow = ncol(pearsonTrimmed), ncol = 5))
  colnames(correlationStats) <- c("inGroup","inGroupPearson","outGroup","outGroupPearson","difference")
  rownames(correlationStats) <- colnames(pearsonTrimmed)
  for (j in 1:nrow(correlationStats)) #for each sample that has a replicate
  {
    inGroup <- pearsonTrimmed[rownames(metadataTrimmed)[(which(metadataTrimmed$bioRepID == bioRepTrimmed[j]))],j, drop=FALSE]
    outGroup <- pearsonTrimmed[rownames(metadataTrimmed)[(which(metadataTrimmed$bioRepID != bioRepTrimmed[j]))],j, drop=FALSE]
    
    #determine the minimum correlation within the group of samples that share the same bioRepID
    minBioGroupCorrelation <- min(inGroup[,1])
    #determine which sample name corresponds to the min pearson in the bio rep group
    minBioGroupSample <- rownames(inGroup)[as.integer(which(inGroup[,1] == minBioGroupCorrelation))]
    #determine the maximum correlation outside of the bioRepID group
    maxOutGroupCorrelation <- max(outGroup[,1])
    #determine which sample name corresponds to the max pearson outside of the bioRepID group
    maxOutGroupSample <- rownames(outGroup)[as.integer(which(outGroup[,1] == maxOutGroupCorrelation))]
    
    #calculate the difference in pearson between the min of the bioRepGroup and the max of the outgroup
    difference <- minBioGroupCorrelation - maxOutGroupCorrelation
    #add values to matrix
    correlationStats[j,] <- c(minBioGroupSample,minBioGroupCorrelation,maxOutGroupSample,maxOutGroupCorrelation,difference)
  }
  
  #collect the sample type information from the metadata file
  correlationStats <- cbind(correlationStats, sampleGroups = metadataTrimmed$condition)
  
  #plot correlations
  #This will remove rows that correspond to samples that dont have Affy SNP 6 array data (currently 8)MatlabColorPanel = colorRampPalette(c("#3361A5", "#248AF3", "#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
  MatlabColorPanel = colorRampPalette(c("#3361A5", "#248AF3", "#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
  t1<-theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(size=.4),
    axis.text.x = element_text(face = "plain", color = "black", size = 12, hjust = 1, angle = 90, vjust = 0.5),
    axis.text.y = element_text(face = "plain", color = "black", size = 12, hjust = 1),
    axis.ticks = element_line(colour = 'black'),
    axis.title.x = element_text(face = "plain", color = "black", size = 14, hjust = 0.5),
    axis.title.y = element_text(face = "plain", color = "black", size = 14, hjust = 0.5),
    plot.title = element_text(face="bold", color = "black", size=12, hjust=0.5)
  )
  
  #create rounded numbers (to 1 decimal place) for plot extremities
  diff_min <- min(as.numeric(correlationStats$difference))
  diff_max <- max(as.numeric(correlationStats$difference))
  plot_min <- min(-0.1, sign(diff_min) * ceiling(abs(diff_min)*10)/10)
  plot_max <- sign(diff_max) * ceiling(abs(diff_max)*10)/10
  
  gg <- ggplot(correlationStats) + geom_jitter(aes(x=sampleGroups,y=as.numeric(difference), colour=as.numeric(inGroupPearson)), width=0.2) +
    t1  +  scale_y_continuous(limits = c(plot_min, plot_max), breaks=seq(from=plot_min, to=plot_max, by=0.1), labels=seq(from=plot_min, to=plot_max, by=0.1)) +
    ylab("Correlation Difference\n(Expected Donor - Next Closest)") + xlab("Sample Type") + scale_colour_gradientn(colours = MatlabColorPanel(256), limits = (c(0,1))) +
    geom_hline(yintercept = 0, color = "red2")
  
  ggsave(paste0(outputDir,"/AllSNPs_Genotype_Correlations.pdf"), gg, device="pdf", width=15, height=6)
  write.table(correlationStats, paste0(outputDir,"/AllSNPs_Genotype_Correlations_STATS.txt"), quote=FALSE, sep = "\t", col.names=NA)
}

################################################################################################
## Input variables

## uses optparse package to read in command line arguments

option_list = list(
  make_option(c("--inFile"), action="store", default=NA, type='character', help="/path/to/inputMatrix.txt"),
  make_option(c("--cores"), action="store", default=1, type='integer', help="number of cores to run if parallel processing is desired"),
  make_option(c("--outDir"), action="store", default=getwd(), type='character', help="/path/to/output/dir/"),
  make_option(c("--manifest"), action="store", default=NA, type='character', help="/path/to/manifest/file.txt")
)

opt = parse_args(OptionParser(option_list=option_list))

if(is.na(opt$inFile) || is.na(opt$manifest)){
  # print usage information if an argument is not provided
  stop("\n Usage: Rscript /share/PI/howchang/users/mcorces/scripts/ATAC/TCGA_correlateBirdseed.R --inFile </path/to/inputMatrix.txt> --tcgaDir </path/to/tcga/birdseed/rds/dir/> --cores <numberOfCores> --outDir </path/to/outDir/> --manifest </path/to/manifest.txt> \n
       All arguments are required except --cores! outDir defaults to current working directory if not supplied. \n\n")
} 

print(paste("input matrix = ",opt$inFile),sep="")
print(paste("cores = ",opt$cores),sep="")
print(paste("output file = ",opt$outDir),sep="")
print(paste("manifest file = ",opt$manifest),sep="")

inputFile <- opt$inFile
cores <- opt$cores
outputDir <- opt$outDir
manifestFile <- opt$manifest

#read in input file (matrix of ATAC birdseed calls)
atacBirdseed <- read.table(inputFile, header=TRUE, row.names=1, sep="\t", as.is=TRUE, check.names=FALSE)

#make empty vector to hold correlations
correlations <- as.data.frame(matrix(0,nrow=ncol(atacBirdseed), ncol=ncol(atacBirdseed)))
colnames(correlations) <- colnames(atacBirdseed)
rownames(correlations) <- colnames(atacBirdseed)

#establish parallel environment for correlations with TCGA RDS files
registerDoParallel(cores = cores)
c1 <- makeCluster(cores, type="FORK")


#for each ATAC sample (each column in atacBirdseed)
for (i in 1:ncol(atacBirdseed))
{
  print(i)
  #find the snp probes that actually have a call (ie not -1)
  atac <- atacBirdseed[which(atacBirdseed[,i] != -1),i,drop=FALSE]
  probes <- rownames(atac)
  #run parallel correlations with each individual TCGA RDS file and combine them into a column with rbind
  currentCor <- foreach(j = 1:ncol(atacBirdseed), .combine=rbind) %dopar%
  {
    #read in TCGA file
    tcga <- atacBirdseed[,j, drop = FALSE]
    tcga <- tcga[probes,, drop = FALSE]
    #subset TCGA file to only contain the desired probes
    cor(x = as.integer(atac[which(tcga != -1),1]), y = as.integer(tcga[which(tcga != -1),1]))
  }
  correlations[,i] <- currentCor[,1]
}

stopCluster(c1)
registerDoSEQ()

write.table(correlations, paste0(outputDir,"/AllSNPs_Genotype_Correlations.txt"),sep="\t", col.names = NA, quote = FALSE)

#run correlation difference calculations and make plot and stats file
manifest_df <- read.table(manifestFile, header=FALSE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
colnames(manifest_df) <- c("BAM","condition","bioRepID")
compareCorrelations(correlations, manifest_df)

