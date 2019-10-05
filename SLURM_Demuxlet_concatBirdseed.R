################################################################################################################
#1/24/18
#Last Updated 11/27/18
#Ryan Corces
#SLURM_Demuxlet_concatBirdseed.R
#
#This program is meant to concatenate ATAC-seq birdseed-style files into a single matrix based on the "call" column.
#Any row which has exclusively "-1" calls will be removed from the downstream matrix to conserve space.
#The script reads in an input directory and processes all files containing an ".affx" suffix
#
#
#NOTES:
#usage - 
#example - 
#Rscript SLURM_Demuxlet_concatBirdseed.R --inDir --outFile
#where:
#--inDir </path/to/inputMpileupMatrix.txt>
#--outFile </path/to/outfile.txt>
#All arguments are required
#
#On Sherlock2 - 
#sbatch --mem=128000 --cpus-per-task=20 --time=24:00:00 --partition=howchang,sfgf --distribution=block --ntasks=1 --job-name=concatBirdseed  --output=/scratch/users/mcorces/temp/logs/concatBirdseed.log --wrap="Rscript /share/PI/howchang/users/mcorces/scripts/ATAC/TCGA_concatenateBirdseed.R --inDir /scratch/users/mcorces/TCGA_ATAC/genotyping/ATAC_birdseed/ --outFile /scratch/users/mcorces/TCGA_ATAC/genotyping/test.outfile.txt"
#
#output files are automatically named by changing the filename suffix from ".vcf" to ".snp" or ".affx" depending
#
#After running this script, you have a matrix of snps by samples. This should then be used as input to
#the script TCGA_correlateBirdseed.R which will determine the similarity between the ATAC data and TCGA SNP array data.
#
#
################################################################################################
## Load dependencies
library(optparse)
################################################################################################
#Global Variables and functions

################################################################################################
## Input variables

print(paste("Start",date(),sep=" - "))

## uses optparse package to read in command line arguments

option_list = list(
  make_option(c("--inDir"), action="store", default=NULL, type='character',
              help="/path/to/input/directory/"),
  make_option(c("--outFile"), action="store", default=NULL, type='character',
              help="/path/to/output/file.txt")
)

opt = parse_args(OptionParser(option_list=option_list))

if(is.null(opt$inDir) || is.null(opt$outFile)){
  # print usage information if an argument is not provided
  stop("\n Usage: Rscript /share/PI/howchang/users/mcorces/scripts/ATAC/TCGA_concatenateBirdseed.R --inDir </path/to/input/directory/> --outFile </path/to/output/file.txt> \n
       All arguments are required! \n\n")
}

inputDir <- opt$inDir
outputFile <- opt$outFile


print(paste("input dir = ",opt$inDir),sep="")
print(paste("output file = ",opt$outFile),sep="")

################################################################################################
setwd(inputDir)
#find all files to be concatenated (suffix = ".affx")
affxFiles <- list.files(inputDir, pattern = "*.affx")

#build matrix from first file
#read in first file
input <- read.table(affxFiles[1], header = TRUE, row.names = 1, check.names = TRUE, sep = "\t",
                    colClasses = c("character","character","character","integer","character","integer","character","integer","integer"))
birdseedMat <- as.data.frame(matrix(0, nrow = nrow(input), ncol = 1))
rownames(birdseedMat) <- input$probeID
birdseedMat$V1 <- input$call
colnames(birdseedMat)[1] <- unlist(strsplit(affxFiles[1], split = "\\."))[1]

#concatenate all other files
for (i in 2:length(affxFiles))
{
  print(paste("Processing file #",i," - ",affxFiles[i]))
  #read in file
  input <- read.table(affxFiles[i], header = TRUE, row.names = 1, check.names = TRUE, sep = "\t",
                      colClasses = c("character","character","character","integer","character","integer","character","integer","integer"))
  temp <- data.frame(row.names = input$probeID, call = input$call)
  birdseedMat <- cbind(birdseedMat, call = temp$call)
  colnames(birdseedMat)[i] <- unlist(strsplit(affxFiles[i], split = "\\."))[1] #take file prefix (before first ".") as the column name
}

#remove all rows where only "-1" calls are made to reduce file size
finalMat <- birdseedMat[-(which(rowSums(birdseedMat) == -(ncol(birdseedMat)))),]

write.table(finalMat, outputFile,sep="\t", col.names = NA, quote = FALSE)

print(paste("Start",date(),sep=" - "))
