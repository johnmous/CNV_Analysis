# Author: Ioannis Moustakas, i.moustakas@uva.nl 
# Title: Calculate and plot the Fold Change in coverage across two samples. Meant for copy number variation discovery
#       Input is a tab separated file with two columns: genomic position, coverage at this position 
#       OR gene ID and number of reads mapping on the gene (htseq-count was used to create this file)


############ $$$$$$$$$$$$$$$$$ ##################
# Plot the average coverage per non overlapping windows

rm(list = ls())
path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/coverageAnalysis/"
plotPath= "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Images/copyNumberVariationAnalysis/slidingWindow1000nt/"

# create a list of sample name to the corresponding sample number, per experiment

sampleIDsForCAZ256 = list(Broth_WT_D30 = 2,
                          Broth_CAZ256_D6 = 4,
                          Broth_CAZ256_D9 = 5,
                          Broth_CAZ256_D12 = 6,
                          Broth_CAZ256_D21 = 7,
                          Broth_CAZ256_WO_D15 = 8
)

sampleIDsForMPN256= list(Broth_WT_D30 = 2,
                        Broth_MPN256_D5 = 9,
                        Broth_MPN256_D10 = 10,
                        Broth_MPN256_D25 = 11,
                        Broth_MPN256_WO_D15 = 12
)

sampleIDsForPT512 = list(Broth_WT_D30 = 2,
                         Broth_PT512_D5 = 14,
                         Broth_PT512_D12 = 15,
                         Broth_PT512_D22 = 16,
                         Broth_PT512_WO_D15 = 17
)

sampleIDsForTBM32 = list(Broth_WT_D30 = 2,
                         Broth_TBM32_D6 = 18,
                         Broth_TBM32_D13 = 19,
                         Broth_TBM32_D30 = 20,
                         Broth_TBM32_WO_D15 = 21
)

sampleIDsForCIP128 = list(Broth_WT_D30 = 2,
                          Broth_CIP128_D6 = 22,
                          Broth_CIP128_D12 = 23,
                          Broth_CIP128_D21 = 24,
                          Broth_CIP128_WO_D15 = 25
)

sampleIDsForTBM1024 = list(Evans_WT_D30 = 3,
                           Evans_TBM1024_D14 = 26,
                           Evans_TBM1024_D30 = 27,
                           Evans_TBM1024_WO_D15 = 28
)

sampleIDsForCIP256 = list(WT=1,
                          Broth_WT_D30 = 2,
                          Evans_WT_D30 = 3,
                          Evans_TCIP256_D11 = 29,
                          Evans_CIP256_D19 = 30
)

# load Broth_WT_D30 that the rest of the samples will be compared against
sampleOne <- read.delim(paste0(path, "002_Cover.txt"), header=F)
covOne <- sampleOne$V2


drawPlots <- function(sampleNames, covOne) {
  for (name in names(sampleNames)){
    i = as.integer(sampleNames[[name]])
    window=1000
    sampleTwo <- read.delim(paste0(path, sprintf("%03d_Cover.txt", i)), header=F)
    covTwo <- sampleTwo$V2
    coverageDiff <- c()
    startPoint <- c()
    # average coverage in a window and calculate log(c1)-log(c2)
    for (start in seq(1, length(covOne), by=window)){
      diff <- mean(log2(covTwo[start:start+window]+1))- mean(log2(covOne[start:start+window]+1))
      coverageDiff <- c(coverageDiff, diff) 
      startPoint <- c(startPoint, start) 
    }
    filename=paste0(plotPath, "Broth_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(startPoint, coverageDiff, main=paste0("Broth_WT_D30 vs ", name))
    dev.off()
}
}

drawPlots(sampleIDsForCAZ256, covOne)
drawPlots(sampleIDsForMPN256, covOne)
drawPlots(sampleIDsForPT512, covOne)
drawPlots(sampleIDsForTBM32, covOne)
drawPlots(sampleIDsForCIP128, covOne)

# load Evans_WT_D30 that will be used for comparison
sampleOne <- read.delim(paste0(path, "003_Cover.txt"), header=F)
covOne <- sampleOne$V2

drawPlots <- function(sampleNames, covOne) {
  for (name in names(sampleNames)){
    i = as.integer(sampleNames[[name]])
    window=1000
    sampleTwo <- read.delim(paste0(path, sprintf("%03d_Cover.txt", i)), header=F)
    covTwo <- sampleTwo$V2
    coverageDiff <- c()
    startPoint <- c()
    # average coverage in a window and calculate log(c1)-log(c2)
    for (start in seq(1, length(covOne), by=window)){
      diff <- mean(log2(covTwo[start:start+window]+1))- mean(log2(covOne[start:start+window]+1))
      coverageDiff <- c(coverageDiff, diff) 
      startPoint <- c(startPoint, start) 
    }
    filename=paste0(plotPath, "Evan's-WT-D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(startPoint, coverageDiff, main=paste0("Evan's-WT-D30 vs ", name))
    dev.off()
  }
}

drawPlots(sampleIDsForTBM1024, covOne)
drawPlots(sampleIDsForCIP256, covOne)

######## Draw 004 vs 005
sampleOne <- read.delim(paste0(path, "004_Cover.txt"), header=F)
covOne <- sampleOne$V2

i=5
window=1000
sampleTwo <- read.delim(paste0(path, sprintf("%03d_Cover.txt", i)), header=F)
covTwo <- sampleTwo$V2
coverageDiff <- c()
startPoint <- c()
# average coverage in a window and calculate log(c1)-log(c2)
for (start in seq(1, length(covOne), by=window)){
  diff <- mean(log2(covTwo[start:start+window]+1))- mean(log2(covOne[start:start+window]+1))
  coverageDiff <- c(coverageDiff, diff) 
  startPoint <- c(startPoint, start) 
}
filename=paste0(plotPath, "Broth-CAZ 256-D6vsBroth-CAZ 256-D9.png")
png(filename, width = 1024, height = 768)
plot(startPoint, coverageDiff, main="Broth-CAZ 256-D6vsBroth-CAZ 256-D9")
dev.off()



############ $$$$$$$$$$$$$$$$$ ##################
# Plot the read count per gene using htseq-count table

rm(list = ls())
path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/htseq-count/"
plotPath= "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/htseq-count/Plots/"

sampleOne <- read.delim(paste0(path, "002_ReadCountPerGene.txt"), header=F)
# remove the last 5 lines from the df
control <- sampleOne[c(0:(nrow(sampleOne)-5)), ]
covControl <- control$V2

# sampleTwo <- read.delim(paste0(path, "007_ReadCountPerGene.txt"), header=F)
# sample <- sampleTwo[c(0:(nrow(sampleOne)-5)), ]
# covSample <- sample$V2

# diff <- log(covSample+1)-log(covControl+1)
# filename=paste0(plotPath, "WT-D30.vs.CAZ256-D21.png")
# png(filename, width = 1024, height = 768)
# plot(c(1:length(diff)), diff)
# dev.off()
# index <- which(diff>2)
# sample[index, ]
#hist(diff, breaks = 100)
# Prokka_01665

drawPlots <- function(sampleNames, covControl) {
  for (name in names(sampleNames)){
    i = as.integer(sampleNames[[name]])
    sampleTwo <- read.delim(paste0(path, sprintf("%03d_ReadCountPerGene.txt", i)), header=F)
    sample <- sampleTwo[c(0:(nrow(sampleOne)-5)), ]
    covSample <- sample$V2

    coverageDiff <- log(covSample+1)-log(covControl+1)
      
    filename=paste0(plotPath, "geneLevel_Broth_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(c(1:length(coverageDiff)), coverageDiff, main=paste0("geneLevel_Broth_WT_D30 vs ", name))
    dev.off()
  }
}

drawPlots(sampleIDsForCAZ256, covControl)
drawPlots(sampleIDsForMPN256, covControl)
drawPlots(sampleIDsForPT512, covControl)
drawPlots(sampleIDsForTBM32, covControl)
drawPlots(sampleIDsForCIP128, covControl)


