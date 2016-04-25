# Author: Ioannis Moustakas, i.moustakas@uva.nl 
# Title: Draw plots to facilitate the descovery of CNVs using per gene read count
#        Requires per gene read counts as calculated by HTSeq-count
rm(list = ls())

## Load all sample counts in a data frame

path = "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/htseq-count/"
countFiles <- list.files(pattern="_ReadCountPerGene.txt$", path = path)
samples <- read.delim(paste0(path, countFiles[1]), header=F)

for (file in countFiles[2:length(countFiles)]){
  sample <- read.delim(paste0(path, file), header=F)
  samples <- cbind(samples, sample[,2])
}
nameSplit <- strsplit(countFiles, "_ReadCountPerGene.txt", fixed=T)
sampleName <- sapply(nameSplit, function(x) x[1])
colnames(samples) <- c("GeneName", paste0("S", sampleName) )
# gene names as row names, then remove from DF
rownames(samples) <- samples$GeneName
samples <- samples[-1]
# print and shech the htseq summary lines
tail(samples)
# remove the last 5 lines of the htseq-count
samples <- samples[c(1:(nrow(samples)-5)), ]

# SF normalize the table
sizeFactors.mad <- function (counts, locfunc = median){
  loggeomeans <- rowMeans(log(counts))
  apply(counts, 2, function(cnts) exp(locfunc((log(cnts) -
                                                 loggeomeans)[is.finite(loggeomeans)])))
}
sf <- sizeFactors.mad(samples)
#divide countdata by sizefactors#
CountTable.scaled <- samples
for(i in 1:ncol(CountTable.scaled)){
  CountTable.scaled[,i] <- CountTable.scaled[,i]/sf[i]
} 


#### 
# Remove lines (read counts for asingle gene) averageing below 50
boolean <- apply(CountTable.scaled, 1, mean)>50 
CountTable.scaled <- CountTable.scaled[boolean, ]

sampleIDsForCAZ256 = list(WT= 1,
                          Broth_WT_D30 = 2,
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
                          Evans_CIP256_D11 = 29,
                          Evans_CIP256_D19 = 30
)

plotPath= "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/htseq-count/PlotsFiltered/"

# control Sample: Broth-WT-D30
covControl <- CountTable.scaled$S002

drawPlots <- function(sampleNames, covControl) {
  for (name in names(sampleNames)){
    i = as.integer(sampleNames[[name]])
    covSample <- CountTable.scaled[,i]
    
    coverageDiff <- log(covSample+1, 2)-log(covControl+1, 2)
    
    filename=paste0(plotPath, "geneLevelNormal_Broth_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(c(1:length(coverageDiff)), coverageDiff, main=paste0("geneLevel_Broth_WT_D30 vs ", name), 
         xlab = "Gene ID (Ordered accodring to genome coordinates)", ylab = "log2 Fold Change Coverage")
    dev.off()
    
    filename=paste0(plotPath, "geneLevelFixedScaleNormal_Broth_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(c(1:length(coverageDiff)), coverageDiff,  main=paste0("geneLevel_Broth_WT_D30 vs ", name), ylim=c(-3, 3), 
         xlab = "Gene ID (Ordered accodring to genome coordinates)", ylab = "log2 Fold Change Coverage")
    dev.off()
  }
}

drawPlots(sampleIDsForCAZ256, covControl)
drawPlots(sampleIDsForMPN256, covControl)
drawPlots(sampleIDsForPT512, covControl)
drawPlots(sampleIDsForTBM32, covControl)
drawPlots(sampleIDsForCIP128, covControl)


# save Prokka_01665 read counts in a table
write.table(samples[1665, ], file=paste0(path, "Prokka_01665_Counts.txt"), sep="\t", row.names = T, quote = F)

# Control sample Evan's-WT-D30
covControl <- CountTable.scaled$S003

drawPlots <- function(sampleNames, covControl) {
  for (name in names(sampleNames)){
    i = as.integer(sampleNames[[name]])
    covSample <- CountTable.scaled[,i]
    
    coverageDiff <- log(covSample+1, 2)-log(covControl+1, 2)
    
    filename=paste0(plotPath, "geneLevelNormal_Evans_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(c(1:length(coverageDiff)), coverageDiff, main=paste0("geneLevel_Evans_WT_D30 vs ", name), 
         xlab = "Gene ID (Ordered accodring to genome coordinates)", ylab = "log2 Fold Change Coverage")
    dev.off()
    
    filename=paste0(plotPath, "geneLevelFixedScaleNormal_Evans_WT_D30 vs ", name, ".png")
    png(filename, width = 1024, height = 768)
    plot(c(1:length(coverageDiff)), coverageDiff, main=paste0("geneLevel_Evans_WT_D30 vs ", name), ylim=c(-3, 3), 
         xlab = "Gene ID (Ordered accodring to genome coordinates)", ylab = "log2 Fold Change Coverage")
    dev.off()
  }
}

drawPlots(sampleIDsForTBM1024, covControl)
drawPlots(sampleIDsForCIP256, covControl)





# 
# geneName <- as.character(samples$GeneName)
# nameSplit <- strsplit(geneName, ".", fixed=T)
# geneNameNotransc <- sapply(nameSplit, function(x) x[1])
# samples$GeneName <- geneNameNotransc
# 
# write.table(samples, file=paste0(path, "allSamplesGeneCounts.txt"), sep="\t", row.names = F, quote = F)
# 
# 
# 
# 
# plotPath= "/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/htseq-count/Plots/"
