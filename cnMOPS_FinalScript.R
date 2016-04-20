# Author: Ioannis Moustakas, i.moustakas@uva.nl
# Title: CNV analysis for Yanfang data. Use cn.mops package. Final working version


########## $$$$$$$$$ ##########
# run it over all samples
library(cn.mops)
setwd("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/")

BamFiles <- list.files(pattern=".bam$")
# BamFiles <- c(BamFiles[2], BamFiles[4:8])
bamDatRanges <- getReadCountsFromBAM(BamFiles, mode = "unpaired", parallel=5)
# haplcn.mops is for haploid organisms 
res <- haplocn.mops(bamDatRanges)
res <- calcIntegerCopyNumbers(res)
plot(res, which=1)

#this saves the results in a dataframe
CNVRegions <- as.data.frame(cnvr(res))
# sample names are randomly thrown into the table. Put them in order (from 001 to 030)
sortedNames <- sort(colnames(CNVRegions)[6:35])
names <- colnames(CNVRegions)[1:5]
sortedNames <- c(names, sortedNames)
CNVRegions <- CNVRegions[sortedNames]

# Take experiment one (CAZ256) and remove lines with identical content (no variation between samples identified)
CAZ256 <- CNVRegions[ ,c(1:5, 7, 9:13)]
head(CAZ256)
lineBoolean <- apply(CAZ256, 1, function(x) {
                    print(x[6:11])
                    uniqCNVsInRow <- unique(x[6:11])
                    if (length(uniqCNVsInRow)==1) F
                    else T
})

CNVsCAZ256 <- CAZ256[lineBoolean, ]
write.table(CNVsCAZ256, file="CNVsCAZ256FromWhole.txt", sep="\t", row.names=F, quote = F)
# this approach results in more CNVs identified. Better to proceed with checking each experiment individually

write.table(CNVRegions, file="allExperiments.txt", sep="\t", row.names=F, quote = F)

# for (i in 1:nrow(CNVRegions)){
#   plot(res, which=i)
# )


######### $$$$$$$$ #########
# run it over each of the experiments individually 
library(cn.mops)
setwd("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/")

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

sampleIDsForCIP256 = list(Evans_WT_D30 = 3,
                          Evans_TCIP256_D11 = 29,
                          Evans_CIP256_D19 = 30
)


library(cn.mops)
setwd("/zfs/datastore0/group_root/MAD-RBAB/02_Collaborators/MAD1117-Ter_Kuile/MAD1117-P001-gDNA_sequencing/MAD1117-P001-E001_2014_32x_gDNA_Yanfang_svleeuw1/Results/outputBaseClearReference/bamFiles/")

CallCNVs <- function(sample, experimentName){
  sampleIndex <- unlist(sample)
  sampleNames <- names(sample)

  BamFiles <- list.files(pattern=".bam$")
  BamFiles <- BamFiles[sampleIndex]
  bamDatRanges <- getReadCountsFromBAM(BamFiles, mode = "unpaired", parallel=5)
  # haplcn.mops is for haploid organois 
  res <- haplocn.mops(bamDatRanges)
  res <- calcIntegerCopyNumbers(res)
  #this saves the results in a dataframe
  CNVRegions <- as.data.frame(cnvr(res))
  # sample names are randomly thrown into the table. Put them in order 
  #numberColumns <- ncol(CNVRegions)
  sortedNames <- paste0("X", BamFiles)
  names <- colnames(CNVRegions)[1:5]
  sortedNames <- c(names, sortedNames)
  CNVRegions <- CNVRegions[sortedNames]
  # Rename the sample names so the tables are easy to read
  colnames(CNVRegions) <- c(names, sampleNames)
  
  write.table(CNVRegions, file=paste0("CNVsWindow100_", experimentName, ".txt"), sep="\t", row.names=F, quote=F)
  
}

CallCNVs(sampleIDsForCAZ256, "CAZ256")
CallCNVs(sampleIDsForMPN256, "MPN256")
CallCNVs(sampleIDsForPT512, "PT512")
CallCNVs(sampleIDsForTBM32, "TMB32")
CallCNVs(sampleIDsForCIP128, "CIP128")
CallCNVs(sampleIDsForTBM1024, "TBM1024")
CallCNVs(sampleIDsForCIP256, "CIP256")


