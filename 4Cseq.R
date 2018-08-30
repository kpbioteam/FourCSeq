system.file("extdata/python/demultiplex.py", package="FourCSeq")
library(FourCSeq)
referenceGenomeFile = system.file("extdata/dm3_chr2L_1-6900.fa",
                                  package="FourCSeq")
referenceGenomeFile
## [1] "/tmp/RtmpZKRQPk/Rinst5d63188339ad/FourCSeq/extdata/dm3_chr2L_1-6900.fa"
bamFilePath = system.file("extdata/bam",
                          package="FourCSeq")
bamFilePath
## [1] "/tmp/RtmpZKRQPk/Rinst5d63188339ad/FourCSeq/extdata/bam"
primerFile = system.file("extdata/primer.fa",
                         package="FourCSeq")
primerFile
## [1] "/tmp/RtmpZKRQPk/Rinst5d63188339ad/FourCSeq/extdata/primer.fa"

writeLines(readLines(primerFile))

#Initialization of the FourC object

metadata <- list(projectPath = "exampleData",
                 fragmentDir = "re_fragments",
                 referenceGenomeFile = referenceGenomeFile,
                 reSequence1 = "GATC",
                 reSequence2 = "CATG",
                 primerFile = primerFile,
                 bamFilePath = bamFilePath)

colData <- DataFrame(viewpoint = "testdata",
                     condition = factor(rep(c("WE_68h", "MESO_68h", "WE_34h"),
                                            each=2),
                                        levels = c("WE_68h", "MESO_68h", "WE_34h")),
                     replicate = rep(c(1, 2),
                                     3),
                     bamFile = c("CRM_ap_ApME680_WE_6-8h_1_testdata.bam",
                                 "CRM_ap_ApME680_WE_6-8h_2_testdata.bam",
                                 "CRM_ap_ApME680_MESO_6-8h_1_testdata.bam",
                                 "CRM_ap_ApME680_MESO_6-8h_2_testdata.bam",
                                 "CRM_ap_ApME680_WE_3-4h_1_testdata.bam",
                                 "CRM_ap_ApME680_WE_3-4h_2_testdata.bam"),
                                  sequencingPrimer="first")
fc <- FourC(colData, metadata)

save(fc,file = 'fcdata.rdata')

#Fragment reference

fc <- addFragments(fc)
rowRanges(fc)
findViewpointFragments(fc)
fc <- addViewpointFrags(fc)

colData(fc)$chr = "chr2L"
colData(fc)$start = 6027
colData(fc)$end = 6878

fc <- countFragmentOverlaps(fc, trim=4, minMapq=30)
fc <- combineFragEnds(fc)
assays(fc)
head(assay(fc, "counts"))
data(fc)
metadata(fc)$projectPath
metadata(fc)$projectPath <- "exampleData"
writeTrackFiles(fc)
writeTrackFiles(fc, format='bedGraph')
fc <- smoothCounts(fc)

plotScatter(fc[,c("ap_WE_68h_1", "ap_WE_68h_2")],
            xlab="Replicate1", ylab="Replicate2", asp=1)

fcf <- getZScores(fc)

zScore <- assay(fcf, "zScore")

hist(zScore[,"ap_MESO_68h_1"], breaks=100)

#qqnorm(zScore[,"ap_MESO_68h_1"],
  #     main="Normal Q-Q Plot - ap_MESO_68h_1")
#abline(a=0, b=1)

fcf <- addPeaks(fcf, zScoreThresh=3, fdrThresh=0.01)
plotFits(fcf[,1], main="")
###########################################################
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
plotZScores(fcf[,c("ap_WE_68h_1", "ap_WE_68h_2")],
            txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene)
## [1] "ap"
## Successfully plotted results.

###Detecting differences

fcf <- getDifferences(fcf,
                      referenceCondition="WE_68h")
plotDispEsts(fcf)
plotNormalizationFactors(fcf)
plotMA(results(fcf, contrast=c("condition", "WE_68h", "MESO_68h")),
       alpha=0.01,
       xlab="Mean 4C signal",
       ylab="log2 fold change",
       ylim=c(-3.1,3.1))
results <- getAllResults(fcf)
dim(results)
## [1] 1872 16
head(results)[,1:6]
plotDifferences(fcf,
                txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                plotWindows = 1e+05,
                textsize=16)

apId <- "FBgn0000099" #flybase gene id of the ap gene is "FBgn0000099"
apGene <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene, 
                filter=list(gene_id=apId))

apPromotor <- promoters(apGene, upstream = 500, downstream=100)

frags <- rowRanges(fcf)
if(length(frags) != nrow(results))
  stop("Number of rows is not the same for the fragment data and results table.")
ov <- findOverlaps(apPromotor, frags)                
results[subjectHits(ov),1:6]                
