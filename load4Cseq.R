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
setwd("/Users/katarzynamurat/Documents/4Cseq/test-data")
save(fc,file = 'fcdata.rdata')
