require("FourCSeq", quietly = TRUE)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

fcf <- getZScores(fc)

zScore <- assay(fcf, "zScore")

save(fcf,file = '$output1')
write.table(zScore, file= '$output2', quote = FALSE,row.names = FALSE, sep = "\t")

pdf('$output3')

hist(zScore[,"ap_MESO_68h_1"], breaks=100)

dev.off()

pdf('$output4')

qqnorm(zScore[,"ap_MESO_68h_1"], main="Normal Q-Q Plot - ap_MESO_68h_1")

abline(a=0, b=1)

dev.off()

fcf <- addPeaks(fcf, zScoreThresh=3, fdrThresh=0.01)

pdf('$output5')

plotFits(fcf[,1], main="")

dev.off()

#gene <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
pdf('$output6')

plotZScores(fcf[,c("ap_WE_68h_1", "ap_WE_68h_2")],
            txdb=gene)
dev.off()
