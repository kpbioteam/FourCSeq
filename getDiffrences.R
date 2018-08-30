require("FourCSeq", quietly = TRUE)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

###Detecting differences
fcf <- get(load('$input'))

fcf <- getDifferences(fcf,
                      referenceCondition="WE_68h")
pdf('$output1')

plotDispEsts(fcf)

dev.off()

pdf('$output2')

plotNormalizationFactors(fcf)

dev.off()

pdf('$output3')

plotMA(results(fcf, contrast=c("condition", "WE_68h", "MESO_68h")),
       alpha=0.01,
       xlab="Mean 4C signal",
       ylab="log2 fold change",
       ylim=c(-3.1,3.1))

dev.off()

results <- getAllResults(fcf)
write.table(results, file= '$output4', quote = FALSE,row.names = FALSE, sep = "\t")

pdf('$output5')

plotDifferences(fcf,
                txdb=TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                plotWindows = 1e+05,
                textsize=16)
dev.off()
