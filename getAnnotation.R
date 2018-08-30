require("FourCSeq", quietly = TRUE)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

fcf <- get(load('$input'))

apId <- "FBgn0000099" #flybase gene id of the ap gene is "FBgn0000099"
apGene <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene, 
                filter=list(gene_id=apId))

apPromotor <- promoters(apGene, upstream = 500, downstream=100)

frags <- rowRanges(fcf)
if(length(frags) != nrow(results))
  stop("Number of rows is not the same for the fragment data and results table.")
ov <- findOverlaps(apPromotor, frags)                

results <- results[subjectHits(ov),1:6]                

write.table(results, file= '$output4', quote = FALSE,row.names = FALSE, sep = "\t")
