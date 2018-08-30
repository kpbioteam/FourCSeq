require("FourCSeq", quietly = TRUE)

fc <- get(load('fcdata.rdata'))
#exampleData (path)

#fc <- get(load('$input1'))

fc <- combineFragEnds(fc)

assays(fc)
head(assay(fc, "counts"))
data(fc)

metadata(fc)$projectPath <- "exampleData"
writeTrackFiles(fc)
writeTrackFiles(fc, format='bedGraph')

save(fc,file = '$output1')
write.table(writeTrackFiles(fc, format='bedGraph'), file= '$output2', quote = FALSE, sep = "\t")
