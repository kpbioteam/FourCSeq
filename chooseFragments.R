require("FourCSeq", quietly = TRUE)

#fc <- get(load('$input'))
colData(fc)$chr = "chr2L"
colData(fc)$start = 6027
colData(fc)$end = 6878

save(fc,file = '$output1')
