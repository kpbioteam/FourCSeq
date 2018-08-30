require("FourCSeq", quietly = TRUE)

fc <- get(load('fcdata.rdata'))

#fc <- get(load('$input1'))

fc <- countFragmentOverlaps(fc, trim=4, minMapq=30)

save(fc,file = '$output1')

