require("FourCSeq", quietly = TRUE)

#fc <- get(load('$input'))

fc <- smoothCounts(fc)

save(fc,file = '$output1')
