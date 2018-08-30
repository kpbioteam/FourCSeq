require("FourCSeq", quietly = TRUE)

fc <- get(load('fcdata.rdata'))

#fc <- get(load('$input'))

fc <- addFragments(fc)

findViewpointFragments(fc)

fc <- addViewpointFrags(fc)

save(fc,file = 'addFragments.rdata')
