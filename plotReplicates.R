require("FourCSeq", quietly = TRUE)

#fc <- get(load('$input'))

pdf('$output1')

plotScatter(fc[,c("ap_WE_68h_1", "ap_WE_68h_2")],
            xlab="Replicate1", ylab="Replicate2", asp=1)

dev.off()
