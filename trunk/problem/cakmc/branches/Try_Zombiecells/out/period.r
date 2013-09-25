require('TSA')
library("lattice")
library("grid")
jet.colors <- colorRampPalette(c("white", "red", "blue", "green"))

parainfo <- read.csv("control.csv", strip.white=TRUE)
L <- parainfo$VALUE[parainfo$PARAMETER=='L']
H <- parainfo$VALUE[parainfo$PARAMETER=='H']
.tend <- parainfo$VALUE[parainfo$PARAMETER=='tend']
.tpinc <- parainfo$VALUE[parainfo$PARAMETER=='tpinc']
.divide <- which(parainfo$PARAMETER=='useomp')
ompinfo <- parainfo[.divide:nrow(parainfo),]
parainfo <- parainfo[1:(.divide-1),]
H <- H + 2

N <- 1000
T <- seq(N)
for (i in seq(N) ){
    datafile <- sprintf("%s%05d%s", "m", i, ".dat")
    z <- matrix(scan(datafile, n=L*H, quiet=TRUE),
                L, H, byrow=TRUE)
    per <- periodogram(z[,401])
    T[i] <- 1/per$freq[which.max(per$spec)]
}
