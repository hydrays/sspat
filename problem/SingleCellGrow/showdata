# Kernel Density Plot

par(mfrow=c(3,1))
mcell <- as.numeric(data7[8,3:ncol(data7[1,])])
mcell[is.na(mcell)]<-0
d1 <- density(mcell,from=0)
mcell <- as.numeric(data6[8,3:ncol(data6[1,])])
mcell[is.na(mcell)]<-0
d2 <- density(mcell,from=0)
mcell <- as.numeric(data5[8,3:ncol(data5[1,])])
mcell[is.na(mcell)]<-0
d3 <- density(mcell,from=0)

plot(d1,ylim=c(0, max(d1$y, d2$y, d3$y)))
lines(d2,col='green')
lines(d3, col='red')

mcell <- as.numeric(data4[8,3:ncol(data4[1,])])
mcell[is.na(mcell)]<-0
d1 <- density(mcell, from=0)
mcell <- as.numeric(data3[8,3:ncol(data3[1,])])
mcell[is.na(mcell)]<-0
d2 <- density(mcell, from=0)
mcell <- as.numeric(data2[8,3:ncol(data2[1,])])
mcell[is.na(mcell)]<-0
d3 <- density(mcell, from=0)

plot(d1,ylim=c(0, max(d1$y, d2$y, d3$y)))
lines(d2,col='green')
lines(d3, col='red')


mcell <- as.numeric(data1[9,4:ncol(data1[1,])])
mcell[is.na(mcell)]<-0
d1 <- density(mcell, from=0)
mcell <- as.numeric(data0[9,4:ncol(data1[1,])])
mcell[is.na(mcell)]<-0
d2 <- density(mcell, from=0)

plot(d1,ylim=c(0, max(d1$y, d2$y)))
lines(d2,col='green')

