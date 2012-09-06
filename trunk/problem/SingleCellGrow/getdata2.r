data <- read.csv('dataA700.csv')
mcella1 <- as.numeric(data[8,3:ncol(data[8,])])
mcella1[is.na(mcella1)] <- 0

data <- read.csv('dataA5000.csv')
mcella2 <- as.numeric(data[8,3:ncol(data[8,])])
mcella2[is.na(mcella2)] <- 0

data <- read.csv('dataA130000.csv')
mcella3 <- as.numeric(data[8,3:ncol(data[8,])])
mcella3[is.na(mcella3)] <- 0

data <- read.csv('dataC3500.csv')
mcellc1 <- as.numeric(data[8,3:ncol(data[8,])])
mcellc1[is.na(mcellc1)] <- 0

data <- read.csv('dataC34000.csv')
mcellc2 <- as.numeric(data[8,3:ncol(data[8,])])
mcellc2[is.na(mcellc2)] <- 0

data <- read.csv('dataC368000.csv')
mcellc3 <- as.numeric(data[8,3:ncol(data[8,])])
mcellc3[is.na(mcellc3)] <- 0

x = c(mcella1, mcella2, mcella3, mcellc1, mcellc2, mcellc3)
write.table(x, file='cellsum3.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella1, mcella2, mcellc1, mcellc2)
write.table(x, file='cellsum2.csv', sep = ",", row.names = FALSE, col.names=FALSE)


x = c(mcella1, mcellc1)
write.table(x, file='cellsum1.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcellc1)
write.table(x, file='y.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella1, mcella2, mcella3)
write.table(x, file='cellsuma.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella1, mcella2)
write.table(x, file='cellsuma2.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcellc1, mcellc2, mcellc3)
write.table(x, file='cellsumc.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella1)
write.table(x, file='x1.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella2)
write.table(x, file='x2.csv', sep = ",", row.names = FALSE, col.names=FALSE)

x = c(mcella3)
write.table(x, file='x3.csv', sep = ",", row.names = FALSE, col.names=FALSE)