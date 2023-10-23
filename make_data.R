###############
## Make data ##
###############

fname <- "./data/daikatsu/all_rockfall_data.csv"
dat <- as.matrix(read.csv(fname, header = T))

ny <- 9
x <-  dat[(1:162)[c(-50,-100,-150)],c(3,4,5,6,7,8)]
y <- dat[(1:162)[c(-50,-100,-150)],ny]
xv <- dat[c(50,100,150),c(3,4,5,6,7,8)]
yv <- dat[c(50,100,150),ny]

fname <- "./data/training.txt"
write.table(cbind(x,y),fname,row.names = F, col.names = F)
fname <- "./data/test.txt"
write.table(cbind(xv,yv),fname,row.names = F, col.names = F)

#---------------------------------------


x <- c(seq(-30, -20, 0.25), seq(-19, -15, 2), 
       seq(-9, -5, 1), seq(-2.9, 1, 0.5), 
       seq(2, 10, 3), seq(11, 14, 3), 
       seq(14.1, 20, 0.5), seq(22, 30, 3))
y <- sin(x) / x + rnorm(length(x), sd=0.05)

nn <- seq(25,length(x),25)
xv <- x[nn] ; yv <- y[nn]
x <- x[-nn] ; y <- y[-nn]

fname <- "./data/training.txt"
write.table(cbind(x,y),fname,row.names = F, col.names = F)
fname <- "./data/test.txt"
write.table(cbind(xv,yv),fname,row.names = F, col.names = F)






