setwd("~/GitLab/stats/gpr/")


###############
## Functions ## ---------
###############

library(MASS)

f.ker <- function(x1,x2,bt){
  y1 <- exp(-bt*sum((x1-x2)^2)) 
  return(y1)
}

f.gpr <- function(x1,y1,x2,bt,sig){
  ymean <- mean(y1)
  y1 <- y1 - ymean
  
  kmat <- matrix(0,length(y1),length(y1))
  for(j in 1:length(y1)){
    for(i in 1:length(y1)){
      kmat[i,j] <- f.ker(x1[i,],x1[j,],bt)
    }
  }
  kmat <- kmat + diag(sig^2,length(y1),length(y1))
  
  Kinv <- ginv(kmat)
  
  res <- matrix(0,nrow(x2),2)
  for(k in 1:nrow(x2)){
    k1 <- c()
    for(i in 1:length(y1)){
      k1[i] <- f.ker(x2[k,],x1[i,],bt)
    }
    kk <- f.ker(x2[k,],x2[k,],bt) + sig^2
    
    res[k,] <- c(k1%*%Kinv%*%y1 + ymean, sqrt(max(0, kk - k1%*%Kinv%*%k1)))
  }
  
  return(res)
}

f.cv <- function(x,y,bt,sig){
  fold <- 3
  
  nr <- fold ; nc <- ceiling(length(y)/fold)
  nn <- matrix(c(1:length(y),rep(0,length=nr*nc-length(y))),nr,nc)
  error <- c()
  for(kk in 1:nrow(nn)){
    n0 <- (nn[kk,])[nn[kk,]>0]
    x1 <- x[-n0,] ; x2 <- x[n0,]
    if(is.null(ncol(x1))){
      x1 <- matrix(x1,length(x1),1)
      x2 <- matrix(x2,length(x2),1)
    }
    y1 <- y[-n0] ; y2 <- y[n0]
    
    res <- f.gpr(x1,y1,x2,bt,sig)
    
    error <- c(error,res[,1] - y2)
  }
  return(sqrt(mean(error^2)))
}

decide.hpar <- function(x,y){
  r1 <- -7:4
  r2 <- -7:4
  jj <- 0
  er <- matrix(0,length(r1),length(r2))
  mn <- 10^10
  for(j in r2){
    jj <- jj + 1 ; ii <- 0
    for(i in r1){
      ii <- ii + 1
      er[ii,jj] <- f.cv(x,y,10^i,10^j)
      if(er[ii,jj]<mn){mn <- er[ii,jj] ; rr <- c(i,j)}
    }
  }
  min(er) ; c(10^rr[1], 10^rr[2])
  r1 <- seq(rr[1]-1,rr[1]+1,length=10) ; r2 <- seq(rr[2]-1,rr[2]+1,length=10)
  
  er <- 0 ; mn <- 10^10
  for(j in r2){
    for(i in r1){
      er <- f.cv(x,y,10^i,10^j)
      if(er<mn){mn <- er ; rr <- c(10^i,10^j)}
    }
  }
  mn ; rr
  
  return(rr)
}

##-----------------------


###############
## Read data ## ------------
###############

## Training data ##
fname <- "./data/training.txt"
dat <- as.matrix(read.table(fname,header = F))
xx <- dat[,1:(ncol(dat)-1)]
if(is.null(ncol(xx))){xx <- matrix(xx,length(xx),1)}
if(ncol(xx)>1){
  xra <- matrix(0,ncol(xx),2)
  for(i in 1:ncol(xx)){
    xra[i,] <- range(xx[,i])
    xx[,i]  <- (xx[,i] - xra[i,1])/(xra[i,2] - xra[i,1])
  }
}
yy <- dat[,ncol(dat)]


## Test data ##
fname.test <- "./data/test.txt"
if(file.exists(fname.test)){
  fname <- fname.test
  dat <- as.matrix(read.table(fname,header = F))
  xv <- dat[,1:(ncol(dat)-1)]
  if(is.null(ncol(xv))){xv <- matrix(xv,length(xv),1)}
  yv <- dat[,ncol(dat)]
  if(ncol(xx)>1){
    for(i in 1:ncol(xx)){
      xv[,i]  <- (xv[,i] - xra[i,1])/(xra[i,2] - xra[i,1])
    }
  }
}
##--------------------------


###########################
## Find hyper-parameters ## -------
###########################

rr <- decide.hpar(xx,yy)
bt <- rr[1] ; sig <- rr[2]

## Save data 
fname <- "./res/hpar.txt"
write.table(c(bt, sig),fname, row.names = F, col.names = F)
##---------------------------------


####################
## Reconstruction ## ---------------------
####################

res1 <- f.gpr(xx,yy,xx,bt,sig)

sqrt(mean((res1[,1]-yy)^2)) # Error (RMSE)

fname <- "./res/res_training.txt"
colnames(res1) <- c("Mean", "Sd")
write.table(res1,fname,row.names = F)

##----------------------------------------


##########################
## Result for test data ## ------
##########################

if(file.exists(fname.test)){
  res2 <- f.gpr(xx,yy,xv,bt,sig)
  cbind(res2[,1],yv)
  
  fname <- "./res/res_test.txt"
  colnames(res2) <- c("Mean", "Sd")
  write.table(res2,fname,row.names = F)
}

##-------------------------------



##########################################
## Compare GPR results with data (Plot) ## ------
##########################################

fname <- "./res/compare_output.png" ## figure name

ndata <- 2
type <- rep("p", ndata) # p, b or l ## Plot type 
lty <- rep(0,ndata) ## Line type
lwd <- rep(3,ndata) ## Line width
pch <- rep(16,ndata) ## Plot shape
clr <- 1:ndata ## Colors
cex <- rep(2.5,ndata) ## Plot size
font <- "serif" # serif, sans or mono
cex.lab <- 3.0 ## size of labels
cex.axis <- 2.8 ## size of axes

f.grid <- T # T or F
cex.leg <- 3
logxy <- "" # log axis ("x", "y" or "xy")

if(file.exists(fname.test)){
  xlim <- ylim <- range(yy,yv,res1[,1], res2[,1])
}
if(file.exists(fname.test)==F){
  xlim <- ylim <- range(yy,res1[,1])
}


p.ratio <- 1.0 ## plot ratio
f.legend <- T # T or F
legend <- c("Training data", "Test data") # legend 
p.leg <- "topleft" #Position for legend

xlab <- "Training (Test) data" # label for x axis
ylab <- "GPR results" # label for y axis


png(fname, width = 1000, height = p.ratio*1000)
par(mar = c(6,6,3,3), mgp = c(4,1.5,0))
par(family=font)

plot("", "", xlim= xlim, ylim = ylim, xlab = xlab, ylab = ylab,
     cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, log=logxy)
if(f.grid){grid()}
abline(0, 1, lwd=lwd, col=8)
par(new=T)
plot(yy, res1[,1], xlim= xlim, ylim = ylim, 
     xlab = "", ylab = "", 
     type=type[1], lty=lty[1], lwd=lwd[1], pch=pch[1], col=clr[1],
     cex = cex[1], xaxt="n", yaxt="n",
     log=logxy)
if(file.exists(fname.test)){
  par(new=T)
  plot(yv, res2[,1], xlim= xlim, ylim = ylim, 
       xlab = "", ylab = "", 
       type=type[2], lty=lty[2], lwd=lwd[2], pch=pch[2], col=clr[2],
       cex = cex[2], xaxt="n", yaxt="n",
       log=logxy)
}
if(f.legend){
  legend(p.leg,legend = legend, lty=lty, lwd=lwd, pch=pch, col=clr, cex=cex.leg)
}
dev.off()

##-----------------------------------------------


##########################################
## Graph ( only if ncol(xx) == 1 or 2 ) ##  -----
##########################################

if(ncol(xx)==1){
  
  fname <- "./res/output_graph.png" ## figure name
  
  png(fname, width = 1000, height = 800)
  par(mar = c(6,6,3,3), mgp = c(4,1.5,0))
  par(family=font)
  
  xp <- seq(range(xx)[1], range(xx)[2],length=100)
  yp <- f.gpr(xx,yy,matrix(xp,length(xp),1),bt,sig)
  
  xlim <- range(xx)
  ylim <- c(min(yp[,1]-2*yp[,2],yy),max(yp[,1]+2*yp[,2],yy))

  plot(xx, yy, col="blue", pch=19, xlim=xlim, ylim=ylim,cex.lab=cex.lab,
       cex.axis=cex.axis, xlab="x", ylab="y", cex=1.5)
  grid()
  lines(xp, yp[,1], col="red", lwd=2)
  lines(xp, yp[,1]  + 2 * yp[,2], col="salmon")
  lines(xp, yp[,1]  - 2 * yp[,2], col="salmon")
  polygon(c(xp, rev(xp)),
          c(yp[,1] + 2 * yp[,2], rev(yp[,1] - 2 * yp[,2])),
          col=adjustcolor("salmon", alpha.f=0.2), border=NA)
  dev.off()
}

if(ncol(xx)==2){
  
  s1 <- seq(min(xx[,1]), max(xx[,1]),length=50)
  s2 <- seq(min(xx[,2]), max(xx[,2]),length=50)
  ss <- matrix(0, length(s1)*length(s2), ncol(xx)) 
  ii <- 0
  for(i2 in 1:length(s2)){
    for(i1 in 1:length(s1)){
      ii <- ii + 1 ; ss[ii, ] <- c(s1[i1], s2[i2])
    }
  }
  res <- f.gpr(xx, yy, ss , hpar[1:3], hpar[4])
  res.m <- matrix(res[,1], length(s1), length(s2)) 
  res.v <- matrix(res[,2], length(s1), length(s2)) 
  
  rgl::persp3d(s1, s2, res.m, alpha=0.95, col=2)
  rgl::persp3d(s1, s2, res.m + 2*res.v, col="gray50", alpha=0.2, add=T)
  rgl::persp3d(s1, s2, res.m - 2*res.v, col="gray50", alpha=0.2, add=T)
  rgl::plot3d(xx[,1], xx[,2], yy, col=1, size=0.75, pch=20, type="s", add=T)
  rgl::aspect3d(1,1,1)
  rgl::axes3d()
  
  load("./data/pv_rs.txt")
  par3d(windowRect=pv$windowRect)  
  par3d(zoom=pv$zoom)              
  par3d(userMatrix=pv$userMatrix)  
  
  fname <- "./res/output_graph_3d.png"
  rgl.snapshot(fname, fmt="png")
  
}


##----------------------------------------

