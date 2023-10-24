
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

## Method for finding hyper-parameters
method <- "MCMC"
##--------------------------


###########################
## Find hyper-parameters ## -------
###########################
source("GPR_functions.R")

lhmax <- -10^10
for(knum in 1:length(vec.npar)){
  print("======================")
  print(paste("== Kernel number: ", knum, " ==", sep=""))
  print("======================")
  npar <- vec.npar[knum]
  if(method=="MCMC"){
    ## MCMC ##
    bounds <- t(matrix(rep(c(10^-7,10^4),npar),2, npar))
    hp_mcmc <- hpar.mcmc(xx,yy,bounds, knum)
    lh_mcmc <- f.lh(xx, yy, hp_mcmc, knum)
    print("Parameters:")
    print(hp_mcmc)
    print(paste("Likelihood(MCMC): ",lh_mcmc))
    if(lh_mcmc>lhmax){
      hpar <- hp_mcmc
      lhmax <- lh_mcmc
      kk <- knum
    }
  }
  
  if(method=="grad"){
    ## Gradient method ##
    hp_ini <- matrix(rep(c(0.01, 1, 100), npar)
                     ,npar,3)

    lh_grad <- -10^15
    for(i in 1:nrow(hp_ini)){
      print("============================")
      print(paste("Initial value case:",i))
      print("============================")
      hp0 <- hpar.grad(xx, yy, ini_th = hp_ini[i,], knum)
      lh0 <- f.lh(xx, yy, hp0, knum)
      print("Parameters:")
      print(hp0)
      print(paste("Likelihood(Grad): ",lh0))
      if(lh0 > lh_grad){
        hp_grad <- hp0
        lh_grad <- lh0
      }
      print("- - - - - - - - - - - - -")
    }
    if(lh_grad>lhmax){
      hpar <- hp_grad
      lhmax <- lh_grad
      kk <- knum
    }
  }
  
}

hpar
knum <- kk

## Save data 
fname <- "./res/hpar.txt"
write.table(c(knum,hpar),fname, row.names = F, col.names = F)
##---------------------------------



####################
## Reconstruction ## ---------------------
####################

res1 <- f.gpr(xx, yy, xx, hpar, knum)

hist(res1[,1]-yy)
er_rec <- sqrt(mean((res1[,1]-yy)^2)) # Error (RMSE)
print(paste("RMSE :", er_rec))

fname <- "./res/res_training.txt"
colnames(res1) <- c("Mean", "Sd")
write.table(res1,fname,row.names = F)

##----------------------------------------


##########################
## Result for test data ## ------
##########################

if(file.exists(fname.test)){
  res2 <- f.gpr(xx,yy,xv,hpar, knum)
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
  
  dx <- (range(xx)[2] - range(xx)[1])*0.05
  xp <- seq(range(xx)[1]-dx, range(xx)[2]+dx,length=1000)
  yp <- f.gpr(xx,yy,matrix(xp,length(xp),1),hpar, knum)
  
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
  res <- f.gpr(xx, yy, ss , hpar, knum)
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

