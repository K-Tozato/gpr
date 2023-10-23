#setwd("~/GitLab/stats/gpr/")


###############
## Functions ## ---------
###############

library(MASS)

f.ker1 <- function(x1, x2, th, f.diag=F, grad=F){
  kk <- th[1] * exp(-th[2]*sum((x1-x2)^2))
  dlt <- 0
  if(f.diag){
    kk <- kk + th[3]
    dlt <- 1
  }
  if(grad==T){
    kk <- c()
    kk[1] <- f.ker(x1,x2, th) - th[3]*dlt 
    kk[2] <- (f.ker(x1,x2,th) - th[3]*dlt) * (-th[2]*sum((x1-x2)^2))
    kk[3] <- th[3]*dlt
  }
  return(kk)
}

f.ker2 <- function(x1, x2, th, f.diag=F, grad=F){
  nn <- length(x1)
  kk <- th[1] * exp( -sum( th[2:(nn+1)]*((x1-x2)^2) ) )
  dlt <- 0
  if(f.diag){
    kk <- kk + th[nn+2]
    dlt <- 1
  }
  if(grad==T){
    kk <- c()
    kk[1] <- f.ker(x1, x2, th) - th[nn+2]*dlt 
    for(i in 1:nn){
      kk[i+1] <- (f.ker(x1,x2,th) - th[nn+2]*dlt) * (-th[i+1]*((x1[i]-x2[i])^2))
    }
    kk[nn+2] <- th[nn+2]*dlt
  }
  return(kk)
}

f.gpr <- function(x1,y1,x2,th){
  ymean <- mean(y1)
  y1 <- y1 - ymean
  
  kmat <- matrix(0,length(y1),length(y1))
  for(j in 1:length(y1)){
    for(i in 1:length(y1)){
      kmat[i,j] <- f.ker(x1[i,], x1[j,], th, i==j)
    }
  }
  Kinv <- solve(kmat,tol=10^-50)
  
  res <- matrix(0,nrow(x2),2)
  for(k in 1:nrow(x2)){
    k1 <- c()
    for(i in 1:length(y1)){
      k1[i] <- f.ker(x2[k,],x1[i,],th)
    }
    kk <- f.ker(x2[k,],x2[k,],th, T)
    
    res[k,] <- c(k1%*%Kinv%*%y1 + ymean, sqrt(max(0, kk - matrix(k1,1,length(k1))%*%Kinv%*%k1)))
  }
  
  return(res)
}

f.lh <- function(x1,y1,th){
  ymean <- mean(y1)
  y1 <- y1 - ymean
  
  kmat <- matrix(0,length(y1),length(y1))
  for(j in 1:length(y1)){
    for(i in 1:length(y1)){
      kmat[i,j] <- f.ker(x1[i,], x1[j,], th, i==j)
    }
  }
  Kinv <- solve(kmat,tol=10^-50)
  detK <- determinant(kmat)
  logdetK <- as.numeric(detK$modulus) * detK$sign
  
  if(detK$sign == 1){
    ll <- - 0.5*( logdetK + matrix(y1,1,length(y1))%*%Kinv%*%y1 + length(y1)*log(2*pi) )
  }
  if(detK$sign == -1){ll <- "Error"}
  
  return(ll)
}

hpar.mcmc <- function(x1, y1, bounds){
  ## (M-H (Metropolis-Hastings) method)
  
  n.itr <- 10000
  logb <- log(bounds)
  
  logp <- (logb[,1] + logb[,2])/2
  para <- exp(logp)
  loglh <- f.lh(x1,y1,para)
  
  par_list <- exp(logp)
  lh_list  <- loglh
  
  r0 <- 0.1 ; ii <- 0
  for(i in 1:n.itr){
    Sig <- matrix(0, npar, npar)
    diag(Sig) <- logb[,2]-logb[,1]
    move <- r0 * mvrnorm(1, rep(0,npar), Sig)
    
    fr <- c(min(logp + move - logb[,1]), max(logp + move - logb[,2]))
    if(fr[1] < 0 || fr[2] > 0){
      ii <- 0
      while(fr[1] < 0 || fr[2] > 0){
        ii <- ii + 1
        move <- (r0/ii) * mvrnorm(1, rep(0,npar), Sig)
        fr <- c(min(logp + move - logb[,1]), max(logp + move - logb[,2]))
        if(ii == 10^5){break}
      }
    }
    if(ii == 10^5){
      print("Error !! (ii > 10^5)")
      break
    }
    
    logp.n  <- logp + move
    para.n  <- exp(logp.n)
    loglh.n <- f.lh(x1,y1,para.n)
    
    if(loglh.n == "Error"){
      print("Error !! (det(K) <= 0)")
      break
    }
    
    rr <- exp(loglh.n - loglh)
    
    if(rr >= 1 || rr > runif(1,0,1)){
      para  <- para.n
      logp  <- logp.n
      loglh <- loglh.n
    }

    par_list <- rbind(par_list, exp(logp))
    lh_list[i]  <- loglh
    
    if(i%%1000==1){
      print("- - - - - - - - - - - - - - -")
      print( paste("Iteration",i) )
      print("Parameters: ")
      print(exp(logp))
      print(paste("Likelihood: ", max(lh_list)))
    }
  
  }

  nmx <- (1:length(lh_list))[lh_list==max(lh_list)]
  rr_mcmc <- par_list[max(nmx),]
  return(rr_mcmc)
}

hpar.grad <- function(x1, y1, ini_th=c()){
  ymean <- rep(0,length(y1))
  #ymean <- mean(y1)
  #y1 <- y1 - ymean
  
  if(length(ini_th)!=npar){bb <- rep(0, npar)}
  if(length(ini_th)==npar){bb <- log(ini_th)}
  eta  <- 0.005
  fer <- 0 ; dd0 <- 100
  for(it in 1:100000){
    th <- exp(bb) 
    
    kmat  <- matrix(0,length(y1),length(y1))
    dkmat <- matrix(0,length(y1)*length(y1), npar)
    ii <- 0
    for(j in 1:length(y1)){
      for(i in 1:length(y1)){
        ii <- ii + 1
        kmat[i,j]  <- f.ker(x1[i,],x1[j,], th, i==j)
        dkmat[ii,] <- f.ker(x1[i,],x1[j,], th, i==j, grad = T)
      }
    }
    Kinv <- solve(kmat,tol=10^-50)
    
    dfi <- rep(0,length=npar)
    for(i in 1:npar){
      K1 <- matrix(dkmat[,i],length(y1),length(y1))
      dfi[i] <- - 0.5*sum(diag(Kinv%*%K1)) + 0.5*matrix(y1,1,length(y1))%*%Kinv%*%K1%*%Kinv%*%y1
    }
    bb <- bb + eta * dfi

    if(it > 1){dd0 <- dd}
    dd <- sqrt(sum(dfi^2))
    aa <- exp(bb)
    
    if(it%%500==1){
      print("- - - - - - - - - - - - - - -")
      print( paste("Iteration",it) )
      print( paste("||df|| =", dd) )
      print("Parameters: ")
      print(aa)
      print(paste("Likelihood: ", f.lh(x1,y1,aa)))
    }
    
    if(dd0<dd){
      fer <- fer + 1
      if(fer >100){print("Error !!");break}
    }
    if(dd<dd0){fer <- 0}
    if(dd < 10^-3){break}
  }
  
  th <- exp(bb)
  return(th)
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

n.gpr <- 1 ## 1or2
if(n.gpr==1){f.ker <- f.ker1 ; npar <- 3}
if(n.gpr==2){f.ker <- f.ker2 ; npar <- 2+ncol(xx)}
##--------------------------


###########################
## Find hyper-parameters ## -------
###########################

## MCMC ##
bounds <- t(matrix(c(10^-7,10^4,10^-7,10^4,
                     10^-7,10^4),2,3))
hp_mcmc <- hpar.mcmc(xx,yy,bounds)
lh_mcmc <- f.lh(xx, yy, hp_mcmc[1:3])
print("Parameters:")
print(hp_mcmc)
print(paste("Likelihood(MCMC): ",lh_mcmc))


## Gradient method ##
hp_ini <- matrix(c(0.01,0.01,0.01,
                   10,1,0.1,
                   0.01,0.1,10)
                 ,npar,3)
hp_ini <- t(hp_ini)

lh_grad <- -10^15
for(i in 1:nrow(hp_ini)){
  print("============================")
  print(paste("Initial value case:",i))
  print("============================")
  hp0 <- hpar.grad(xx, yy, ini_th = hp_ini[i,])
  lh0 <- f.lh(xx, yy, hp0)
  print("Parameters:")
  print(hp0)
  print(paste("Likelihood(Grad): ",lh0))
  if(lh0 > lh_grad){
    hp_grad <- hp0
    lh_grad <- lh0
  }
  print("- - - - - - - - - - - - -")
}

if(lh_grad >= lh_mcmc){
  hpar <- hp_grad
}else{
  hpar <- hp_mcmc
}

## Save data 
fname <- "./res/hpar.txt"
write.table(hpar,fname, row.names = F, col.names = F)
##---------------------------------



####################
## Reconstruction ## ---------------------
####################

res1 <- f.gpr(xx,yy,xx,hpar)

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
  res2 <- f.gpr(xx,yy,xv,hpar)
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
  yp <- f.gpr(xx,yy,matrix(xp,length(xp),1),hpar)
  
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
  res <- f.gpr(xx, yy, ss , hpar)
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

