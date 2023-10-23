setwd("~/GitLab/stats/gpr/")
## Example 1 (1D-GPR) ## -----------------------------

xx <- c(seq(-3,-1,0.25), 0, seq(0.1,1.5,0.05), 2.4,2.6,2.8,3)
xx <- matrix(xx,length(xx),1)
yy <- sin(2*xx) + 2/3*xx + 1/10*(xx^2) + rnorm(length(xx),0,0.5^2) 
plot(xx,yy)

## MCMC ##
bounds <- t(matrix(c(10^-7,10^4,10^-7,10^4,
                     10^-7,10^4),2,3))
hp_mcmc <- hpar.mcmc(xx,yy,bounds)
lh_mcmc <- f.lh(xx, yy, hp_mcmc)
print("Parameters:")
print(hp_mcmc)
print(paste("Likelihood(MCMC): ",lh_mcmc))
hpar <- hp_mcmc[nrow(hpar),1:3]

hpar <- c(1, 1, 0.5)


fname <- "./res/output_graph_ppt.png" ## figure name

png(fname, width = 1000, height = 800)
par(mar = c(6,6,3,3), mgp = c(4,1.5,0))
par(family="serif")

dx <- (range(xx)[2] - range(xx)[1])*0.05
xp <- seq(range(xx)[1]-dx, range(xx)[2]+dx,length=1000)
yp <- f.gpr(xx,yy,matrix(xp,length(xp),1),hpar)

xth <- seq(range(xx)[1]-dx, range(xx)[2]+dx,length=1000)
yth <- sin(2*xth) + 2/3*xth + 1/10*(xth^2) 

xlim <- range(xx)
ylim <- c(min(yp[,1]-2*yp[,2],yy),max(yp[,1]+2*yp[,2],yy))

plot(xx, yy, col="blue", pch=19, xlim=xlim, ylim=ylim,cex.lab=3,
     cex.axis=3, xlab="x", ylab="y", cex=2)
grid()
lines(xp, yp[,1]  + 2 * yp[,2], col="salmon")
lines(xp, yp[,1]  - 2 * yp[,2], col="salmon")
polygon(c(xp, rev(xp)),
        c(yp[,1] + 2 * yp[,2], rev(yp[,1] - 2 * yp[,2])),
        col=adjustcolor("salmon", alpha.f=0.2), border=NA)
lines(xth, yth, col=1, lwd=2)
lines(xp, yp[,1], col="red", lwd=2.5)
dev.off()

## Example 2 (2D-GPR) ## ----------------------------

f.ref <- function(x){
  return(sin(2*x[1]) + 2/3*x[1] + 0.1*x[2]^2)
}

x1 <- c(seq(-3,0,0.5),1,seq(1.2,2.6,0.2),3)
x2 <- c(seq(-3,-1,0.25),0,seq(0.2,2.2,0.4),3)
ii <- 0
xx <- matrix(0,length(x1)*length(x2),2)
for(j in 1:length(x2)){
  for(i in 1:length(x1)){
    ii <- ii + 1
    xx[ii,] <- c(x1[i], x2[j])
  }
}
x1 <- seq(23,267, 4)
xx <- xx[-x1,]
xx <- xx[-(c(70:110, 170:195)),]
xx <- xx[-c(10:29,seq(50,130,4)),]
yy <- rep(0,nrow(xx))
for(i in 1:nrow(xx)){
  yy[i] <- f.ref(xx[i,])
}
yy <- matrix(yy+ rnorm(length(yy),0,0.5),length(yy),1) 

hpar <- c(0.1, 0.1, 0.05)


##PLOT
dx <- c(range(xx[,1])[2] - range(xx[,1])[1],
        range(xx[,2])[2] - range(xx[,2])[1])*0.05
x2 <- seq(range(xx[,1])[1]-dx[1], range(xx[,1])[2]+dx[1],length=60)
x1 <- seq(range(xx[,2])[1]-dx[2], range(xx[,2])[2]+dx[2],length=60)
xp <- matrix(0,length(x1)*length(x2),2)
yth <- matrix(0,length(x1),length(x2))
ii <- 0
for(j in 1:length(x2)){
  for(i in 1:length(x1)){
    ii <- ii + 1
    xp[ii,]  <- c(x1[i],x2[j])
    yth[i,j] <- f.ref(c(x1[i],x2[j]))
  }
}
yp <- f.gpr(xx,yy,xp,hpar)
yp.m <- matrix(yp[,1],length(x1),length(x2)) 
yp.s <- matrix(yp[,2],length(x1),length(x2))

rgl::persp3d(x1, x2, yth,  alpha=0.65, col=1)
rgl::persp3d(x1, x2, yp.m, alpha=0.95, col=2, add=T)
rgl::persp3d(x1, x2, yp.m + 2*yp.s, col="pink", alpha=0.4, add=T)
rgl::persp3d(x1, x2, yp.m - 2*yp.s, col="pink", alpha=0.4, add=T)
rgl::plot3d(xx[,1], xx[,2], yy , col=1, size=0.75, pch=20, type="s", add=T)
rgl::aspect3d(1,1,1)
rgl::axes3d()

## Example 3 (ARD with MCS) ## -----------

x1 <- rnorm(100, 0, 1)
x2 <- x1 + rnorm(length(x1), 0, 0.025)
x3 <- rnorm(100, 0, 1)
  
f.ref <- function(x){sin(2*pi*x)}

yy <- matrix(0,length(x1),1)
for(ii in 1:length(x1)){
  yy[ii]  <- f.ref(x1[ii]) + rnorm(1,0,0.2)
}
ymean <- mean(yy)
xx <- cbind(x1,x2,x3)

n.gpr <- 2 ## 1or2
if(n.gpr==1){f.ker <- f.ker1 ; npar <- 3}
if(n.gpr==2){f.ker <- f.ker2 ; npar <- 2+ncol(xx)}

hpar.mcmc <- function(x1, y1, bounds){
  ## (M-H (Metropolis-Hastings) method)
  
  n.itr <- 50000
  logb <- log(bounds)
  
  logp <- (logb[,1] + logb[,2])/2
  para <- exp(logp)
  loglh <- f.lh(x1,y1,para)
  
  par_list <- exp(logp)
  lh_list  <- loglh
  
  r0 <- 0.05 ; ii <- 0
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
    
    par_list <- rbind(par_list, para)
    lh_list[i+1]  <- loglh
    nmx <- (1:length(lh_list))[lh_list==max(lh_list)]
    nmx <- nmx[1]
    par_list[i+1,] <- par_list[nmx,]
    lh_list[i+1]  <- lh_list[nmx]
    
    if(i%%1000==1){
      print("- - - - - - - - - - - - - - -")
      print( paste("Iteration",i) )
      print("Parameters: ")
      print(exp(logp))
      print(paste("Likelihood: ", max(lh_list)))
    }
  }
  
  nmx <- (1:length(lh_list))[lh_list==max(lh_list)]
  nmx <- nmx[1]
  rr_mcmc <- exp(par_list[max(nmx),])
  return(cbind(par_list,lh_list))
}

## MCMC ##
bounds <- t(matrix(c(10^-7,10^3,10^-7,10^3,
                     10^-7,10^3,10^-7,10^3,
                     10^-7,10^3),2,5))
res <- hpar.mcmc(xx,yy,bounds)
hp_mcmc <- res[nrow(res),1:npar]
lh_mcmc <- f.lh(xx, yy, hp_mcmc)
print("Parameters:")
print(hp_mcmc)
print(paste("Likelihood(MCMC): ",lh_mcmc))



plot(res[,5], log="y")

fname <- "./res/ex3/contribution_ex3.png"
png(fname, width = 800, height = 600)
par(mar = c(6,6,3,6), mgp = c(4,1.5,0))
par(family="serif")

plot(res[,2], xlim=c(0,nrow(res)), ylim=range(res[,2:4]), log="y", col=1,
     lwd=2, type="l", cex.lab=2.5, cex.axis=2.5, xlab="Iteration", ylab="eta_i")
par(new=T)
plot(res[,3], xlim=c(1,nrow(res)), ylim=range(res[,2:4]), log="y", col=2,
     lwd=2, type="l", cex.lab=2.5, cex.axis=2.5, xlab="", ylab="")
par(new=T)
plot(res[,4], xlim=c(1,nrow(res)), ylim=range(res[,2:4]), log="y", col=3, 
     lwd=2, type="l", cex.lab=2.5, cex.axis=2.5, xlab="", ylab="")
grid()
dev.off()


(cr <- sqrt(hpar[2:3])/sum(sqrt(hpar[2:3])))

## calculate cr by MCS 
n <- 100000
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)
#plot(x1,x2)

xp <- cbind(x1,x2)
yp <- f.gpr(xx,yy,xp,hpar)
sig.all <- sd(yp[,1])

xp <- cbind(x1,rep(0,length(x1)))
yp <- f.gpr(xx,yy,xp,hpar)
sig1 <- sd(yp[,1])

xp <- cbind(rep(0,length(x1)),x2)
yp <- f.gpr(xx,yy,xp,hpar)
sig2 <- sd(yp[,1])

cr <- c(sig1,sig2) / sig.all
(cr <- cr / sum(cr))

## Example 4 (1D GPR + adding training data) ##-------------

f.ref <- function(xx){
  return(sin(2*xx) + 2/3*xx + 1/10*(xx^2))
}

# Initial training data #
xx <- c(seq(-3,-2.5,0.1), -1.5, seq(0.5,1.0,0.1), seq(2.5,3,0.05))
xx <- matrix(xx,length(xx),1)
yy <- f.ref(xx) + rnorm(length(xx),0,0.5^2) 
plot(xx,yy)


xlim <- range(xx)
ylim <- c(-3,4)


for(n.add in 1:10){
  print("=======================")
  print(paste("n.add:",n.add))
  print("=======================")
  
  bounds <- t(matrix(c(10^-7,10^4,10^-7,10^4,
                       10^-7,10^4),2,3))
  hp_mcmc <- hpar.mcmc(xx,yy,bounds)
  lh_mcmc <- f.lh(xx, yy, hp_mcmc[1:3])
  print("Parameters:")
  print(hp_mcmc)
  print(paste("Likelihood(MCMC): ",lh_mcmc))
  
  hpar <- hp_mcmc
  
  
  dx <- (range(xx)[2] - range(xx)[1])*0.05
  xp <- seq(range(xx)[1]-dx, range(xx)[2]+dx,length=600)
  yp <- f.gpr(xx,yy,matrix(xp,length(xp),1),hpar)
  
  xp1 <- seq(range(xx)[1], range(xx)[2],length=1000)
  yp1 <- f.gpr(xx,yy,matrix(xp1,length(xp1),1),hpar)
  np <- (1:nrow(yp1))[yp1[,2]==max(yp1[,2])]
  flag <- length((1:length(xx))[xp1[np]==xx])
  print(paste("FLAG:", flag))
  ii <- 0
  while(flag==1){
    ii <- ii + 1
    if(ii==1){np1 <- np}
    if(ii> 1){np1 <- c(np1,np)}
    np <- (1:nrow(yp1[-np1,]))[yp1[-np1,2]==max(yp1[-np1,2])]
    flag <- length((1:length(xx))[xp1[np]==xx])
  }
  xnew <- xp1[np] ; ynew <- f.ref(xnew) + rnorm(1, 0, 0.5^2)
  
  xth <- seq(range(xx)[1]-dx, range(xx)[2]+dx,length=1000)
  yth <- f.ref(xth) 
  
  #xlim <- range(xx)
  #ylim <- c(min(yp[,1]-2*yp[,2],yy),max(yp[,1]+2*yp[,2],yy))
  
  if(n.add==1){
    fname <- paste("./res/ex4/output_add_","0-1",".png",sep="") ## figure name
    png(fname, width = 1000, height = 800)
    par(mar = c(6,6,3,6), mgp = c(4,1.5,0))
    par(family="serif")
    
    plot(xx, yy, col="blue", main=paste("Step: ", 0), 
         pch=19, xlim=xlim, ylim=ylim,cex.lab=2.5,
         cex.axis=2.5, xlab="x", ylab="y", cex=2, cex.main = 3)
    grid()
    lines(xp, yp[,1]  + 2 * yp[,2], col="salmon")
    lines(xp, yp[,1]  - 2 * yp[,2], col="salmon")
    polygon(c(xp, rev(xp)),
            c(yp[,1] + 2 * yp[,2], rev(yp[,1] - 2 * yp[,2])),
            col=adjustcolor("salmon", alpha.f=0.2), border=NA)
    lines(xth, yth, col=1, lwd=2)
    lines(xp, yp[,1], col="red", lwd=2.5)
    dev.off()
    fname <- paste("./res/ex4/output_add_","0-2",".png",sep="") ## figure name
    png(fname, width = 1000, height = 800)
    par(mar = c(6,6,3,6), mgp = c(4,1.5,0))
    par(family="serif")
    
    plot(xx, yy, col="blue", main=paste("Step: ", 0), 
         pch=19, xlim=xlim, ylim=ylim,cex.lab=2.5,
         cex.axis=2.5, xlab="x", ylab="y", cex=2, cex.main = 3)
    grid()
    lines(xp, yp[,1]  + 2 * yp[,2], col="salmon")
    lines(xp, yp[,1]  - 2 * yp[,2], col="salmon")
    polygon(c(xp, rev(xp)),
            c(yp[,1] + 2 * yp[,2], rev(yp[,1] - 2 * yp[,2])),
            col=adjustcolor("salmon", alpha.f=0.2), border=NA)
    lines(xth, yth, col=1, lwd=2)
    lines(xp, yp[,1], col="red", lwd=2.5)
    res <- f.gpr(xx, yy, matrix(xnew,1,1), hpar)
    abline(v=xnew, col="green", lwd=3)
    dev.off()
  }
  
  fname <- paste("./res/ex4/output_add_",n.add,".png",sep="") ## figure name
  
  png(fname, width = 1000, height = 800)
  par(mar = c(6,6,3,6), mgp = c(4,1.5,0))
  par(family="serif")
  
  plot(xx, yy, col="blue", main=paste("Step: ", n.add), 
       pch=19, xlim=xlim, ylim=ylim,cex.lab=2.5,
       cex.axis=2.5, xlab="x", ylab="y", cex=2, cex.main = 3)
  grid()
  lines(xp, yp[,1]  + 2 * yp[,2], col="salmon")
  lines(xp, yp[,1]  - 2 * yp[,2], col="salmon")
  polygon(c(xp, rev(xp)),
          c(yp[,1] + 2 * yp[,2], rev(yp[,1] - 2 * yp[,2])),
          col=adjustcolor("salmon", alpha.f=0.2), border=NA)
  lines(xth, yth, col=1, lwd=2)
  lines(xp, yp[,1], col="red", lwd=2.5)
  par(new=T)
  plot(xnew, ynew, col="green", 
       pch=19, xlim=xlim, ylim=ylim,
       cex.axis=2.5, xlab="", ylab="", cex=2.8)
  dev.off()
  
  xx <- rbind(xx,xnew) ; yy <- rbind(yy,ynew)
  print(paste("xnew = ",xnew))
}



