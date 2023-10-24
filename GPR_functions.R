setwd("~/analysis/stats/gpr/")


###############
## Functions ## ---------
###############

library(MASS)

f.ker <- function(x1, x2, th, knum, f.diag=F, grad=F){
  if(knum==1){ # RBF kernel
    kk <- th[1] * exp(-th[2]*sum((x1-x2)^2))
  }
  if(knum==2){ # Linear kernel
    kk <- th[1] * sum(x1*x2)
  }
  if(knum==3){ # exponantial kernel
    kk <- th[1] * exp(-th[2] * sqrt(sum((x1-x2)^2)) )
  }
  if(knum==4){ # Matern3 kernel
    rr <- sqrt(sum((x1-x2)^2))
    kk <- (1 + sqrt(3)*rr/th[1])*exp(-sqrt(3)*rr/th[1])
  }
  if(knum==5){ # Matern5 kernel
    rr <- sqrt(sum((x1-x2)^2))
    kk <- (1 + sqrt(5)*rr/th[1] + sqrt(5)*rr^2/(3*th[1]^2))*exp(-sqrt(5)*rr/th[1])
  }
  
  dlt <- 0
  if(f.diag){
    if(knum==1 || knum==3){
      kk <- kk + th[3]
    }else{
      kk <- kk + th[2]
    }
    dlt <- 1
  }
  if(grad==T){
    kk <- c()
    if(knum==1){
      kk[1] <- f.ker(x1,x2, th, knum) - th[3]*dlt 
      kk[2] <- (f.ker(x1,x2,th, knum) - th[3]*dlt) * (-th[2]*sum((x1-x2)^2))
      kk[3] <- th[3]*dlt
    }
    if(knum==2){
      kk[1] <- th[1] * sum(x1*x2)
      kk[2] <- th[2] * dlt
    }
    if(knum==3){
      kk[1] <- f.ker(x1,x2, th, knum) - th[3]*dlt 
      kk[2] <- (f.ker(x1,x2,th, knum) - th[3]*dlt) * (-th[2]*sqrt(sum((x1-x2)^2)))
      kk[3] <- th[3]*dlt
    }
    if(knum==4){
      kk[1] <- exp(-sqrt(3)*rr/th[1])*3*rr/(th[1]^2)
      kk[2] <- th[2]*dlt
    }
    if(knum==5){
      kk[1] <- exp(-sqrt(5)*rr/th[1]) * (5*rr^2/3/(th[1]^2)+5*sqrt(5)/3*rr^3/th[1]^3)
      kk[2] <- th[2]*dlt
    }
  }
  return(kk)
}

f.gpr <- function(x1,y1,x2,th, knum){
  ymean <- mean(y1)
  y1 <- y1 - ymean
  
  kmat <- matrix(0,length(y1),length(y1))
  for(j in 1:length(y1)){
    for(i in j:length(y1)){
      kmat[i,j] <- kmat[j,i] <- f.ker(x1[i,], x1[j,], th, knum, i==j)
    }
  }
  Kinv <- solve(kmat,tol=10^-50)
  
  res <- matrix(0,nrow(x2),2)
  for(k in 1:nrow(x2)){
    k1 <- c()
    for(i in 1:length(y1)){
      k1[i] <- f.ker(x2[k,],x1[i,],th, knum)
    }
    kk <- f.ker(x2[k,],x2[k,],th, knum, T)
    
    res[k,] <- c(k1%*%Kinv%*%y1 + ymean, sqrt(max(0, kk - matrix(k1,1,length(k1))%*%Kinv%*%k1)))
  }
  
  return(res)
}

f.lh <- function(x1,y1,th, knum){
  ymean <- mean(y1)
  y1 <- y1 - ymean
  
  kmat <- matrix(0,length(y1),length(y1))
  for(j in 1:length(y1)){
    for(i in j:length(y1)){
      kmat[i,j] <- kmat[j,i] <- f.ker(x1[i,], x1[j,], th, knum, i==j)
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

hpar.mcmc <- function(x1, y1, bounds, knum){
  ## (M-H (Metropolis-Hastings) method)
  npar <- nrow(bounds)
  n.itr <- 10000
  logb <- log(bounds)
  
  logp <- (logb[,1] + logb[,2])/2
  para <- exp(logp)
  loglh <- f.lh(x1,y1,para, knum)
  
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
    loglh.n <- f.lh(x1,y1,para.n, knum)
    
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
    #  print("- - - - - - - - - - - - - - -")
      print( paste("Iteration",i, " | ", Sys.time()) )
    #  print("Parameters: ")
    #  print(exp(logp))
    #  print(paste("Likelihood: ", max(lh_list)))
    }
    
  }
  
  nmx <- (1:length(lh_list))[lh_list==max(lh_list)]
  rr_mcmc <- par_list[max(nmx),]
  return(rr_mcmc)
}

hpar.grad <- function(x1, y1, ini_th=c(), knum){
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
        kmat[i,j]  <- f.ker(x1[i,],x1[j,], th, knum, i==j)
        dkmat[ii,] <- f.ker(x1[i,],x1[j,], th, knum, i==j, grad = T)
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
    
    if(it%%1000==1){
      #print("- - - - - - - - - - - - - - -")
      print( paste("Iteration",it) )
      #print( paste("||df|| =", dd) )
      #print("Parameters: ")
      #print(aa)
      #print(paste("Likelihood: ", f.lh(x1,y1,aa, knum)))
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

vec.npar <- c(3,2,3,2,2)
##-----------------------

