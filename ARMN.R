#code to run the ARMN model from Section 4 of the main text
#set working directory to source file location
rm(list=ls())
library(nimble)
library(coda)
library(spdep)
library(maptools)
library(tmap)    # for static and interactive maps
library(leaflet)
library(MASS)

memory.limit(size=1E10)

#load in the data 
load("dzc_data.RData")

#clean the data
tmax<-matrix(temp_bairro$tmax,nrow=160,byrow=TRUE)
tmin<-matrix(temp_bairro$tmin,nrow=160,byrow=TRUE)



z<-matrix(casos_pop_bairro$zika,nrow=160,byrow=TRUE)
d<-matrix(casos_pop_bairro$dengue,nrow=160,byrow=TRUE)
c<-matrix(casos_pop_bairro$chik,nrow=160,byrow=TRUE)
total <- z+d+c

#everything needs to be aggregated to 2 wks
tmax_2wk <- matrix(nrow=160,ncol=52)
tmin_2wk <- matrix(nrow=160,ncol=52)
z_2wk <- matrix(nrow=160,ncol=52)
d_2wk <- matrix(nrow=160,ncol=52)
c_2wk <- matrix(nrow=160,ncol=52)
total_2wk <- matrix(nrow=160,ncol=52)
for(i in 1:160){
  for(t in 1:52){
    total_2wk[i,t] <- total[i,2*t-1]+total[i,2*t]
    c_2wk[i,t] <- c[i,2*t-1]+c[i,2*t]
    d_2wk[i,t] <- d[i,2*t-1]+d[i,2*t]
    z_2wk[i,t] <- z[i,2*t-1]+z[i,2*t]
    tmin_2wk[i,t] <- (tmin[i,2*t-1]+tmin[i,2*t])/2
    tmax_2wk[i,t] <- (tmax[i,2*t-1]+tmax[i,2*t])/2
  }
}

sum(total_2wk==c_2wk+d_2wk+z_2wk) == 52*160

#we don't need originals just rename
total <- total_2wk
rm(total_2wk)
c <- c_2wk
rm(c_2wk)
d <- d_2wk
rm(d_2wk)
z <- z_2wk
rm(z_2wk)
tmin <- tmin_2wk
rm(tmin_2wk)
tmax <- tmax_2wk
rm(tmax_2wk)

#covariates

#maximum temperature (celcius)
hist(tmax)
#reasonably normal
mtmax <- mean(tmax)

#green area
verde <- casos.s49$verde
hist(verde)
hist(log(verde))
# this is a percentage so makes sense to multiply by 100
verde <- verde*100
mverde <- mean(verde)

#favela
favela <- casos.s49$favela
hist(favela)
hist(log(favela+1))
favela <- favela*100
mfavela <- mean(favela)

#ids_p 
#this is a standardized social development index
ids_p <- casos.s49$ids_p
mean(ids_p)
sd(ids_p)


##population density
pop <- casos.s49$pop_bairro
#we will need a map
shape_rj$oo <- seq(1:160)
library(raster)
area(shape_rj)
#the area is already there, so we will use the existing
#this should give # of people per sq km
shape_rj$pop_dens <- shape_rj$Npop/(shape_rj$SHAPEarea/1000000)
map_data <- shape_rj@data
map_data$CODBAIRRO <- as.numeric(levels(map_data$CODBAIRRO))[map_data$CODBAIRRO]
map_data <- map_data[order(map_data$CODBAIRRO),]
merge_data <- map_data[,c("CODBAIRRO","pop_dens")]
colnames(merge_data)[1] <- "CodBairro"
casos.s49 <- merge(casos.s49,merge_data,by=c("CodBairro"),all.x=TRUE)
casos.s49$pop_dens[13] <- casos.s49$pop_bairro[13]/1.2
#it is 100 people per sq km
pop_dens <- casos.s49$pop_dens/100
hist(pop_dens)
mpop_dens <- mean(pop_dens)

#need a neighborhood matrix
shape_rj2$oo <- seq(1:160)
shape_rj2$CODBAIRRO <- as.numeric(levels(shape_rj2$CODBAIRRO))[shape_rj2$CODBAIRRO]
NList <- poly2nb(shape_rj2,row.names=shape_rj2$CODBAIRRO)
NMatrix <- nb2mat(NList,zero.policy = TRUE,style = "B")
NMatrix <- rbind(NMatrix,rep(0,160))
NMatrix <- cbind(NMatrix,rep(0,161))
row.names(NMatrix)[161] <- "13"
#order properly
NMatrix <- NMatrix[order(as.numeric(row.names(NMatrix))),order(as.numeric(row.names(NMatrix)))]
sum(NMatrix[1:160,1:160] != t(NMatrix[1:160,1:160]))
min(rowSums(NMatrix))
#no need for 161
NMatrix <- NMatrix[1:160,1:160]
sum(NMatrix[1:160,1:160] != t(NMatrix[1:160,1:160]))
min(rowSums(NMatrix))
#have to deal with areas with no neighbors
#add centro to paqueta
NMatrix[13,5] <- 1
NMatrix[5,13] <- 1
#Connect University and Mare
NMatrix[105,157] <- 1
NMatrix[157,105] <- 1
sum(NMatrix[1:160,1:160] != t(NMatrix[1:160,1:160]))
min(rowSums(NMatrix))
#so symmetric and now every area has at least one neighbor


sum(total==c+d+z) == 52*160

y <- array(dim=c(3,160,52))
y[1,,] <- d
y[2,,] <- z
y[3,,] <- c

#now make neighborhood prevalence covariate
nei_prev <- array(dim=c(3,160,52))

for(k in 1:3){
  for(i in 1:160){
    nei <- which(NMatrix[i,]==1)
    for(t in 1:52){
      popnei <- 0
      casenei <- 0
      for(j in nei){
        popnei <- popnei+pop[j]
        casenei <- casenei+y[k,j,t]
      }
      nei_prev[k,i,t] <- casenei/(popnei/10000)
    }
  }
}

hist(nei_prev[2,,])
hist(log(nei_prev[3,,]+1))

nei_prev[1,,] <- log(nei_prev[1,,]+1)
nei_prev[2,,] <- log(nei_prev[2,,]+1)
nei_prev[3,,] <- log(nei_prev[3,,]+1)

hist(nei_prev[3,,])
hist(nei_prev[1,,])
hist(nei_prev[2,,])
mnei_prev1 <- mean(nei_prev[1,,])
mnei_prev2 <- mean(nei_prev[2,,])
mnei_prev3 <- mean(nei_prev[3,,])

#indicator far at least one case in neighboring area
nei_1case <- array(dim=c(3,160,52))
for(k in 1:3){
  for(i in 1:160){
    nei <- which(NMatrix[i,]==1)
    for(t in 1:52){
      casenei <- 0
      for(j in nei){
        casenei <- casenei+y[k,j,t]
      }
      nei_1case[k,i,t] <- ifelse(casenei>0,1,0)
    }
  }
}

lpsi <- log(y+1)
mlpsi1 <- mean(lpsi[1,,])
mlpsi2 <- mean(lpsi[2,,])
mlpsi3 <- mean(lpsi[3,,])

dengeeConsts <- list(mT=52,N=160,total=total,psi=y,
                     tmax=tmax,mtmax=mean(tmax),
                     verde=verde,mverde=mean(verde),
                     favela=favela,mfavela=mean(favela),
                     ids_p=ids_p,pop_dens=pop_dens,
                     mpop_dens=mean(pop_dens),
                     nei_prev=nei_prev,
                     mnei_prev1=mnei_prev1,mnei_prev2=mnei_prev2,
                     mnei_prev3=mnei_prev3,
                     lpsi=lpsi,mlpsi1=mlpsi1,
                     mlpsi2=mlpsi2,mlpsi3=mlpsi3)

dengeeData <- list(y=y)

dengeeCode <- nimbleCode({
  
  #likelihood y
  for(i in 1:N){
    for(t in 2:mT){
      
      y[1:3,i,t] ~ dmulti(size=total[i,t],prob=pi[1:3,i,t])
      
      lp2p1[i,t] <- phi[1,i,t]
      lp3p1[i,t] <- phi[2,i,t]
      
      phi[1:2,i,t] ~ dmnorm(mean=mu[1:2,i,t],prec=Prec[1:2,1:2])
      
      mu[1,i,t] <- alpha02[i]+alpha12*(tmax[i,t-1]-mtmax)+
        alpha92*(lpsi[2,i,t-1]-mlpsi2)+alpha10*(lpsi[1,i,t-1]-mlpsi1)+
        alpha62*(nei_prev[2,i,t-1]-mnei_prev2)+
        alpha7*(nei_prev[1,i,t-1]-mnei_prev1)+
        alpha82*(lpsi[3,i,t-1]-mlpsi3)
      
      mu[2,i,t] <- alpha03[i]+alpha13*(tmax[i,t-1]-mtmax)+
        alpha93*(lpsi[3,i,t-1]-mlpsi3)+alpha10*(lpsi[1,i,t-1]-mlpsi1)+
        alpha63*(nei_prev[3,i,t-1]-mnei_prev3)+
        alpha7*(nei_prev[1,i,t-1]-mnei_prev1)+
        alpha83*(lpsi[2,i,t-1]-mlpsi2)
      
      pi[1,i,t] <- 1/(1+exp(lp2p1[i,t])+exp(lp3p1[i,t]))
      pi[2,i,t] <- exp(lp2p1[i,t])/(1+exp(lp2p1[i,t])+exp(lp3p1[i,t]))
      pi[3,i,t] <- exp(lp3p1[i,t])/(1+exp(lp2p1[i,t])+exp(lp3p1[i,t]))
    }
  }
  #random intercepts
  for(i in 1:N){
    alpha02[i]~dnorm(malpha02+alpha22*(verde[i]-mverde)+
                       alpha32*(favela[i]-mfavela)+
                       alpha42*ids_p[i]+
                       alpha52*(pop_dens[i]-mpop_dens),prec_02)
    alpha03[i]~dnorm(malpha03+alpha23*(verde[i]-mverde)+
                       alpha33*(favela[i]-mfavela)+
                       alpha43*ids_p[i]+
                       alpha53*(pop_dens[i]-mpop_dens),prec_03)
  }
  
  ##priors
  #Precision matrix
  Prec[1:2,1:2] ~ dwish(R=R[1:2,1:2],df=3)
  R[1,1:2] <- c(1,0)
  R[2,1:2] <- c(0,1)
  CovM[1:2,1:2] <- inverse(Prec[1:2,1:2])
  sdz <- sqrt(CovM[1,1])
  sdc <- sqrt(CovM[2,2])
  rhozc <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
  #spatial intercepts
  prec_02~dgamma(.1,.1)
  sigma_02 <- 1/sqrt(prec_02)
  prec_03~dgamma(.1,.1)
  sigma_03 <- 1/sqrt(prec_03)
  #regression priors
  malpha02 ~ dnorm(0,sd=5)
  malpha03 ~ dnorm(0,sd=5)
  alpha12 ~ dnorm(0,sd=5)
  alpha13 ~ dnorm(0,sd=5)
  alpha22 ~ dnorm(0,sd=5)
  alpha23 ~ dnorm(0,sd=5)
  alpha32 ~ dnorm(0,sd=5)
  alpha33 ~ dnorm(0,sd=5)
  alpha42 ~ dnorm(0,sd=5)
  alpha43 ~ dnorm(0,sd=5)
  alpha52 ~ dnorm(0,sd=5)
  alpha53 ~ dnorm(0,sd=5)
  alpha62 ~ dnorm(0,sd=5)
  alpha63 ~ dnorm(0,sd=5)
  alpha7 ~ dnorm(0,sd=5)
  alpha82 ~ dnorm(0,sd=5)
  alpha83 ~ dnorm(0,sd=5)
  alpha92 ~ dunif(0,1)
  alpha93 ~ dunif(0,1)
  alpha10 ~ dunif(-1,0)
  
})


inits <- list(phi=array(data=rnorm(n=2*160*52,mean=0,sd=1),dim=c(2,160,52)),
              Prec=diag(2),
              malpha02 = rnorm(n=1,mean=0,sd=1),
              malpha03 = rnorm(n=1,mean=0,sd=1),
              alpha12 = rnorm(n=1,mean=0,sd=1),
              alpha13 = rnorm(n=1,mean=0,sd=1),
              alpha22 = rnorm(n=1,mean=0,sd=1),
              alpha23 = rnorm(n=1,mean=0,sd=1),
              alpha32 = rnorm(n=1,mean=0,sd=1),
              alpha33 = rnorm(n=1,mean=0,sd=1),
              alpha42 = rnorm(n=1,mean=0,sd=1),
              alpha43 = rnorm(n=1,mean=0,sd=1),
              alpha52 = rnorm(n=1,mean=0,sd=1),
              alpha53 = rnorm(n=1,mean=0,sd=1),
              alpha62 = rnorm(n=1,mean=0,sd=1),
              alpha63 = rnorm(n=1,mean=0,sd=1),
              alpha7 = rnorm(n=1,mean=0,sd=1),
              alpha82 = rnorm(n=1,mean=0,sd=1),
              alpha83 = rnorm(n=1,mean=0,sd=1),
              alpha92 = runif(n=1,min = 0,max=1),
              alpha93 = runif(n=1,min = 0,max=1),
              alpha10 = runif(n=1,min = -1,max=0),
              alpha02 = rnorm(n=160,mean=0,sd=1),
              alpha03 = rnorm(n=160,mean=0,sd=1),
              prec_02 = runif(n=1,min=0,max=10),
              prec_03 = runif(n=1,min=0,max=10))

dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)

Cdengee <- compileNimble(dengeemodel)


dengeeConf <- configureMCMC(dengeemodel, print = TRUE)

#these parameters exhibit slow mixing and high posterior correlations
dengeeConf$removeSampler(c("alpha42","alpha52","alpha43","alpha53"))
dengeeConf$addSampler(target = c("alpha42"),type = "slice")
dengeeConf$addSampler(target = c("alpha52"),type = "slice")
dengeeConf$addSampler(target = c("alpha43"),type = "slice")
dengeeConf$addSampler(target = c("alpha53"),type = "slice")

dengeeConf$addMonitors(c("sigma_02","sigma_03","alpha02","alpha03","phi",
                         "sdz","sdc","rhozc"))

print(dengeeConf)

dengeeMCMC <- buildMCMC(dengeeConf)

CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)

initsFunction <- function()  list(phi=array(data=rnorm(n=2*160*52,mean=0,sd=1),dim=c(2,160,52)),
                                  Prec=diag(2),
                                  malpha02 = rnorm(n=1,mean=0,sd=1),
                                  malpha03 = rnorm(n=1,mean=0,sd=1),
                                  alpha12 = rnorm(n=1,mean=0,sd=1),
                                  alpha13 = rnorm(n=1,mean=0,sd=1),
                                  alpha22 = rnorm(n=1,mean=0,sd=1),
                                  alpha23 = rnorm(n=1,mean=0,sd=1),
                                  alpha32 = rnorm(n=1,mean=0,sd=1),
                                  alpha33 = rnorm(n=1,mean=0,sd=1),
                                  alpha42 = rnorm(n=1,mean=0,sd=1),
                                  alpha43 = rnorm(n=1,mean=0,sd=1),
                                  alpha52 = rnorm(n=1,mean=0,sd=1),
                                  alpha53 = rnorm(n=1,mean=0,sd=1),
                                  alpha62 = rnorm(n=1,mean=0,sd=1),
                                  alpha63 = rnorm(n=1,mean=0,sd=1),
                                  alpha7 = rnorm(n=1,mean=0,sd=1),
                                  alpha82 = rnorm(n=1,mean=0,sd=1),
                                  alpha83 = rnorm(n=1,mean=0,sd=1),
                                  alpha92 = runif(n=1,min = 0,max=1),
                                  alpha93 = runif(n=1,min = 0,max=1),
                                  alpha10 = runif(n=1,min = -1,max=0),
                                  alpha02 = rnorm(n=160,mean=0,sd=1),
                                  alpha03 = rnorm(n=160,mean=0,sd=1),
                                  prec_02 = runif(n=1,min=0,max=10),
                                  prec_03 = runif(n=1,min=0,max=10))

samples <- runMCMC(CdengeeMCMC,  niter =300000,nchains = 3,nburnin=50000
                   ,samplesAsCodaMCMC = TRUE,thin=25,inits = initsFunction)


#check convergence
param_list <- c("malpha02","malpha03","sigma_02","sigma_03",
                "alpha12","alpha13","alpha22","alpha23",
                "alpha32","alpha33",
                "alpha42","alpha43",
                "alpha52","alpha53",
                "alpha62","alpha63","alpha7",
                "alpha82","alpha83",
                "alpha92","alpha93","alpha10",
                "Prec[1, 1]","Prec[1, 2]",
                "Prec[2, 2]")
#should all be less then 1.05
gelman.diag(samples[,param_list])
alpha02_list <- paste0("alpha02[",1:160,"]")
alpha03_list <- paste0("alpha03[",1:160,"]")
param_list <- c(param_list,alpha02_list,alpha03_list)
gelman.diag(samples[,param_list])
#should all be greater than 1000
effectiveSize(samples[,param_list])

#check the space-time residuals
phi_gr <- array(dim=c(2,160,52))
for(k in 1:2){
  for(i in 1:160){
    print(i)
    for(t in 2:52){
      phi_gr[k,i,t] <- gelman.diag(samples[,c(paste0("phi[",k,", ",i,", ",t,"]"))])$psrf[1,2]
    }
  }
}
max(phi_gr,na.rm = TRUE)

#view traceplots
plot(samples[,c("malpha02","malpha03")])
plot(samples[,c("sigma_02","sigma_03")])
plot(samples[,c("Prec[1, 1]","Prec[1, 2]")])
plot(samples[,c("Prec[2, 1]","Prec[2, 2]")])
plot(samples[,c("sdz","sdc","rhozc")])
plot(samples[,c("alpha12","alpha13")])
plot(samples[,c("alpha22","alpha23")])
summary(samples[,c("alpha22","alpha23")])
plot(samples[,c("alpha32","alpha33")])
summary(samples[,c("alpha32","alpha33")])
plot(samples[,c("alpha42","alpha43")])
summary(samples[,c("alpha42","alpha43")])
plot(samples[,c("alpha52","alpha53")])
summary(samples[,c("alpha52","alpha53")])
plot(samples[,c("alpha82","alpha83")])
plot(samples[,c("alpha02[101]","alpha03[101]")])
plot(samples[,c("alpha02[2]","alpha03[2]")])
plot(samples[,c("alpha02[53]","alpha03[53]")])
plot(samples[,c("alpha02[150]","alpha03[150]")])
plot(samples[,c("alpha62","alpha63","alpha7")])
plot(samples[,c("alpha92","alpha93","alpha10")])

#look at fit and calc WAIC
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
psi <- y
mtmax <- mean(tmax)

#WAIC
lppd <- 0 
pwaic <- 0

for(i in 1:160){
  print(i)
  for(t in 2:52){
    
    if(total[i,t]>0){
      
      lp2p1 <- samps[,paste0("phi.",1,"..",i,"..",t,".")]
      lp3p1 <- samps[,paste0("phi.",2,"..",i,"..",t,".")]
      
      pi1 <- 1/(1+exp(lp2p1)+exp(lp3p1))
      pi2 <- exp(lp2p1)/(1+exp(lp2p1)+exp(lp3p1))
      pi3 <- exp(lp3p1)/(1+exp(lp2p1)+exp(lp3p1))
      
      likef <- rep(NA,30000)
      for(s in 1:30000){
        likef[s] <- dmultinom(x=y[1:3,i,t],size=total[i,t],prob=c(pi1[s],pi2[s],pi3[s]))
      }
      
      lppd <- lppd + log(mean(likef))
      pwaic <- pwaic + var(log(likef))
    }
  }
}

waic <- -2*(lppd-pwaic)
#should be around 23,974

#draw the fitted values as described in SM Section 5
psi <- y
fittedy <- array(dim = c(3,160,52,30000))

#store the Cov matrix
Covm <- array(dim=c(30000,2,2))
for(m in 1:30000){
  Precm <- matrix(nrow=2,ncol=2)
  Precm[1, 1] <- samps[m,paste0("Prec.",1,"..",1,".")]
  Precm[1, 2] <- samps[m,paste0("Prec.",1,"..",2,".")]
  Precm[2, 1] <- samps[m,paste0("Prec.",2,"..",1,".")]
  Precm[2, 2] <- samps[m,paste0("Prec.",2,"..",2,".")]
  Covm[m,,] <- inverse(Precm)
  
}

for(i in 1:160){
  print(i)
  for(t in 2:52){
    
    
    
    psi1it_star <- rep(0,30000)
    psi2it_star <- rep(0,30000)
    for(m in 1:30000){
      errr <- mvrnorm(n=1,mu=c(0,0),Sigma=Covm[m,,])
      psi1it_star[m] <- errr[1]
      psi2it_star[m] <- errr[2]
    }
    
    lp2p1 <- as.numeric(unlist(samps[paste0("alpha02.",i,".")]+
                                 samps["alpha12"]*(tmax[i,t-1]-mtmax)+
                                 samps["alpha92"]*(lpsi[2,i,t-1]-mlpsi2)+samps["alpha10"]*(lpsi[1,i,t-1]-mlpsi1)+
                                 samps["alpha62"]*(nei_prev[2,i,t-1]-mnei_prev2)+
                                 samps["alpha7"]*(nei_prev[1,i,t-1]-mnei_prev1)+
                                 samps["alpha82"]*(lpsi[3,i,t-1]-mlpsi3)+
                                 psi1it_star))
    lp3p1 <- as.numeric(unlist(samps[paste0("alpha03.",i,".")]+
                                 samps["alpha13"]*(tmax[i,t-1]-mtmax)+
                                 samps["alpha93"]*(lpsi[3,i,t-1]-mlpsi3)+samps["alpha10"]*(lpsi[1,i,t-1]-mlpsi1)+
                                 samps["alpha63"]*(nei_prev[3,i,t-1]-mnei_prev3)+
                                 samps["alpha7"]*(nei_prev[1,i,t-1]-mnei_prev1)+
                                 samps["alpha83"]*(lpsi[2,i,t-1]-mlpsi2)+
                                 psi2it_star))
    
    
    pi1 <- 1/(1+1*exp(lp2p1)+1*exp(lp3p1))
    pi2 <- (1*exp(lp2p1))/(1+1*exp(lp2p1)+1*exp(lp3p1))
    pi3 <- (1*exp(lp3p1))/(1+1*exp(lp2p1)+1*exp(lp3p1))
    
    for(s in 1:30000){
      fittedy[1:3,i,t,s] <- rmultinom(n=1,size=total[i,t],prob=c(pi1[s],pi2[s],pi3[s])) 
    }
    
  }
}

#plot the fitted values 
d <- y[1,,]
i <- 3

med_d <- apply(fittedy[1,i,,],MARGIN = 1,mean)
d_lower <- apply(fittedy[1,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
d_upper <- apply(fittedy[1,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
med_z <- apply(fittedy[2,i,,],MARGIN = 1,mean)
z_lower <- apply(fittedy[2,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
z_upper <- apply(fittedy[2,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))
med_c <- apply(fittedy[3,i,,],MARGIN = 1,mean)
c_lower <- apply(fittedy[3,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.025),na.rm=TRUE))
c_upper <- apply(fittedy[3,i,,],MARGIN = 1,function(x) quantile(x,probs=c(.975),na.rm=TRUE))

#compare fit to empirical estimates
pi1_empirical <- sum(d[i,])/sum(total[i,])
pi2_empirical <- sum(z[i,])/sum(total[i,])
pi3_empirical <- sum(c[i,])/sum(total[i,])


par(mfrow=c(1,1))
plot(d[i,],ylim = c(0,max(d_upper,d[i,],pi1_empirical*total[i,],na.rm = TRUE)),
     main=paste0(i))
lines(med_d,col="red")
lines(d_lower,lty = 2,col="red")
lines(d_upper,lty = 2,col="red")
#lines(pi1_empirical*total[i,])
par(mfrow=c(1,1))
plot(c[i,],ylim = c(0,max(c_upper,c[i,],na.rm = TRUE)),
     main=paste0(i))
lines(med_c,col="red")
lines(c_lower,lty = 2,col="red")
lines(c_upper,lty = 2,col="red")
par(mfrow=c(1,1))
plot(z[i,],ylim = c(0,max(z_upper,z[i,],na.rm = TRUE)),
     main=paste0(i))
lines(med_z,col="red")
lines(z_lower,lty = 2,col="red")
lines(z_upper,lty = 2,col="red")








