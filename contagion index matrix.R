#define the binary matrix
tmat = { matrix(c( 0,0,1,0,0,0,0,0,0,0,
                   0,0,1,0,1,0,0,0,0,1,
                   1,0,0,0,0,1,0,0,1,1,
                   1,0,1,0,0,1,0,0,0,0,
                   0,0,1,0,0,0,0,0,1,0,
                   1,1,1,1,1,0,0,0,1,1,
                   0,0,1,0,0,1,0,1,1,1,
                   0,0,0,0,0,0,0,1,0,0,
                   1,0,0,0,0,1,0,0,0,1,
                   0,0,0,0,0,0,0,1,0,0),nr=10,byrow=TRUE) }

#run permutations of the matrix
library(vegan)
set.seed(4)
sim <- permatfull(tmat, fixedmar = "none", shuffle = "both", times = 9999, burnin = 20000, thin = 500, mtype = "prab")
summary(sim)

#generate list of matrices
randomat<- list(x=sim$perm)
summary(randomat)

#convert matrices to rasters
library(raster)
library(sp)
library(rgdal)
library(rlist)

#in the list of simulated matrices
a<-list()
for (i in 1:length(randomat$x)){
  a<-list.append(a,raster(randomat$x[[i]]))
}
summary(a)
length(a)

#for the observed data
observed<-raster(tmat)
image(tmat)


#estimate contagion index
library(landscapemetrics)
contagion.index <- lsm_l_contag(observed,verbose = TRUE)

#for the simulated matrices
for (i in 1:length(a)){
  contagion.index<- rbind(contagion.index,lsm_l_contag(a[[i]],verbose = TRUE))
}

#removed the observed
contagion.index.sim <- contagion.index[-c(1),] 
summary(contagion.index.sim)

mean.contagio<-mean(contagion.index.sim$value)
mean.contagio
sd.contagio<-sd(contagion.index.sim$value)
sd.contagio

x <- seq(-4,4,length=100)*sd.contagio + mean.contagio
hx <- dnorm(x,mean.contagio,sd.contagio)

#plot
plot(x, hx, type="n", xlab="Contagion index", ylab="",
     main="Normal Distribution")
i <- x <= contagion.index$value[1]
polygon(c(x[i],contagion.index$value[1]), c(hx[i],0), col="red")
area <- 1-pnorm(contagion.index$value[1], mean.contagio, sd.contagio)
result <- paste("P(Observed)=",
                signif(area, digits=3))
mtext(result,3)
lines(x, hx)
lines(density(contagion.index$value), col="blue")
abline(v=contagion.index$value[1])

