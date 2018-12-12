source("~/src/git/ProjectiveDecomposition/Projective.R")
###
# Author:  Max Robinson, Institute for Systems Biology
# Version:  1
# Date:    11 Dec 2018
###

###
# Radial pattern of small circles:
#
# k radial lines
# m distances from the origin
# n points per circle
###
k     <- 7  # number of radial gridlines
m     <- 7  # number of circular gridlines
n     <- 60 # points per small circle
#rel.r <- 1/16  # relative radius of small circle
offset.r <- 1/16  # relative radius of small circle

right.angle <- pi/2
theta <- (1:k)*right.angle/(k+1) # direction of radial line (k)
primary.r <- c(1:m)              # distance from origin to center of small circle (m)
phi       <- sort(2*pi*runif(n)) # offset angles common to all small circles (n)

gridded.wheels <- function(i,aspect=1) {
  # decode i into 1 + i.n + n*(i.m + m*i.k)
  j   <- i-1
  i.n <- j %% n
  j   <- round((j - i.n) / n) # just in case floating point division is used
  i.m <- j %% m
  i.k <- round((j - i.m) / m)
  alpha <- theta[1+i.k]
  r     <- primary.r[1+i.m]
  beta  <- phi[1+i.n]
#  r * c(cos(alpha) + rel.r*cos(beta), aspect*sin(alpha) + rel.r*sin(beta))
  c(r + offset.r*cos(beta), aspect*alpha*(k+1)/right.angle + offset.r*sin(beta))
}

wheels.within.wheels <- function(i,aspect=1) {
  # decode i into 1 + i.n + n*(i.m + m*i.k)
  j   <- i-1
  i.n <- j %% n
  j   <- round((j - i.n) / n) # just in case floating point division is used
  i.m <- j %% m
  i.k <- round((j - i.m) / m)
  alpha <- theta[1+i.k]
  r     <- primary.r[1+i.m]
  beta  <- phi[1+i.n]
  #  r * c(cos(alpha) + rel.r*cos(beta), sin(alpha) + rel.r*sin(beta))
  c(r*cos(alpha) + offset.r*cos(beta), aspect*r*sin(alpha) + offset.r*sin(beta))
}

#g <- (1+sqrt(5))/2
N <- k * m * n
Cs <- hsv((2/3)*c(1:N)/(N+1),1,5/6)
M1 <- matrix(sapply(c(1:N),gridded.wheels),nrow=N,ncol=2,byrow=TRUE)
M2 <- matrix(sapply(c(1:N),wheels.within.wheels,3),nrow=N,ncol=2,byrow=TRUE)
s.M1    <- scale(M1,center=TRUE,scale=TRUE)
s.M2    <- scale(M2,center=TRUE,scale=TRUE)
s.logM1 <- scale(log(M1),center=TRUE,scale=TRUE)
s.logM2 <- scale(log(M2),center=TRUE,scale=TRUE)
pd.M1   <- projective.decomposition(M1,tol = 1e-12,verbose=TRUE)
pd.M2   <- projective.decomposition(M2,tol = 1e-12,verbose=TRUE)
W1 <- M1
W1 <- M1 / (pd.M1$rms * pd.M1$row.factor[row(M1)] * pd.M1$col.factor[col(M1)])
W2 <- M2
W2 <- M2 / (pd.M2$rms * pd.M2$row.factor[row(M2)] * pd.M2$col.factor[col(M2)])

#rms <- function(x) { sqrt(mean(x*x)) }
xyatan2 <- function(x) { atan2(x[2],x[1])}

M1r <- M1
M1r[,2] <- unlist(apply(M1,1,rms))
M1r[,1] <- unlist(apply(M1,1,xyatan2))

M2r <- M2
M2r[,2] <- unlist(apply(M2,1,rms))
M2r[,1] <- unlist(apply(M2,1,xyatan2))


###
### Gridded data
###

# Figure 1
par(mfcol=c(2,2))
xyr <- range(c(0,M1))
plot(M1,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Data',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(a=0,b=1,lty=2)

xyr <- range(c(0,W1))
plot(W1,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Projective Decomposition',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(a=0,b=1,lty=2)

xyr <- range(c(0,s.M1))
plot(s.M1,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='z-scores',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(a=0,b=1,lty=2)

xyr <- range(c(0,s.logM1))
plot(s.logM1,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='z-scores of log',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(a=0,b=1,lty=2)

# Figure 2
par(mfcol=c(2,2))
yr <- range(c(0,unlist(apply(M1,1,rms))))
plot(unlist(apply(M1,1,xyatan2)),unlist(apply(M1,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi,  ylim=yr, main='Data',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

yr <- range(c(0,pd.M1$row.factor))
plot(unlist(apply(W1,1,xyatan2)), pd.M1$row.factor, main='Projective Decomposition',
     xlab='Angle (Rad)', ylab='Row scaling factor',
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi,  ylim=yr)
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

yr <- range(c(0,unlist(apply(s.M1,1,rms))))
plot(unlist(apply(s.M1,1,xyatan2)),unlist(apply(s.M1,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='z-scores',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

yr <- range(c(0,unlist(apply(s.logM1,1,rms))))
plot(unlist(apply(s.logM1,1,xyatan2)),unlist(apply(s.logM1,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='z-scores of log',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

###
### Radially-arrayed data
###

# Figure 2a
par(mfcol=c(2,2))
xyr <- range(c(0,M2))
plot(M2,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Data',
     xlab='x', ylab='y')
abline(v=5*c(round(xyr[1]/5):round(xyr[2]/5)),lty=3)
abline(h=5*c(round(xyr[1]/5):round(xyr[2]/5)),lty=3)

xyr <- range(c(0,W2))
plot(W2,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Projective Decomposition',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)

xyr <- range(c(0,s.M2))
plot(s.M2, pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='z-scores of log',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)

xyr <- range(c(0,s.logM2))
plot(s.logM2,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='z-scores of log',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)

# Figure 2b
par(mfcol=c(2,2))
yr <- range(c(0,unlist(apply(M2,1,rms))))
plot(unlist(apply(M2,1,xyatan2)),unlist(apply(M2,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='Data',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=5*c(round(min(yr)/5):round(max(yr)/5)),lty=3)

yr <- range(c(0,pd.M2$row.factor))
plot(unlist(apply(W2,1,xyatan2)), pd.M2$row.factor, main='Projective Decomposition',
     xlab='Angle (Rad)', ylab='Row scaling factor',
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr)
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

yr <- range(c(0,unlist(apply(s.M2,1,rms))))
plot(unlist(apply(s.M2,1,xyatan2)),unlist(apply(s.M2,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='z-scores',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

yr <- range(c(0,unlist(apply(s.logM2,1,rms))))
plot(unlist(apply(s.logM2,1,xyatan2)),unlist(apply(s.logM2,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='z-scores of log',
     xlab='Angle (Rad)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)


pd.sM1 <- projective.decomposition(s.M1,tol=1e-12,verbose=TRUE)
W3 <- s.M1
W3 <- s.M1 / (pd.sM1$rms * pd.sM1$row.factor[row(W3)] * pd.sM1$col.factor[col(W3)])

# Fig 3
par(mfrow=c(2,2))
xyr <- range(c(0,s.M1))
plot(s.M1, pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Data',
     xlab='x', ylab='y')
abline(v=c(round(min(xyr)):round(max(xyr))),lty=3)
abline(h=c(round(min(xyr)):round(max(xyr))),lty=3)

yr <- range(c(0,unlist(apply(s.M1,1,rms))))
plot(unlist(apply(s.M1,1,xyatan2)),unlist(apply(s.M1,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='Data (polar)',
     xlab='Angle', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)

xyr <- range(c(0,W3))
plot(W3,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Projective Decomposition',
     xlab='x', ylab='y')
abline(v=c(round(min(xyr[1])):round(max(xyr[2]))),lty=3)
abline(h=c(round(min(xyr[1])):round(max(xyr[2]))),lty=3)

yr <- range(c(0,pd.sM1$row.factor))
plot(unlist(apply(W3,1,xyatan2)), pd.sM1$row.factor, main='Projective Decomposition (polar)',
     xlab='Angle (Rad)', ylab='Row scaling factor',
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr)
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)


pd.sM2 <- projective.decomposition(s.M2,tol=1e-12,verbose=TRUE)
W4 <- s.M2
W4 <- s.M2 / (pd.sM2$rms * pd.sM2$row.factor[row(W4)] * pd.sM2$col.factor[col(W4)])

par(mfcol=c(2,2))
xyr <- range(c(0,s.M2))
plot(s.M2, pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Z-scores',
     xlab='x', ylab='y')
abline(v=c(round(min(xyr)):round(max(xyr))),lty=3)
abline(h=c(round(min(xyr)):round(max(xyr))),lty=3)
abline(a=0,b=tan(0),lty=2)
abline(a=0,b=tan(pi/4),lty=2)
abline(a=0,b=tan(pi/2),lty=2)
abline(a=0,b=tan(3*pi/4),lty=2)

yr <- range(c(0,unlist(apply(s.M2,1,rms))))
plot(unlist(apply(s.M2,1,xyatan2)),unlist(apply(s.M2,1,rms)),
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr, main='Polar view',
     xlab='atan2(y,x)', ylab='RMS')
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=5*c(round(min(yr)/5):round(max(yr)/5)),lty=3)
abline(v=c(-1,1)*pi/4,lty=2)
abline(v=c(-1,1)*pi/2,lty=2)
abline(v=c(-1,1)*3*pi/4,lty=2)
abline(v=c(-1:1)*pi,lty=2)

xyr <- range(c(0,W4))
plot(W4,pch='.',cex=3,col=Cs, xlim=xyr, ylim=xyr, main='Projective Decomposition',
     xlab='x', ylab='y')
abline(v=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(h=c(round(xyr[1]):round(xyr[2])),lty=3)
abline(a=0,b=tan(0),lty=2)
abline(a=0,b=tan(pi/4),lty=2)
abline(a=0,b=tan(pi/2),lty=2)
abline(a=0,b=tan(3*pi/4),lty=2)

yr <- range(c(0,pd.sM2$row.factor))
plot(unlist(apply(W4,1,xyatan2)), pd.sM2$row.factor, main=' ',
     xlab='atan2(y,x)', ylab='Row scaling factor',
     pch='.',cex=3,col=Cs, xlim=c(-1,1)*pi, ylim=yr)
abline(v=c(-4:4)*pi/4,lty=3) ; abline(h=c(round(min(yr)):round(max(yr))),lty=3)
abline(v=c(-1,1)*pi/4,lty=2)
abline(v=c(-1,1)*pi/2,lty=2)
abline(v=c(-1,1)*3*pi/4,lty=2)
abline(v=c(-1:1)*pi,lty=2)

