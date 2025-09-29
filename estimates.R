# estimates the parameters from the low-dimensional 
# auxiliary models in the supplementary material
# uses R
# requires libraries: coda, readr, FME

library(readr)
library(deSolve)
library(rootSolve)
library(coda)
library(FME)
library(ggplot2)

# section S3.1

# create data frame

Data<-data.frame(matrix(nrow=5,ncol=4))
names(Data) <- c("time","U0", "U", "A")

t<-c(0,28,56,12*7,26*7)
Data$time[1]<-t[1]
Data$time[2:5]<-t[2:5]-0.1
# adalimumab
# mean - no ADA
Data$U0<- c(0,  5.09734375, 6.470375,7.65090625, 9.250875)
# mean - with ADA
Data$U<-c(0, 3.839333, 2.972167, 1.829833, 0.830500)
# mean - antibody (converted from AU/mL)
Data$A<-c(0,0.3008, 0.4360, 0.6280, 0.9560)

# time points for integration of ODE
tpts = seq(0, 200,by=0.01)

# parameter values
 sigma0=1e-3

# initial guess
parms0 <- c(k_101 = 0.1, sigma1=0.1,Ka=1.0,b_U = 0.23,	mu_U=0.026)

model0 <- function(t, x, parms) {
  with(as.list(parms), {
    U0=x[1]; U=x[2];   A=x[3]
    dU0dt= b_U -mu_U*U0
    dUdt = b_U -k_101*A*U -mu_U*U
    dAdt = sigma0*U +sigma1*A*(1-A/Ka)
    return(list(c(dU0dt,dUdt,dAdt)))
  })
}


ModelCost0 <- function(P) {
  
  # solve ODE for a given set of parameters
  out <- ode(y = c(0,0,0), func = model0, parms = P, times = tpts)
  outdf=data.frame(out[,])
  names(outdf)<-c("time","U0","U","A")
  outdf=outdf[outdf$time %in% Data$time,]
  return(modCost(outdf, Data)) #
}

Fit0<-modFit(f=ModelCost0,p=parms0,
             lower = c(0.06, 0,0.21, 0.16, 0.015),
             upper = c(.5, 1, 2.3,0.32, 0.035)) 

summary(Fit0)
residuals(Fit0)


out <- ode(y = c(0,0,0), func = model0, parms = Fit0$par, times = tpts)


# plot
filename<-"sectionS31.pdf"
pdf(file = filename )
par(mfrow=c(2,1))
plot(t, Data$U,ylim = c(0, max(Data$U)+1), pch = 16, col = "red",
     main = "adalimumab", xlab = "time", ylab = expression( mu*"g/ml"))
lines(out[,1],out[,3], col = "red", lty = 1)
plot(t, Data$A/0.012, ylim = c(0, (max(Data$A/0.012)+15)), pch = 16, col = "red",
     main = "antibody", xlab = "time", ylab = "AU/ml")
lines(out[,1],out[,4]/0.012, col = "red", lty = 1)
dev.off()


# section S3.2
# estimate parameters gamma4, c4 from data for T cells and MTX
# taken from Genestier et al. (1998)

# independent - time (hours)
t<-c(24,48,72,96); 
# covert to days
t<-t/24;

# observed values (T cells), no MTX
Y<-c(1.35,1.66,1.76,1.83); 

# data with MTX concentration (molar)
MTX<-c(0.01,0.1,1,10); 
# convert using molar mass
MTX<-MTX*0.4544 #microg/mL

# observed values with the respective MTX concentration
Y10<-c(0.719,0.415,0.263,0.241)
Y1<-c(0.826,0.638,0.455,0.357)
Y01<-c(0.888,0.723,0.567,0.415) 
Y001<-c(0.933,0.799,0.679,0.701)

# time range for integration
tpts <- seq(0, 96,by=0.1)

# create data frame
Data <- data.frame(matrix(nrow=4,ncol=4))

names(Data) <- c("time", "Y","Y1","Y2")
	Data$time<-t
	Data$Y<-Y
	Data$Y1<-Y01
	Data$Y2<-Y1

alpha=2
# initial guess
parms1 <- c(a=0.5,b=2.2,gamma4=0.4,c4=0.1)


modelTcell<-function(t, x, parms) {
  with(as.list(parms), {
    Y=x[1]
    Y1=x[2]
    Y2=x[3]
    dYdt =a*Y-b*Y^2
    dY1dt=a*Y1-b*Y1^2-gamma4*MTX[2]^alpha/(c4^alpha+MTX[2]^alpha)*Y1
    dY2dt=a*Y2-b*Y2^2-gamma4*MTX[3]^alpha/(c4^alpha+MTX[3]^alpha)*Y2
    return(list(c(dYdt,dY1dt,dY2dt)))
  })
}

ModelCost1 <- function(P) {
  
  # solve ODE for a given set of parameters
  out <- ode(y = c(1.0,1.0,1.0), func = modelTcell, parms = P, times = tpts)
  
  outdf=data.frame(out[,1:4])
  names(outdf)<-c("time","Y","Y1","Y2")
  outdf=outdf[outdf$time %in% Data$time,]
  return(modCost(outdf, Data))
}

fit1 <- modFit(f = ModelCost1, p = parms1, lower = c(0,0,0,0), upper = c(1,5,15,1))


# export result
summary(fit1)
plot(residuals(fit1))

out <- ode(y = c(1.0,1.0,1.0), func = modelTcell, parms = fit1$par, times = tpts)

plot(t, Data[,"Y"], ylim = c(0,2), pch = 16, col = "red",
     main = "T cells (MTX=0)", xlab = "time", ylab = "Y(t)")
lines(out[,1],out[,2], col = "red", lty = 2)

plot(t, Data[,"Y1"], ylim = c(0,1), pch = 16, col = "red",
     main = "T cells (MTX=0.1)", xlab = "time", ylab = "Y(t)")
lines(out[,1],out[,3], col = "red", lty = 2)

plot(t, Data[,"Y2"], ylim = c(0,1), pch = 16, col = "red",
     main = "T cells (MTX=1.0)", xlab = "time", ylab = "Y(t)")
lines(out[,1],out[,4], col = "red", lty = 2)


# section S3.3
# data on macrophage
# taken from Seitz et al. (1998)

# independent (time)
t<-c(0, 1, 3, 7)

# observed data
data1 <- c(0.2, 0.43, 1.53, 2.3)
data2 <- c(0.2, 0.24, 0.28, 0.31) 

# used dose of MTX
X<-9.08 # microg/ml
m_D<-0.2


# create data frame
Data <- data.frame(matrix(nrow=4,ncol=3))
names(Data) <- c("time", "Y0","Y")
    Data$time<-t
    Data$Y0<-data1
    Data$Y<-data2
# normalise data to initial value at t=0

    Data$Y0<-data1/data1[1]
    Data$Y<-data2/data2[1]
# initial guess
     parms2<-c(gamma1=1,K=2.5)
# initial value for the ODE
z0=c(0.2,0.2)/data1[1]
                   
modelMPH <- function(t, z, parms) {
  with(as.list(parms), {
    y0=z[1]
    y=z[2]
    dY0 = m_D*(K-y0)
    dY = m_D*(K-y)-gamma1*X*y
    
   return(list(c(dY0,dY)))
  })
}

ModelCost2 <- function(P) {
  
  # solve ODE for a given set of parameters
  out <- ode(y = z0, func = modelMPH, parms = P, times = t)
  
  #  # Filter data that contains time points where data is available
  outdf=data.frame(out[,1:3])
  names(outdf)<-c("time","Y0","Y")
  outdf=outdf[outdf$time %in% Data$time,]
  return(modCost(outdf, Data)) # object of class modCost
}

Fit2 <- modFit(f = ModelCost2, p = parms2, lower = c(0,1) , upper = c(15,30))


# export result
summary(Fit2)
plot(residuals(Fit2))

tpts=seq(0,50,by=0.01)

z0=c(14.8296,14.8296  )

out1 <- ode(y = z0, func = modelMPH, parms = Fit2$par, times = tpts)

plot(Data[,"time"], Data[,"Y0"], ylim = c(0, 15),xlim=c(0,max(T)), pch = 16, col = "red",
     main = "no MTX", xlab = "time", ylab = "cells")
lines(out1[,1],out1[,2], col = "red", lty = 2)

plot(Data[,"time"], Data[,"Y"], ylim = c(0, 4),xlim=c(0,max(T)), pch = 16, col = "blue",
     main = "MTX", xlab = "time", ylab = "cells")
lines(out1[,1],out1[,3], col = "blue", lty = 2)


# Section S3.4
# effect on B10 cells

# create data frame

Data <- data.frame(matrix(nrow=4,ncol=2))
names(Data) <- c("time", "B10")


# parameters
mu10<-0.032 
tau_X<-7
# parameter values
dose_MTX=10; V_X=10.03*0.8
a_X= 8.64
d_X=(2.08+3.33)/2

# all patients
dataB10<-c(4744,6387,6059,7071)
# normalise data
dataB10<-dataB10/dataB10[1]

# independents
t<-c(0,28+2*tau_X,12*7+2*tau_X,26*7+2*tau_X)

tpts <- seq(0, 200,by=0.01)


Data$time<-t
# use normalised values relative to baseline
Data$B10<-c(1,1.59,1.57,1.54)

# use averaged MTX concentration over 14 days
MTX<-0.066
model30 <-function(t, x, parms) {
  with(as.list(parms), {
    
    # variables
    B10=x[1]; 
    
    dB10=mu10*(K_10+K_10*a_10*MTX-B10)*B10
    return(list(c(dB10)))
  })
}

model31<-function(t, x, parms) {
  with(as.list(parms), {
    if (t < 2*tau_X) { MTX0 <- dose_MTX*exp(-a_X*(t))} 
      else { MTX0 <-dose_MTX*exp(-a_X*(t%%tau_X)) }
    # variables
    B10=x[1]; MTX =x[2]
    
    dMTX = a_X*MTX0/V_X-d_X*MTX
    dB10=mu10*(K_10+K_10*a_10*MTX-B10)*B10
    return(list(c(dB10,dMTX)))
  })
}

# initial guess
parms30 <- c(a_10 = 10)

ModelCost30 <- function(P) {
  
  # solve ODE for a given set of parameters
  out <- ode(y = c(1), func = model30, parms = P, times = t)
  
  #  # Filter data that contains time points where data is available
  outdf=data.frame(out[,1:2])
  names(outdf)<-c("time","B10")
  outdf=outdf[outdf$time %in% Data$time,]
  return(modCost(outdf, Data)) # object of class modCost
}

Fit3 <- modFit(f = ModelCost30, p = parms30, 
               lower = c(.06), upper = c(50))

# export result
summary(Fit3)
plot(residuals(Fit3))

out <- ode(y = c(1), func = model30, parms = Fit3$par, times = tpts)
out1 <- ode(y = c(1,0), func = model31, parms = Fit3$par, times = tpts)

# plot

plot(Data[,"time"], Data[,"B10"], ylim = c(0.95,1.9), pch = 16, col = "red",
     main = "B10 cells", xlab = "time", ylab = "B_10(t)")
lines(out[,1],out[,2], col = "blue", lty = 1)
# full model
lines(out1[,1],out1[,2], col = "red", lty = 2)



# Section S3.5
# create data frame

Data <- data.frame(matrix(nrow=4,ncol=3))

names(Data) <- c("time", "NK_tot","NK1")


# total NK cells: MTX +, ADA negative

dataNK_tot<-c(119303, 91196, 80890, 80554 )

# activated NK1 cells:  MTX +, ADA negative

dataNK1<-c(84075, 53372, 53372, 46692 )

# independent variable (time)

t<-c(0,28,12*7,26*7); 
# shift times to avoid non-smooth points
t<-c(0,28+2*tau_X,12*7+2*tau_X,26*7+2*tau_X)
t[2:4]<-t[2:4]-0.1

tpts = seq(0, 200,by=0.01)
# parameter values
dose_MTX=10; V_X=10.03*0.8
a_X= 8.64
d_X=(2.08+3.33)/2

mNK = 0.05; # 1/2-life 14 days
tau_X=7


Data$time<-t
Data$NK_tot<-dataNK_tot
Data$NK1<-dataNK1

# define initial value and carrying capacity
NK1_0<-Data$NK1[1]
NKtot_0<-Data$NK_tot[1]
K_NK<-NKtot_0
a_2<-mNK*NK1_0/(NKtot_0-NK1_0)

# use averaged MTX concentration over 14 days
MTX<-0.066
model40<-function(t, x, parms) {
  with(as.list(parms), {

    # variables
    NKtot=x[1]; NK1 =x[2] 
    
    dMTX = a_X*MTX0/V_X-d_X*MTX
    # NK cells 
    dNKtot = mNK*(K_NK-NKtot) - gamma2*MTX*NKtot
    dNK1= a_2*(NKtot-NK1) - mNK*NK1 - gamma3*MTX*NK1
    
    return(list(c(dNKtot,dNK1)))
  })
}

model41<-function(t, x, parms) {
  with(as.list(parms), {
    if (t < 2*tau_X) { MTX0 <- dose_MTX*exp(-a_X*(t))} 
        else { MTX0 <-dose_MTX*exp(-a_X*(t%%tau_X)) }

    # variables
    NKtot=x[1]; NK1 =x[2]; MTX=x[3] 
    
    dMTX = a_X*MTX0/V_X-d_X*MTX
    # NK cells
    dNKtot = mNK*(K_NK-NKtot) - gamma2*MTX*NKtot
    dNK1= a_2*(NKtot-NK1) - mNK*NK1 - gamma3*MTX*NK1
    
    return(list(c(dNKtot,dNK1,dMTX)))
  })
}

## initial guess
parms40 <- c(gamma2=2.02,gamma3=0.21)

# model cost
ModelCost40 <- function(P) {
  
  # solve ODE for a given set of parameters
  out <- ode(y = c(NKtot_0,NK1_0), func = model40, parms = P, times = tpts)
  
  #  # Filter data that contains time points where data is available
  outdf=data.frame(out[,1:3])
  names(outdf)<-c("time","NK_tot","NK1")
  outdf=outdf[outdf$time %in% Data$time,]
  return(modCost(outdf, Data)) # object of class modCost
}

Fit4 <- modFit(f = ModelCost40, p = parms40, 
               lower = c(0.01,0.05), 
               upper = c(5,4))
# summary
summary(Fit4)
plot(residuals(Fit4))

out <- ode(y = c(NKtot_0,NK1_0), func = model40, parms = Fit4$par,
           times = tpts)

out1 <- ode(y = c(NKtot_0,NK1_0,0), func = model41, parms = Fit4$par,
           times = tpts)
           
# rerun for gamma3=0 (n.s.)           
out <- ode(y = c(NKtot_0,NK1_0,0), func = model4, parms = c(gamma2=0.5105,gamma3=0),
           times = tpts)

plot(Data[,"time"], Data[,"NK_tot"], ylim = c(K_NK/2,K_NK), 
     pch = 16, col = "red",
     main = "NK cells", xlab = "time", ylab = "NK total")
lines(out[,1],out[,2], col = "blue", lty = 2)
lines(out1[,1],out1[,2], col = "red", lty = 2)

plot(t, Data[,"NK1"], ylim = c(K_NK/3,K_NK), pch = 16, col = "red",
     main = "NK1 cells", xlab = "time", ylab = "activ. NK")
lines(out[,1],out[,3], col = "blue", lty = 2)
lines(out1[,1],out1[,3], col = "red", lty = 2)

