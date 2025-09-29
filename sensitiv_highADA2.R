# sensitivity analysis
# high ADA-titre (>1000 AU/mL)
# with memory T and memory B cells

library(deSolve)
library(FME)
# model without MTX
# B10 reg. cells are assumed constant

kappa0 =1e-6; # rogers
kappa1 = 1e-9; # rogers pg/ml -> microg/ml
kappa3 = 1e-9; # rogers
kappa5=2.5e-8; # assumed

dose_U=50
kabsorpU=0.28
tau_U=14
d_U = 0.03

d_a=161
k_100=1010
K_B0=5000
K_H0=24000
K_NK=1e5
d_A=2.07
p_A=0.004

KD0=5e4
KBm=100
KHm=100
k_F=9.7e-6
mD=0.2
mD1=0.75
mB=0.032
mBa=0.032

mHa=0.05
mHe=0.113
mNK=0.5

mP=0.01
mPs=0.13

alpha_1=0.5
alpha_2=0.064

aB_T=1
aBm_T=1.5
aT_APC=0.004
aTm_APC=0.004

r_Hm=0.047
r_Bm=0.5

k_101=0.166

model1<-function(pars){

dxdt <- function(t, x, pars) {
  with(as.list(c(x,pars)), {
    # variables
    U=x[1]; A=x[2]; 
    Ha=x[3]; Hm=x[4]; Hma =x[5]; 
    Ba = x[6]; Bm=x[7]; Bma=x[8]; 
    Ps=x[9]; P=x[10];
    NK0 =x[11]; NK1=x[12];
    D0 = x[13]; Dast = x[14]; D1 = x[15];
    H0=x[16];  He=x[17];
    B0 =x[18];
    TNFa=x[19]; 
    
    # sums cell total
    Bsum=B0+Ba+Bm+Bma;
    Hsum=H0+Ha+Hm+Hma+He;
    
    # APC
    # immature APC
    dD0=mD*KD0*(1-D0/KD0)-aAPC_U*D0*U;
    
    dDast = aAPC_U*D0*U-alpha_1*(1+TNFa/(TNFa+k_F))*Dast;
    
    ## mature APC
    
    dD1 = alpha_1*(1+TNFa/(TNFa+k_F))*Dast-mD1*D1/(1+TNFa/(TNFa+k_F));
    
    ## CD4 T cells
    
    dH0 = -aT_APC*D1/(D1+h_n*Hsum)*(TNFa/(TNFa+k_F)+1)*H0;
    
    dHa = aT_APC*D1/(D1+h_n*Hsum)*(TNFa/(TNFa+k_F)+1)*H0 + q_m*r_Ha*(D1/(D1+h_n*Hsum))*Ha - mHa*Ha/(1+TNFa/(TNFa+k_F));
    
    dHm = (1-q_m)*p_m*((D1/(D1+h_n*Hsum))*r_Ha*Ha + (D1/(D1+h_m*Hsum))*r_Hm*Hma)-aTm_APC*D1/(D1+h_m*Hsum)*(TNFa/(TNFa+k_F)+1)*Hm + aTm_APC*(TNFa/(TNFa+k_F)+1)*Hm*(1-Hm/KHm);
    
    dHma = aTm_APC*D1/(D1+h_m*Hsum)*(TNFa/(TNFa+k_F)+1)*Hm + q_m*r_Hm*(D1/(D1+h_m*Hsum))*Hma - mHa*Hma/(1+TNFa/(TNFa+k_F));
    
    dHe = (1-p_m)*(1-q_m)*((D1/(D1+h_n*Hsum))*r_Ha*Ha + (D1/(D1+h_m*(Hsum)))*r_Hm*Hma)- mHe*He/(1+TNFa/(TNFa+k_F)) ;
    
    
    # B cells
    
    dB0 = -aB_T*He/(He+hB_n*Bsum)*B0 ;
    
    dBa = aB_T*He/(He+hB_n*Bsum)*B0 +s_m*r_Ba*(He/(He+hB_n*Bsum))*Ba -mBa*Ba -r_Ps*Ba;
    
    dBm = r_m*(1-s_m)*(r_Ba*(He/(He+hB_n*Bsum))*Ba + r_Bm*(He/(He+hB_m*Bsum))*Bma) - aBm_T*He/(He+hB_m*Bsum)*Bm + aBm_T*Bm*(1-Bm/KBm);
    
    dBma = aBm_T*He/(He+hB_m*Bsum)*Bm +s_m*r_Bm*(He/(He+hB_m*Bsum))*Bma-mBa*Bma;
    
    # short-lived plasma cells
    dPs = r_Ps*Ba -mPs*Ps;
    
    dP= (1-r_m)*(1-s_m)*(r_Ba*(He/(He+hB_n*Bsum))*Ba + r_Bm*(He/(He+hB_m*Bsum))*Bma)-mP*P;
    
    # NK cells
    
    dNK0 = mNK*(K_NK-NK0);
    
    dNK1= alpha_2*(1+TNFa/(TNFa+k_F))*(NK0-NK1)-mNK*NK1;
    
    # TNFa
    
    dTNFa = kappa0+kappa1*(D0+D1)+kappa3*Ha+kappa5*NK1-d_a*TNFa-k_100*TNFa*U;
    
    # adalimumab [anti-TNFa]
    U0 =dose_U*(1-exp(-kabsorpU*(floor(t/tau_U)+1)*tau_U))/(1-exp(-kabsorpU*tau_U))*exp(-kabsorpU*(t%%tau_U))
    
    dU = kabsorpU*U0/V_U-k_101*A*U-d_U*U;
    
    #  anti-drug antibody
    
    dA = p_A*(Ps +P)-d_A*A;
    
    return(list(c(dU, dA, dHa, dHm,dHma, dBa, dBm,dBma,dPs,
                  dP, dNK0, dNK1, dD0, dDast, 
                  dD1, dH0, dHe, dB0, dTNFa)))
  })
}

# initial conditions
y<- c(U=0, A=0, Ha=0, Hm=0,Hma=0, Ba=0,Bm=0, Bma=0, Ps=0, P=0, NK0=K_NK, NK1=6e4, D0=KD0, Dast=0,
      D1=0, H0= K_H0, He=0, B0=K_B0, TNFa=9.629e-6)

timepts <- seq(0, 185,by=0.1)

ii<- which(timepts %in% seq(3.5, 185, by = 14))
out <- ode(y = y, parms = pars, times = timepts, func = dxdt)

as.data.frame(out[ii,])
}

parms1<-c(h_n=0.01,
          h_m=0.01,
	hB_n=0.1,
	hB_m=0.07,
	V_U=10,
	p_m=0.5,
	q_m=0.5,
	r_m=0.5,
	s_m=0.5,
	r_Ps=0.5,
	r_Ha=0.05,
	r_Ba=0.35,
	aAPC_U=0.025,
	aB_T=1,
	aBm_T=1.5,
	aT_APC=0.004,
	aTm_APC=0.004
)

# plot results
out1<-model1(pars=parms1)
par(mfrow = c(1, 2))
plot(out1$time,out1$A,main="ADA",ylab="A(t)",xlab="time",type="l",lty=1,col="red")
plot(out1$time,out1$U,main="adalimumab",ylab="U(t)",xlab="time", type="l",lty=1,col="blue")
par(mfrow = c(1, 1))

# export sensitivity indices
Sfun1<-sensFun(func=model1,parms=parms1,
               sensvar=c("A","U"),tiny=1e-5)
summary(Sfun1,var=TRUE)

plot(Sfun1, which = c("U"), xlab="time", lwd = 2)
png(filename="RplotADAlow_lowcoh_pairs.png",width=2000,height=1500,units="px")
pairs(Sfun1, which = c("A", "U"), col = c("blue", "green"))
dev.off()
write.csv(as.data.frame(Sfun1),"sens_ADAlow2.csv")


