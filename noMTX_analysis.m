
%%-------------------------------------------------------------------------------
%% use MATLAB
%% Code computes the trajectories for the virtual patient cohorts 
%% without co-medication of MTX
%% and exports the results
%%----------------------------------------------------------------------------------

% generate data
function noMTX_analysis

tend=200;

%%% Set error tolerances for the integration
options = odeset('RelTol',1e-9,'AbsTol',1e-6,'MaxStep',1/30);

alpha_1=0.5;

d_a=161;
d_U=0.03;
d_A=2.07;

k_F = 9.7e-6; %TNFa concentration = 9.7 pg/ml -> micro g/ML [Deveci]

% production of TNFa
k0_61 =1e-6; % rogers
k1_61 = 1e-9; % rogers pg/ml -> microg/ml
k3_61 = 1e-9; % rogers
k5_61=2.5e-8; % assumed

kabsorpU = 0.28; %  TYPICAL MONOCLONAL: ref. Ternant et al. ...
d_U = 0.03; %day^-1 (21 days 1/2-life)
tau_U=14;

KD0=5e4; % estimated from Dombrecht
K_B=5e3; % FROM DATA
K_H = 24000; % from data

k_100=1010;
k_101=0.166;

mD=0.2;
mD1=0.75;

K_NK=1e5; %from DATA
mNK = 0.05; % 1/2-life 14 days

mBa=0.032;

mHa=0.05;
mHe=0.113;

mP=0.01;
mPs=0.13;

dose_U=50;

% initial condition

x0 =zeros(19,1);
% D0 
x0(1) = KD0; 
% Dast = 
x0(2) =0; 
% D1 = 
x0(3) = 0;
% H0=
x0(4)=K_H; 
% Ha = Hm = Ham= He=0;
x0(5:8)=0; 

% B0 
x0(9)=K_B;
% Ba = Bm = Bam = 0;
x0(10:12)=0; 
% NK
x0(13)=K_NK;
x0(14)=60000;
% TNFa
x0(15)= k_F; % mg/L 

alpha_2=2/3*mNK*(x0(14))/(K_NK-x0(14)); % MTX- group

% titles for files:

patientno =('# virtual patient no. %g \n');
formatSpec=('%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n');

% variant with rates from estim_main.xlsx
% 

% use LH sample with 15%

filename1="LHsimul_lowADAtitre15.xlsx";

filename2="trajLHsimul_lowADAtitre15.txt";

% % values for low ADA both sub-cohorts

% ADA-low subcohort <30 AU/mL
M1=readmatrix('LH15_lowADA1.csv','Range','A1:I130');

% read values from LH
aAPC_U=M1(:,1);

aB_T=M1(:,2);
aBm_T=0.001;

aT_APC=M1(:,3);
aTm_APC=0.001;

h_n=M1(:,4);
h_m=1;
b_n=M1(:,5);
b_m=0.05;

p=M1(:,6);
q=M1(:,7);
r=0.5;
s=0.5;

r_Ba=0.35;
r_Bm=0.4;
r_Ha=M1(:,8);
r_Hm=0.02;

r_Ps=0.9;

V_U=M1(:,9);

p_A=0.00269;

KBm=10;
KHm=10;

% number of trajectories
N_cohort=size(M1,1);

% set p=0, r=0
% let kBm=KHm=1 to avoid division by zero
% simulate


% initial condition for number of samples
y0=zeros(length(x0)*N_cohort,1);
for ind=1:length(x0)
    y0(1+(ind-1)*N_cohort:ind*N_cohort)=x0(ind)*ones(N_cohort,1);
end


[tpts,out] = ode15s(@(t,x)model0(t,x,N_cohort,dose_U,alpha_1,aT_APC,aTm_APC,...
   aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s, r_Ps, r_Ha,r_Hm,...
   r_Ba,r_Bm,V_U,h_m,h_n,b_n,b_m), (0:1e-2:tend),y0,options);

% find indices

ind4=find(tpts==28);
ind8=find(tpts==54);
ind12=find(tpts==84);
ind26=find(tpts==182);
% 
drug_TNFi = out(:,15*N_cohort+1:16*N_cohort);

ADA_sol = out(:,16*N_cohort+1:17*N_cohort);

% take out values at week=4,8,12,26
sampleADAL = drug_TNFi([ind4,ind8,ind12,ind26],:);
sampleANTIB=ADA_sol([ind4,ind8,ind12,ind26],:);

sampleANTIB=sampleANTIB/0.012;
sampleANTIB=sampleANTIB';
sampleADAL=sampleADAL';


D1_sol=out(:,2*N_cohort+1:3*N_cohort);

Ha_sol=out(:,4*N_cohort+1:5*N_cohort); Hm_sol = out(:,5*N_cohort+1:6*N_cohort); 
Hma_sol=out(:,6*N_cohort+1:7*N_cohort); Heff_sol = out(:,7*N_cohort+1:8*N_cohort);
% 
Ba_sol=out(:,9*N_cohort+1:10*N_cohort); Bm_sol=out(:,10*N_cohort+1:11*N_cohort); 
Bma_sol=out(:,11*N_cohort+1:12*N_cohort); 
% 
P_sol = out(:,17*N_cohort+1:18*N_cohort);
Ps_sol = out(:,18*N_cohort+1:19*N_cohort);

NK0_sol=out(:,12*N_cohort+1:13*N_cohort); 
NK1_sol=out(:,13*N_cohort+1:14*N_cohort);

cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
    max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];
% 

writematrix(sampleANTIB,filename1,'Sheet','ADAtitre')
writematrix(sampleADAL,filename1,'Sheet','TNFi')
writematrix(cellmax,filename1,'Sheet','cellmax')

% % save trajectories for the plot
% % choose every 10th patient

% % % 
% 
fileID = fopen(filename2,'w');

fprintf(fileID,'# t D1 Ha Hm Hma Heff Ba Bm Bma P Ps NK0 NK1 TNFi ADA \n');
% % write to file time series for the virtual patients
% % 
% 
for j = (1:10:N_cohort)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:10:end), D1_sol(1:10:end,j), Ha_sol(1:10:end,j), Hm_sol(1:10:end,j),...
        Hma_sol(1:10:end,j), Heff_sol(1:10:end,j), Ba_sol(1:10:end,j), ...
        Bm_sol(1:10:end,j), Bma_sol(1:10:end,j), P_sol(1:10:end,j), ...
        Ps_sol(1:10:end,j),NK0_sol(1:10:end,j),NK1_sol(1:10:end,j), ...
 drug_TNFi(1:10:end,j),ADA_sol(1:10:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);


%
% % low ADA: subcohort >30 AU/mL <100 AU/mL

M1=readmatrix('LH15_lowADA2.csv','Range','A1:F30');

% % read values from LH
aAPC_U=0.025;

aB_T=M1(:,1);
aBm_T=1.5;

aT_APC=M1(:,2);
aTm_APC=0.004;

h_n=0.005;
h_m=0.005;
b_n=M1(:,3);
b_m=2.5;

p=M1(:,4);
q=0.5;
r=0.5;
s=0.5;

r_Ba=0.35;
r_Bm=0.5;
r_Ha=M1(:,5);
r_Hm=0.047;

r_Ps=0.5;

V_U=M1(:,6);

KBm=100;
KHm=100;

% % number of patients is 30

N_cohort=30;

% initial condition for number of samples
y0=zeros(length(x0)*N_cohort,1);
for ind=1:length(x0)
y0(1+(ind-1)*N_cohort:ind*N_cohort)=x0(ind)*ones(N_cohort,1);
end

% % simulate

[tpts,out] = ode15s(@(t,x)model0(t,x,N_cohort,dose_U,alpha_1,aT_APC,aTm_APC,...
   aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s, r_Ps,...
   r_Ha,r_Hm,r_Ba,r_Bm,V_U,h_m,h_n,b_n,b_m),(0:1e-2:tend),y0,options);

ind4=find(tpts==28);
ind8=find(tpts==54);
ind12=find(tpts==84);
ind26=find(tpts==182);
drug_TNFi = out(:,15*N_cohort+1:16*N_cohort);

ADA_sol = out(:,16*N_cohort+1:17*N_cohort);

% take out values at week=4,8,12,26
sampleADAL = drug_TNFi([ind4,ind8,ind12,ind26],:);
sampleANTIB=ADA_sol([ind4,ind8,ind12,ind26],:);

sampleANTIB=sampleANTIB/0.012;
sampleANTIB=sampleANTIB';
sampleADAL=sampleADAL';

D1_sol=out(:,2*N_cohort+1:3*N_cohort);

Ha_sol=out(:,4*N_cohort+1:5*N_cohort); Hm_sol = out(:,5*N_cohort+1:6*N_cohort); 
Hma_sol=out(:,6*N_cohort+1:7*N_cohort); Heff_sol = out(:,7*N_cohort+1:8*N_cohort);
% 
Ba_sol=out(:,9*N_cohort+1:10*N_cohort); Bm_sol=out(:,10*N_cohort+1:11*N_cohort); 
Bma_sol=out(:,11*N_cohort+1:12*N_cohort); 
% 
P_sol = out(:,17*N_cohort+1:18*N_cohort);
Ps_sol = out(:,18*N_cohort+1:19*N_cohort);

NK0_sol=out(:,12*N_cohort+1:13*N_cohort); 
NK1_sol=out(:,13*N_cohort+1:14*N_cohort);

cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
    max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];

% % % write to file
% 
writematrix(sampleANTIB,filename1,'Sheet','ADAtitre','WriteMode','append')
writematrix(sampleADAL,filename1,'Sheet','TNFi','WriteMode','append');
writematrix(cellmax,filename1,'Sheet','cellmax','WriteMode','append')

% 
% % save trajectories for the plot
% 
% % choose every 5th patient
% 
% % % % 
% % 
fileID = fopen(filename2,'a');

% % write to file time series for the virtual patients
% % 
% 
for j = (1:5:N_cohort)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:10:end), D1_sol(1:10:end,j),Ha_sol(1:10:end,j),Hm_sol(1:10:end,j),...
        Hma_sol(1:10:end,j),Heff_sol(1:10:end,j), Ba_sol(1:10:end,j), ...
        Bm_sol(1:10:end,j),Bma_sol(1:10:end,j),P_sol(1:10:end,j), ...
        Ps_sol(1:10:end,j),NK0_sol(1:10:end,j),NK1_sol(1:10:end,j),...
        drug_TNFi(1:10:end,j),ADA_sol(1:10:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);

%%----------------------------------------------------------------------
%%----------------------------------------------------------------------
%%----------------------------------------------------------------------
% % % transient immunity
% 
% % set p=0, r=0
% % let kBm=KHm=1 to avoid division by zero
% 
filename1="LHsimul_transADA15.xlsx";

filename2="trajLHsimul_transADA15.txt";

M1=readmatrix('LH15_transientADA.csv','Range','A1:E30');

% read values from LH
aAPC_U=0.01;

aB_T=0.015;
aBm_T=0;

aT_APC=0.001;
aTm_APC=0;

h_n=0.001;
b_n=0.001;
b_m=1;
h_m=1;

r_Ba=M1(:,3);
r_Bm=0;
r_Ha=M1(:,4);
r_Hm=0;

r_Ps=M1(:,1);
% 
p=0;
q=M1(:,2);
r=0;
s=0.95;

V_U=M1(:,end);

% number of virtual patients is 30
N_cohort=30;

% % simulate: initial condition
K_B=2e3;
K_H=1e4;
x0(4)=K_H; 
x0(9)=K_B;

% initial condition for number of samples
y0=zeros(length(x0)*N_cohort,1);
for ind=1:length(x0)
    y0(1+(ind-1)*N_cohort:ind*N_cohort)=x0(ind)*ones(N_cohort,1);
end

[tpts,out] = ode15s(@(t,x)model0(t,x,N_cohort,dose_U,alpha_1,aT_APC,aTm_APC,...
   aB_T,aBm_T,1,1,k_101,aAPC_U,p_A,p,q,r,s, r_Ps, r_Ha,r_Hm,r_Ba,r_Bm,...
   V_U,h_m,h_n,b_n,b_m),(0:1e-2:tend),y0,options);

ind4=find(tpts==28);
ind8=find(tpts==54);
ind12=find(tpts==84);
ind26=find(tpts==182);
drug_TNFi = out(:,15*N_cohort+1:16*N_cohort);

ADA_sol = out(:,16*N_cohort+1:17*N_cohort);

% take out values at week=,8,12,26
sampleADAL = drug_TNFi([ind4,ind8,ind12,ind26],:);
sampleANTIB=ADA_sol([ind4,ind8,ind12,ind26],:);

sampleANTIB=sampleANTIB/0.012;
sampleANTIB=sampleANTIB';
sampleADAL=sampleADAL';
D1_sol=out(:,2*N_cohort+1:3*N_cohort);

% H0_sol=out(:,4); 
Ha_sol=out(:,4*N_cohort+1:5*N_cohort); Hm_sol = out(:,5*N_cohort+1:6*N_cohort); 
Hma_sol=out(:,6*N_cohort+1:7*N_cohort); Heff_sol = out(:,7*N_cohort+1:8*N_cohort);
% 
Ba_sol=out(:,9*N_cohort+1:10*N_cohort); Bm_sol=out(:,10*N_cohort+1:11*N_cohort); 
Bma_sol=out(:,11*N_cohort+1:12*N_cohort); 
% 
P_sol = out(:,17*N_cohort+1:18*N_cohort);
Ps_sol = out(:,18*N_cohort+1:19*N_cohort);


NK0_sol=out(:,12*N_cohort+1:13*N_cohort); 
NK1_sol=out(:,13*N_cohort+1:14*N_cohort);
% 
% 

cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
    max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];

% % % write to file
% 
writematrix(sampleANTIB,filename1,'Sheet','ADAtitre','WriteMode','append')
writematrix(sampleADAL,filename1,'Sheet','TNFi','WriteMode','append');
writematrix(cellmax,filename1,'Sheet','cellmax','WriteMode','append') 

% % save trajectories for the plot
% % choose every 5th patient
% 
% % % 
% 
fileID = fopen(filename2,'w');

% % write to file time series for the virtual patients
% % 
% 
for j = (1:5:N_cohort)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:10:end), D1_sol(1:10:end,j),Ha_sol(1:10:end,j),Hm_sol(1:10:end,j),...
        Hma_sol(1:10:end,j),Heff_sol(1:10:end,j), Ba_sol(1:10:end,j), ...
        Bm_sol(1:10:end,j),Bma_sol(1:10:end,j),P_sol(1:10:end,j), ...
        Ps_sol(1:10:end,j),NK0_sol(1:10:end,j),NK1_sol(1:10:end,j),...
        drug_TNFi(1:10:end,j),ADA_sol(1:10:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);
% 
% %%----------------------------------------------------------------------
% %%----------------------------------------------------------------------
% %%----------------------------------------------------------------------
% % % high ADA cohort
% 
% % cohort >100, <1000AU/mL

filename1="LHsimul_highADAtitre15.xlsx";

filename2="trajLHsimul_highADAtitre15.txt";

M1=readmatrix('LH15_highADA1.csv','Range','A1:H70');
% 

aAPC_U=0.025;

aB_T=1;
aBm_T=1.5;

aT_APC=M1(:,1);
aTm_APC=0.004;

h_n=0.01;
h_m=0.01;
b_n=0.8;
b_m=M1(:,2);

p=M1(:,3);
q=M1(:,4);
r=M1(:,5);
s=M1(:,6);

r_Ba=0.35;
r_Bm=0.5;

r_Ha=M1(:,7);
r_Hm=0.047;

r_Ps=0.5;
p_A=0.004;

V_U=M1(:,8);



% 
% % number of virtual patients is 70
N_cohort=70;
% 
% 

% 
% % % simulate: initial condition
K_H=24000;
x0(4)=K_H; 
K_B=5e3;
x0(9)=K_B;

% % initial condition for number of samples
y0=zeros(length(x0)*N_cohort,1);
for ind=1:length(x0)
    y0(1+(ind-1)*N_cohort:ind*N_cohort)=x0(ind)*ones(N_cohort,1);
end


% % % simulate
% 
% % loop
% % 
[tpts,out] = ode15s(@(t,x)model0(t,x,N_cohort,dose_U,alpha_1,aT_APC,aTm_APC,...
   aB_T,aBm_T,100,100,k_101,aAPC_U,p_A,p,q,r,s, r_Ps, r_Ha,r_Hm,r_Ba,r_Bm,...
   V_U,h_m,h_n,b_n,b_m),(0:1e-2:tend),y0,options);


ind4=find(tpts==28);
ind8=find(tpts==54);
ind12=find(tpts==84);
ind26=find(tpts==182);
drug_TNFi = out(:,15*N_cohort+1:16*N_cohort);

ADA_sol = out(:,16*N_cohort+1:17*N_cohort);

% take out values at week=4,8,12,26
sampleADAL = drug_TNFi([ind4,ind8,ind12,ind26],:);
sampleANTIB=ADA_sol([ind4,ind8,ind12,ind26],:);

sampleANTIB=sampleANTIB/0.012;
sampleANTIB=sampleANTIB';
sampleADAL=sampleADAL';

D1_sol=out(:,2*N_cohort+1:3*N_cohort);
Ha_sol=out(:,4*N_cohort+1:5*N_cohort); Hm_sol = out(:,5*N_cohort+1:6*N_cohort); 
Hma_sol=out(:,6*N_cohort+1:7*N_cohort); Heff_sol = out(:,7*N_cohort+1:8*N_cohort);
% 
Ba_sol=out(:,9*N_cohort+1:10*N_cohort); Bm_sol=out(:,10*N_cohort+1:11*N_cohort); 
Bma_sol=out(:,11*N_cohort+1:12*N_cohort); 
% 
P_sol = out(:,17*N_cohort+1:18*N_cohort);
Ps_sol = out(:,18*N_cohort+1:19*N_cohort);

NK0_sol=out(:,12*N_cohort+1:13*N_cohort); 
NK1_sol=out(:,13*N_cohort+1:14*N_cohort);


cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
    max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];

% % % write to file
% % 
writematrix(sampleANTIB,filename1,'Sheet','ADAtitre')
writematrix(sampleADAL,filename1,'Sheet','TNFi');
writematrix(cellmax,filename1,'Sheet','cellmax') ;
% 
% 
% % save trajectories for the plot
fileID = fopen(filename2,'w');

% % write to file time series for the virtual patients
% % 
% 
for j = (1:10:N_cohort)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:10:end), D1_sol(1:10:end,j),Ha_sol(1:10:end,j),Hm_sol(1:10:end,j),...
        Hma_sol(1:10:end,j),Heff_sol(1:10:end,j), Ba_sol(1:10:end,j), ...
        Bm_sol(1:10:end,j),Bma_sol(1:10:end,j),P_sol(1:10:end,j), ...
        Ps_sol(1:10:end,j),NK0_sol(1:10:end,j),NK1_sol(1:10:end,j),...
        drug_TNFi(1:10:end,j),ADA_sol(1:10:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');
end

fclose(fileID);


% cohort > 1000 AU/mL

M1=readmatrix('LH15_highADA2.csv','Range','A1:H30');
% 
aAPC_U=0.025;

aB_T=1;
aBm_T=1.5;

aT_APC=M1(:,1);
aTm_APC=0.004;

h_n=0.01;
h_m=0.01;

b_n=0.1;
b_m=M1(:,2);

KBm=100;
KHm=100;

p=M1(:,3);
q=M1(:,4);
r=M1(:,5);
s=M1(:,6);

r_Ba=0.35;
r_Bm=0.5;

r_Ha=M1(:,7);
r_Hm=0.047;

r_Ps=0.5;
p_A=0.004;

V_U=M1(:,8);

N_cohort=30;

% % simulate: initial condition
K_H=24000;
x0(4)=K_H; 
K_B=5e3;
x0(9)=K_B;

% initial condition for number of samples
y0=zeros(length(x0)*N_cohort,1);
for ind=1:length(x0)
    y0(1+(ind-1)*N_cohort:ind*N_cohort)=x0(ind)*ones(N_cohort,1);
end



% % simulate

% loop


[tpts,out] = ode15s(@(t,x)model0(t,x,N_cohort,dose_U,alpha_1,aT_APC,aTm_APC,...
   aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s, r_Ps, r_Ha,r_Hm,r_Ba,r_Bm,...
   V_U,h_m,h_n,b_n,b_m),(0:1e-2:tend),y0,options);


ind4=find(tpts==28);
ind8=find(tpts==54);
ind12=find(tpts==84);
ind26=find(tpts==182);
drug_TNFi = out(:,15*N_cohort+1:16*N_cohort);

ADA_sol = out(:,16*N_cohort+1:17*N_cohort);

% take out values at week=4,8,12,26
sampleADAL = drug_TNFi([ind4,ind8,ind12,ind26],:);
sampleANTIB=ADA_sol([ind4,ind8,ind12,ind26],:);

sampleANTIB=sampleANTIB/0.012;
sampleANTIB=sampleANTIB';
sampleADAL=sampleADAL';

D1_sol=out(:,2*N_cohort+1:3*N_cohort);
Ha_sol=out(:,4*N_cohort+1:5*N_cohort); Hm_sol = out(:,5*N_cohort+1:6*N_cohort); 
Hma_sol=out(:,6*N_cohort+1:7*N_cohort); Heff_sol = out(:,7*N_cohort+1:8*N_cohort);
% 
Ba_sol=out(:,9*N_cohort+1:10*N_cohort); Bm_sol=out(:,10*N_cohort+1:11*N_cohort); 
Bma_sol=out(:,11*N_cohort+1:12*N_cohort); 
% 
P_sol = out(:,17*N_cohort+1:18*N_cohort);
Ps_sol = out(:,18*N_cohort+1:19*N_cohort);

NK0_sol=out(:,12*N_cohort+1:13*N_cohort); 
NK1_sol=out(:,13*N_cohort+1:14*N_cohort);



cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
    max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];

% % % write to file
% 
writematrix(sampleANTIB,filename1,'Sheet','ADAtitre','WriteMode','append')
writematrix(sampleADAL,filename1,'Sheet','TNFi','WriteMode','append');
writematrix(cellmax,filename1,'Sheet','cellmax','WriteMode','append') 

fileID = fopen(filename2,'a');

% % write to file time series for the virtual patients
% % 
% 
for j = (1:6:N_cohort)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:10:end), D1_sol(1:10:end,j),Ha_sol(1:10:end,j),Hm_sol(1:10:end,j),...
        Hma_sol(1:10:end,j),Heff_sol(1:10:end,j), Ba_sol(1:10:end,j), ...
        Bm_sol(1:10:end,j),Bma_sol(1:10:end,j),P_sol(1:10:end,j), ...
        Ps_sol(1:10:end,j),NK0_sol(1:10:end,j),NK1_sol(1:10:end,j),...
        drug_TNFi(1:10:end,j),ADA_sol(1:10:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = model0(t,x,jj, dose_U, alpha_1, aT_APC,aTm_APC,...
   aB_T,aBm_T, KBm, KHm, k_101, aAPC_U, p_A, p, q, r, s,...
   r_Ps, r_Ha, r_Hm,r_Ba,r_Bm,V_U, h_m, h_n, b_n, b_m)

% variables
D0 = x(1:jj); Dast = x(jj+1:2*jj); D1 = x(2*jj+1:3*jj);
H0=x(3*jj+1:4*jj); Ha=x(4*jj+1:5*jj); Hm = x(5*jj+1:6*jj); 
Hma=x(6*jj+1:7*jj); He=x(7*jj+1:8*jj);
B0 =x(8*jj+1:9*jj); Ba = x(9*jj+1:10*jj); Bm = x(10*jj+1:11*jj); Bma = x(11*jj+1:12*jj); 
NK0 =x(12*jj+1:13*jj); NK1=x(13*jj+1:14*jj);
TNFa=x(14*jj+1:15*jj); U=x(15*jj+1:16*jj); A=x(16*jj+1:17*jj);
P =x(17*jj+1:18*jj); Ps=x(18*jj+1:19*jj);
 
% compute value xi(U)

 % sums
 Bsum=B0+Ba+Bm+Bma;
 Hsum=H0+Ha+Hm+Hma+He;

%   % APC

    dD0= mD*KD0*(1-D0/KD0)-aAPC_U.*D0.*U;

    dDast = aAPC_U.*D0.*U-alpha_1.*(1+TNFa./(TNFa+k_F)).*Dast;

%   % mature APC

    dD1 = alpha_1.*(1+TNFa./(TNFa+k_F)).*Dast-mD1.*D1./(1+TNFa./(TNFa+k_F));

%  % CD4 T cells

    dH0 = -aT_APC.*D1./(D1+h_n.*Hsum).*(TNFa./(TNFa+k_F)+1).*H0;

    dHa = aT_APC.*D1./(D1+h_n.*Hsum).*(TNFa./(TNFa+k_F)+1).*H0 ...
        + q.*r_Ha.*(D1./(D1+h_n.*Hsum)).*Ha - mHa.*Ha./(1+TNFa./(TNFa+k_F));

    dHm = (1-q).*p.*((D1./(D1+h_n.*Hsum)).*r_Ha.*Ha ...
        + (D1./(D1+h_m.*Hsum)).*r_Hm.*Hma)-aTm_APC.*D1./(D1+h_m.*Hsum).*(TNFa./(TNFa+k_F)+1).*Hm ...
        + aTm_APC.*(TNFa./(TNFa+k_F)+1).*Hm.*(1-Hm./KHm);

    dHma = aTm_APC.*D1./(D1+h_m.*Hsum).*(TNFa./(TNFa+k_F)+1).*Hm ...
        + q.*r_Hm.*(D1./(D1+h_m.*Hsum)).*Hma - mHa.*Hma./(1+TNFa./(TNFa+k_F));

    dHe = (1-p).*(1-q).*((D1./(D1+h_n.*Hsum)).*r_Ha.*Ha ...
        + (D1./(D1+h_m.*(Hsum))).*r_Hm.*Hma)- mHe.*He./(1+TNFa./(TNFa+k_F)) ;

    
  %  % B cells
    
    dB0 = -aB_T.*He./(He+b_n.*Bsum).*B0 ;

    dBa = aB_T.*He./(He+b_n.*Bsum).*B0 ...
        +s.*r_Ba.*(He./(He+b_n.*Bsum)).*Ba -mBa.*Ba -r_Ps.*Ba;

    dBm = r.*(1-s).*(r_Ba.*(He./(He+b_n.*Bsum)).*Ba ...
        + r_Bm.*(He./(He+b_m.*Bsum)).*Bma) - aBm_T.*He./(He+b_m.*Bsum).*Bm + aBm_T.*Bm.*(1-Bm./KBm);

    dBma = aBm_T.*He./(He+b_m.*Bsum).*Bm ...
        +s.*r_Bm.*(He./(He+b_m.*Bsum)).*Bma-mBa.*Bma;

   
    % short-lived plasma cells
	dPs = r_Ps.*Ba -mPs.*Ps;
	
	% long-lived plasma cells
	dP= (1-r).*(1-s).*(r_Ba.*(He./(He+b_n.*Bsum)).*Ba ...
        + r_Bm.*(He./(He+b_m.*Bsum)).*Bma)-mP*P;

    % NK cells
    
    dNK0 = mNK*(K_NK-NK0);
    
    dNK1= alpha_2.*(1+TNFa./(TNFa+k_F)).*(NK0-NK1)-mNK.*NK1;
    
    % TNFa
    
    dTNFa = k0_61+k1_61.*(D0+D1)+k3_61.*(Ha+Hma)+k5_61.*NK1-d_a.*TNFa-k_100.*TNFa.*U;
    

% adalimumab [anti-TNFa]
    U0 =dose_U*(1-exp(-kabsorpU*(floor(t/tau_U)+1)*tau_U))/(1-exp(-kabsorpU*tau_U))*exp(-kabsorpU*rem(t,tau_U));

    dU = kabsorpU.*U0./V_U-k_101.*A.*U-d_U.*U;
    
%  anti-drug antibody

    dA = p_A.*Ps + p_A.*P -d_A.*A;
    
    out=[dD0; dDast; dD1; dH0;dHa; dHm; dHma; dHe; ...
        dB0; dBa; dBm; dBma; dNK0; dNK1; dTNFa; dU; dA; dP; dPs];
end

end