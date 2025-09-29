%%-------------------------------------------------------------------------------
%% use MATLAB
%% Code computes the trajectories for the virtual patient cohorts 
%% with delayed co-medication of MTX
%% and exports the results
%%----------------------------------------------------------------------------------


function plots_MTX_delayed

tend=280;

%%% Set error tolerances for the integration
options = odeset('RelTol',1e-9,'AbsTol',1e-6,'MaxStep',1/30);

alpha_1=0.025;

K_NK=1.2e5; %from DATA,  MTX+ group

d_a=161;
d_U=0.03;
d_A=2.07;

k_F = 9.7e-6; %TNFa concentration = 9.7 pg/ml -> micro g/ML [Deveci]

% production of TNFa
kappa0 =1e-6; % rogers
kappa1 = 1e-9; % rogers pg/ml -> microg/ml
kappa3 = 1e-9; % rogers
kappa5 = 2.5e-8; % assumed

kabsorpU = 0.28; %  TYPICAL MONOCLONAL: ref. Ternant et al. ...
d_U = 0.03; %day^-1 (21 days 1/2-life)
tau_U = 14;

kabsorpX=8.64;
d_X=(2.08+3.33)/2; % average value 2.705

% MTX delay  
delayX=12; % 4, 8, 12 weeks

delayX=delayX*7;


KD0=5e4; % estimated from Dombrecht
K_B=5e3; % FROM DATA,
K_H = 24000; % from data

k_100=1010;
k_101=0.166;

mD=0.2;
mD1=0.75;


mNK = 0.05; % 1/2-life 14 days

mBa=0.032;

mHa=0.05;
mHe=0.113;

mP=0.01;
mPs=0.13;


% initial condition

x0 =zeros(21,1);
% D0 
x0(1) = KD0; 
% Dast = 
x0(2) =0; 
% D1 = 
x0(3) = 0;
% H0=
x0(4)=K_H; 
% Ha = Hm = Hma= He=0;
x0(5:8)=0; 

% B0 
x0(9)=K_B;
% Ba = Bm = Bma = 0;
x0(10:12)=0; 
% NK
x0(13)=K_NK;
x0(14)=.7*K_NK;
% TNFa
x0(15)= k_F; % mg/L 

%B10reg cells: normalised
K_10=1;

x0(20)=K_10;


alpha_2=mNK*(.7*K_NK)/(.3*K_NK); % MTX+ group

dose_X=10; %10mg/mL (conv_X=4/5) 
convX=0.8;

counter=0;
% titles for files:


patientno =('# virtual patient no. %g \n');
formatSpec=('%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n');


scenarios =[1,1,1,1,1;
    1,1,1,0,1;
    1,1,0,1,1;
    0,1,1,1,1;
    1,0,1,1,1;
    1,1,0,0,1;
    1,0,1,0,1;
    1,0,0,1,1;
    0,1,1,0,1;
    0,1,0,1,1;
    0,0,1,1,1;
    1,0,0,0,1;
    0,1,0,0,1;
    0,0,1,0,1;
    0,0,0,1,0;
    1,1,1,1,0;
    1,1,1,0,0;
    1,1,0,1,0;
    0,1,1,1,0;
    1,0,1,1,0;
    1,1,0,0,0;
    1,0,1,0,0;
    1,0,0,1,0;
    0,1,1,0,0;
    0,1,0,1,0;
    0,0,1,1,0;
    1,0,0,0,0;
    0,1,0,0,0;
    0,0,1,0,0;
    0,0,0,1,0;
    0,0,0,0,1];

% read parameter values!!!
% matrices MUST NOT include header row!!!

M2=readmatrix('param_MoA_MTX.csv','Range','A1:H290');

% % ADA-low subcohort <30 AU/mL
M1=readmatrix('LH15_lowADA1.csv','Range','A1:I130');

% read values from LH
aAPC_U=M1(1:10:end,1);

aB_T=M1(1:10:end,2);
aBm_T=0.001;

aT_APC=M1(1:10:end,3);
aTm_APC=0.001;

h_n=M1(1:10:end,4);
h_m=1;
b_n=M1(1:10:end,5);
b_m=0.05;

p=M1(1:10:end,6);
q=M1(1:10:end,7);
r=0.5;
s=0.5;

r_Ba=0.35;
r_Bm=0.4;
r_Ha=M1(1:10:end,8);
r_Hm=0.02;

r_Ps=0.9;

V_U=M1(1:10:end,9);

p_A=0.00269;

KBm=10;
KHm=10;

% number of trajectories
N_cohort=size(M1,1);

Ncoh0=N_cohort/10;

% initial condition for number of samples
y0=zeros(length(x0)*Ncoh0,1);
for ind=1:length(x0)
    y0(1+(ind-1)*Ncoh0:ind*Ncoh0)=x0(ind)*ones(Ncoh0,1);
end


% input additional parameters
% read from matrix

gamma_1=M2(counter+1:counter+Ncoh0,1);
gamma_2=M2(counter+1:counter+Ncoh0,2);
gamma_3=M2(counter+1:counter+Ncoh0,3);
gamma_4=M2(counter+1:counter+Ncoh0,4);
c_4=M2(counter+1:counter+Ncoh0,5);
a_10=M2(counter+1:counter+Ncoh0,6);
c_40=M2(counter+1:counter+Ncoh0,7);
c_41=M2(counter+1:counter+Ncoh0,8);



counter=counter+N_cohort;


for ind_scen=[9,13,24,28]
    scenario=scenarios(ind_scen,:);

    
    filename2="simul_delayedMTX"+delayX+"_scenario"+ind_scen+"ADAlow.txt";
    
    
% simulate the entire course following weeks

dose_U=50;
tau_X=7;

[tpts,out] = ode15s(@(t,x)model1(t,x,Ncoh0,dose_U,dose_X,tau_X,alpha_1,aT_APC,aTm_APC,...
    aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s,r_Ps,r_Ha, r_Hm,r_Ba,r_Bm,...
    V_U,h_m,h_n,b_n,b_m,gamma_4,c_4,gamma_1,gamma_2,gamma_3,a_10,c_40,c_41,scenario),(0:2e-1:tend),y0,options);


drug_TNFi = out(:,15*Ncoh0+1:16*Ncoh0);

ADA_sol = out(:,16*Ncoh0+1:17*Ncoh0);

D1_sol=out(:,2*Ncoh0+1:3*Ncoh0);
Ha_sol=out(:,4*Ncoh0+1:5*Ncoh0); Hm_sol = out(:,5*Ncoh0+1:6*Ncoh0); 
Hma_sol=out(:,6*Ncoh0+1:7*Ncoh0); Heff_sol = out(:,7*Ncoh0+1:8*Ncoh0);
% 
Ba_sol=out(:,9*Ncoh0+1:10*Ncoh0); Bm_sol=out(:,10*Ncoh0+1:11*Ncoh0); 
Bma_sol=out(:,11*Ncoh0+1:12*Ncoh0); 
% 
P_sol = out(:,17*Ncoh0+1:18*Ncoh0);
Ps_sol = out(:,18*Ncoh0+1:19*Ncoh0);

NK0_sol=out(:,12*Ncoh0+1:13*Ncoh0); 
NK1_sol=out(:,13*Ncoh0+1:14*Ncoh0);

fileID = fopen(filename2,'w');

fprintf(fileID,'# t D1 Ha Hm Hma Heff Ba Bm Bma P Ps NK0 NK1 TNFi ADA \n');
% % write to file time series for the virtual patients
% % 
% 
for j=(1:Ncoh0)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:3:end), D1_sol(1:3:end,j), Ha_sol(1:3:end,j), Hm_sol(1:3:end,j),...
        Hma_sol(1:3:end,j), Heff_sol(1:3:end,j), Ba_sol(1:3:end,j), ...
        Bm_sol(1:3:end,j), Bma_sol(1:3:end,j), P_sol(1:3:end,j), ...
        Ps_sol(1:3:end,j),NK0_sol(1:3:end,j),NK1_sol(1:3:end,j), ...
 drug_TNFi(1:3:end,j),ADA_sol(1:3:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);

end

%
%
% % low ADA: subcohort >30 AU/mL <100 AU/mL

M1=readmatrix('LH15_lowADA2.csv','Range','A1:F30');

% % read values from LH
aAPC_U=0.025;

aB_T=M1(1:5:end,1);
aBm_T=1.5;

aT_APC=M1(1:5:end,2);
aTm_APC=0.004;

h_n=0.005;
h_m=0.005;
b_n=M1(1:5:end,3);
b_m=2.5;

p=M1(1:5:end,4);
q=0.5;
r=0.5;
s=0.5;

r_Ba=0.35;
r_Bm=0.5;
r_Ha=M1(1:5:end,5);
r_Hm=0.047;

r_Ps=0.5;

V_U=M1(1:5:end,6);

KBm=100;
KHm=100;

% % number of patients is 30
N_cohort=size(M1,1);
Ncoh0=N_cohort/5;

gamma_1=M2(counter+1:counter+Ncoh0,1);
gamma_2=M2(counter+1:counter+Ncoh0,2);
gamma_3=M2(counter+1:counter+Ncoh0,3);
gamma_4=M2(counter+1:counter+Ncoh0,4);
c_4=M2(counter+1:counter+Ncoh0,5);
a_10=M2(counter+1:counter+Ncoh0,6);
c_40=M2(counter+1:counter+Ncoh0,7);
c_41=M2(counter+1:counter+Ncoh0,8);


counter=counter+N_cohort;

% initial condition for number of samples
y0=zeros(length(x0)*Ncoh0,1);
for ind=1:length(x0)
y0(1+(ind-1)*Ncoh0:ind*Ncoh0)=x0(ind)*ones(Ncoh0,1);
end


for ind_scen=[9,13,24,28]
    scenario=scenarios(ind_scen,:);

        
    filename2="simul_delayedMTX"+delayX+"_scenario"+ind_scen+"ADAlow.txt";
    


dose_U=50;
tau_X=7;

[tpts,out] = ode15s(@(t,x)model1(t,x,Ncoh0,dose_U,dose_X,tau_X,alpha_1,aT_APC,aTm_APC,...
    aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s,r_Ps,r_Ha, r_Hm,r_Ba,r_Bm,...
    V_U,h_m,h_n,b_n,b_m,gamma_4,c_4,gamma_1,gamma_2,gamma_3,a_10,c_40,c_41,scenario),(0:2e-1:tend),y0,options);

drug_TNFi = out(:,15*Ncoh0+1:16*Ncoh0);

ADA_sol = out(:,16*Ncoh0+1:17*Ncoh0);

D1_sol=out(:,2*Ncoh0+1:3*Ncoh0);
Ha_sol=out(:,4*Ncoh0+1:5*Ncoh0); Hm_sol = out(:,5*Ncoh0+1:6*Ncoh0); 
Hma_sol=out(:,6*Ncoh0+1:7*Ncoh0); Heff_sol = out(:,7*Ncoh0+1:8*Ncoh0);
% 
Ba_sol=out(:,9*Ncoh0+1:10*Ncoh0); Bm_sol=out(:,10*Ncoh0+1:11*Ncoh0); 
Bma_sol=out(:,11*Ncoh0+1:12*Ncoh0); 
% 
P_sol = out(:,17*Ncoh0+1:18*Ncoh0);
Ps_sol = out(:,18*Ncoh0+1:19*Ncoh0);

NK0_sol=out(:,12*Ncoh0+1:13*Ncoh0); 
NK1_sol=out(:,13*Ncoh0+1:14*Ncoh0);

% 
%

fileID = fopen(filename2,'a');

for j = (1:Ncoh0)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:3:end), D1_sol(1:3:end,j), Ha_sol(1:3:end,j), Hm_sol(1:3:end,j),...
        Hma_sol(1:3:end,j), Heff_sol(1:3:end,j), Ba_sol(1:3:end,j), ...
        Bm_sol(1:3:end,j), Bma_sol(1:3:end,j), P_sol(1:3:end,j), ...
        Ps_sol(1:3:end,j),NK0_sol(1:3:end,j),NK1_sol(1:3:end,j), ...
 drug_TNFi(1:3:end,j),ADA_sol(1:3:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);

end

counter=counter+30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % % high ADA group
% % read parameter values from Latin hypercube
 

% % cohort < 1000 AU/mL


M1=readmatrix('LH15_highADA1.csv','Range','A1:H70');
% 

aAPC_U=0.025;

aB_T=1;
aBm_T=1.5;

aT_APC=M1(1:7:end,1);
aTm_APC=0.004;

h_n=0.01;
h_m=0.01;
b_n=0.8;
b_m=M1(1:7:end,2);

p=M1(1:7:end,3);
q=M1(1:7:end,4);
r=M1(1:7:end,5);
s=M1(1:7:end,6);

r_Ba=0.35;
r_Bm=0.5;

r_Ha=M1(1:7:end,7);
r_Hm=0.047;

r_Ps=0.5;
p_A=0.004;

V_U=M1(1:7:end,8);



% 
% % number of virtual patients is 70
N_cohort=size(M1,1);
Ncoh0=N_cohort/7;
% 

gamma_1=M2(counter+1:counter+Ncoh0,1);
gamma_2=M2(counter+1:counter+Ncoh0,2);
gamma_3=M2(counter+1:counter+Ncoh0,3);
gamma_4=M2(counter+1:counter+Ncoh0,4);
c_4=M2(counter+1:counter+Ncoh0,5);
a_10=M2(counter+1:counter+Ncoh0,6);
c_40=M2(counter+1:counter+Ncoh0,7);
c_41=M2(counter+1:counter+Ncoh0,8);

% % % simulate: initial condition
K_H=24000;
x0(4)=K_H; 
K_B=5e3;
x0(9)=K_B;

% % initial condition for number of samples
y0=zeros(length(x0)*Ncoh0,1);
for ind=1:length(x0)
    y0(1+(ind-1)*Ncoh0:ind*Ncoh0)=x0(ind)*ones(Ncoh0,1);
end

counter=counter+N_cohort;



for ind_scen=[9,13,24,28]
    scenario=scenarios(ind_scen,:);

    
    filename2="simul_delayedMTX"+delayX+"_scenario"+ind_scen+"ADAhigh.txt";
    
    
% simulate the entire course following weeks

dose_U=50;
tau_X=7;

[tpts,out] = ode15s(@(t,x)model1(t,x,Ncoh0,dose_U,dose_X,tau_X,alpha_1,aT_APC,aTm_APC,...
    aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s,r_Ps,r_Ha, r_Hm,r_Ba,r_Bm,...
    V_U,h_m,h_n,b_n,b_m,gamma_4,c_4,gamma_1,gamma_2,gamma_3,a_10,c_40,c_41,scenario),(0:2e-1:tend),y0,options);

% ind4=find(tpts==28);
% ind8=find(tpts==54);
% ind12=find(tpts==84);
% ind26=find(tpts==182);

drug_TNFi = out(:,15*Ncoh0+1:16*Ncoh0);

ADA_sol = out(:,16*Ncoh0+1:17*Ncoh0);


D1_sol=out(:,2*Ncoh0+1:3*Ncoh0);
Ha_sol=out(:,4*Ncoh0+1:5*Ncoh0); Hm_sol = out(:,5*Ncoh0+1:6*Ncoh0); 
Hma_sol=out(:,6*Ncoh0+1:7*Ncoh0); Heff_sol = out(:,7*Ncoh0+1:8*Ncoh0);
% 
Ba_sol=out(:,9*Ncoh0+1:10*Ncoh0); Bm_sol=out(:,10*Ncoh0+1:11*Ncoh0); 
Bma_sol=out(:,11*Ncoh0+1:12*Ncoh0); 
% 
P_sol = out(:,17*Ncoh0+1:18*Ncoh0);
Ps_sol = out(:,18*Ncoh0+1:19*Ncoh0);

NK0_sol=out(:,12*Ncoh0+1:13*Ncoh0); 
NK1_sol=out(:,13*Ncoh0+1:14*Ncoh0);

fileID = fopen(filename2,'w');

fprintf(fileID,'# t D1 Ha Hm Hma Heff Ba Bm Bma P Ps NK0 NK1 TNFi ADA \n');
% % write to file time series for the virtual patients
% % 
% 
for j=(1:Ncoh0)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:3:end), D1_sol(1:3:end,j), Ha_sol(1:3:end,j), Hm_sol(1:3:end,j),...
        Hma_sol(1:3:end,j), Heff_sol(1:3:end,j), Ba_sol(1:3:end,j), ...
        Bm_sol(1:3:end,j), Bma_sol(1:3:end,j), P_sol(1:3:end,j), ...
        Ps_sol(1:3:end,j),NK0_sol(1:3:end,j),NK1_sol(1:3:end,j), ...
 drug_TNFi(1:3:end,j),ADA_sol(1:3:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);

% cellmax=[max(Ha_sol,[],1)',max(Hm_sol,[],1)',max(Hma_sol,[],1)',max(Heff_sol,[],1)',...
%     max(Ba_sol,[],1)',max(Bm_sol,[],1)',max(Bma_sol,[],1)',max(P_sol,[],1)',max(Ps_sol,[],1)'];
% 
% % write to file the values for W4,W8, W12, W26 for the virtual patients
% filename1="simul"+ind_scen+"MTX10_lowADAtitre.xlsx";
% writematrix(sampleANTIB,filename1,'Sheet','ADAtitre')
% writematrix(sampleADAL,filename1,'Sheet','TNFi')
% writematrix(cellmax,filename1,'Sheet','cellmax')
end

% % cohort > 1000 AU/mL


M1=readmatrix('LH15_highADA2.csv','Range','A1:H30');
% 
aAPC_U=0.025;

aB_T=1;
aBm_T=1.5;

aT_APC=M1(1:5:end,1);
aTm_APC=0.004;

h_n=0.01;
h_m=0.01;

b_n=0.1;
b_m=M1(1:5:end,2);

KBm=100;
KHm=100;

p=M1(1:5:end,3);
q=M1(1:5:end,4);
r=M1(1:5:end,5);
s=M1(1:5:end,6);

r_Ba=0.35;
r_Bm=0.5;

r_Ha=M1(1:5:end,7);
r_Hm=0.047;

r_Ps=0.5;
p_A=0.004;

V_U=M1(1:5:end,8);



N_cohort=size(M1,1);
Ncoh0=N_cohort/5;

gamma_1=M2(counter+1:counter+Ncoh0,1);
gamma_2=M2(counter+1:counter+Ncoh0,2);
gamma_3=M2(counter+1:counter+Ncoh0,3);
gamma_4=M2(counter+1:counter+Ncoh0,4);
c_4=M2(counter+1:counter+Ncoh0,5);
a_10=M2(counter+1:counter+Ncoh0,6);
c_40=M2(counter+1:counter+Ncoh0,7);
c_41=M2(counter+1:counter+Ncoh0,8);

% % simulate: initial condition
K_H=24000;
x0(4)=K_H; 
K_B=5e3;
x0(9)=K_B;

% initial condition for number of samples
y0=zeros(length(x0)*Ncoh0,1);
for ind=1:length(x0)
    y0(1+(ind-1)*Ncoh0:ind*Ncoh0)=x0(ind)*ones(Ncoh0,1);
end



for ind_scen=[9,13,24,28]
    scenario=scenarios(ind_scen,:);

        
    filename2="simul_delayedMTX"+delayX+"_scenario"+ind_scen+"ADAhigh.txt";
    


dose_U=50;
tau_X=7;

[tpts,out] = ode15s(@(t,x)model1(t,x,Ncoh0,dose_U,dose_X,tau_X,alpha_1,aT_APC,aTm_APC,...
    aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s,r_Ps,r_Ha, r_Hm,r_Ba,r_Bm,...
    V_U,h_m,h_n,b_n,b_m,gamma_4,c_4,gamma_1,gamma_2,gamma_3,a_10,c_40,c_41,scenario),(0:2e-1:tend),y0,options);

drug_TNFi = out(:,15*Ncoh0+1:16*Ncoh0);

ADA_sol = out(:,16*Ncoh0+1:17*Ncoh0);



D1_sol=out(:,2*Ncoh0+1:3*Ncoh0);
Ha_sol=out(:,4*Ncoh0+1:5*Ncoh0); Hm_sol = out(:,5*Ncoh0+1:6*Ncoh0); 
Hma_sol=out(:,6*Ncoh0+1:7*Ncoh0); Heff_sol = out(:,7*Ncoh0+1:8*Ncoh0);
% 
Ba_sol=out(:,9*Ncoh0+1:10*Ncoh0); Bm_sol=out(:,10*Ncoh0+1:11*Ncoh0); 
Bma_sol=out(:,11*Ncoh0+1:12*Ncoh0); 
% 
P_sol = out(:,17*Ncoh0+1:18*Ncoh0);
Ps_sol = out(:,18*Ncoh0+1:19*Ncoh0);

NK0_sol=out(:,12*Ncoh0+1:13*Ncoh0); 
NK1_sol=out(:,13*Ncoh0+1:14*Ncoh0);


fileID = fopen(filename2,'a');

for j = (1:Ncoh0)
    fprintf(fileID,patientno,j);
    
    outm=[tpts(1:3:end), D1_sol(1:3:end,j), Ha_sol(1:3:end,j), Hm_sol(1:3:end,j),...
        Hma_sol(1:3:end,j), Heff_sol(1:3:end,j), Ba_sol(1:3:end,j), ...
        Bm_sol(1:3:end,j), Bma_sol(1:3:end,j), P_sol(1:3:end,j), ...
        Ps_sol(1:3:end,j),NK0_sol(1:3:end,j),NK1_sol(1:3:end,j), ...
 drug_TNFi(1:3:end,j),ADA_sol(1:3:end,j)/0.012]';

fprintf(fileID,formatSpec, outm);
fprintf(fileID,'\n');

end

fclose(fileID);

end
% 



% model with MTX
function out = model1(t,x,jj,dose_U,dose_X,tau_X,alpha_1,aT_APC,aTm_APC,...
    aB_T,aBm_T,KBm,KHm,k_101,aAPC_U,p_A,p,q,r,s,r_Ps,r_Ha, r_Hm,r_Ba,r_Bm,...
    V_U,h_m,h_n,b_n,b_m,gamma_4,c_4,gamma_1,gamma_2,gamma_3,a_10,c_40,c_41,scenario)

% variables
D0 = x(1:jj); Dast = x(jj+1:2*jj); D1 = x(2*jj+1:3*jj);
H0=x(3*jj+1:4*jj); Ha=x(4*jj+1:5*jj); Hm = x(5*jj+1:6*jj); 
Hma=x(6*jj+1:7*jj); He=x(7*jj+1:8*jj);
B0 =x(8*jj+1:9*jj); Ba = x(9*jj+1:10*jj); Bm = x(10*jj+1:11*jj); Bma = x(11*jj+1:12*jj); 
NK0 =x(12*jj+1:13*jj); NK1=x(13*jj+1:14*jj);
TNFa=x(14*jj+1:15*jj); U=x(15*jj+1:16*jj); A=x(16*jj+1:17*jj);
P =x(17*jj+1:18*jj); Ps=x(18*jj+1:19*jj);
B10=x(19*jj+1:20*jj); MTX=x(20*jj+1:21*jj); 

 % sums
 Bsum=B0+Ba+Bm+Bma;
 Hsum=H0+Ha+Hm+Hma+He;
 
 % injection of methotrexate
doseX=dose_X*(t>=delayX);

% mode of action: methotrexate
%  on naive APC: D0
 phi1=gamma_1.*MTX.*D0;

% on activated T and activated T memory - Ha, Hma
 phi2= [gamma_4.*MTX.^2./(c_4.^2+MTX.^2).*Ha; gamma_4.*MTX.^2./(c_4.^2+MTX.^2).*Hma];
 
% on NK cells
 phi3=[gamma_2.*MTX.*NK0; gamma_3.*MTX.*NK1];

%  on Breg cells
 phi4=a_10.*MTX.*K_10;
 
% on activated B lympocytes
 phi5= [c_40.*MTX.*Ba;c_41.*MTX.*Bma];


dD0=mD*KD0*(1-D0/KD0)-aAPC_U.*D0.*U -scenario(1)*phi1;

dDast = aAPC_U.*D0.*U - alpha_1./(B10/K_10).*(1+TNFa./(TNFa+k_F)).*Dast;

dD1 = alpha_1./(B10/K_10).*(1+TNFa./(TNFa+k_F)).*Dast - mD1*D1./(1+TNFa./(TNFa+k_F));

% CD4 T cells

dH0 = -aT_APC./(B10/K_10).*D1./(D1+h_n.*Hsum).*(TNFa./(TNFa+k_F)+1).*H0;

dHa = aT_APC./(B10/K_10).*D1./(D1+h_n.*Hsum).*(TNFa./(TNFa+k_F)+1).*H0 ...
    +q.*r_Ha./(B10/K_10).*(D1./(D1+h_n.*Hsum)).*Ha ...
    - mHa*Ha./(1+TNFa./(TNFa+k_F))-scenario(2)*phi2(1:jj);

dHm = (1-q).*p./(B10/K_10).*((D1./(D1+h_n.*Hsum)).*r_Ha.*Ha ...
    + (D1./(D1+h_m.*Hsum)).*r_Hm.*Hma) ...
    -aTm_APC./(B10/K_10).*D1./(D1+h_m.*Hsum).*(TNFa./(TNFa+k_F)+1).*Hm ...
    + aTm_APC./(B10/K_10).*(TNFa./(TNFa+k_F)+1).*Hm.*(1-Hm./KHm);

dHma = aT_APC./(B10/K_10).*D1./(D1+h_m.*Hsum).*(TNFa./(TNFa+k_F)+1).*Hm ...
    +q.*r_Hm./(B10/K_10).*(D1./(D1+h_m.*Hsum)).*Hma ...
    - mHa.*Hma./(1+TNFa./(TNFa+k_F))-scenario(2)*phi2(jj+1:end);

dHe = (1-p).*(1-q)./(B10/K_10).*((D1./(D1+h_n.*Hsum)).*r_Ha.*Ha ...
    + (D1./(D1+h_m.*(Hsum))).*r_Hm.*Hma)- mHe.*He./(1+TNFa./(TNFa+k_F)) ;

    
%  % B cells
    dB0 = -aB_T.*He./(He+b_n.*Bsum).*B0 ;

    dBa = aB_T.*He./(He+b_n.*Bsum).*B0 ...
        +s.*r_Ba.*(He./(He+b_n.*Bsum)).*Ba -mBa.*Ba -r_Ps.*Ba -scenario(5)*phi5(1:jj);

    dBm = r.*(1-s).*(r_Ba.*(He./(He+b_n.*Bsum)).*Ba ...
        + r_Bm.*(He./(He+b_m.*Bsum)).*Bma) - aBm_T.*He./(He+b_m.*Bsum).*Bm ...
        + aBm_T.*Bm.*(1-Bm./KBm);

    dBma = aBm_T.*He./(He+b_m.*Bsum).*Bm ...
        +s.*r_Bm.*(He./(He+b_m.*Bsum)).*Bma-mBa.*Bma -scenario(5)*phi5(jj+1:end);

    
 %   % B10 cells
    
    dB10=mBa*(K_10+scenario(4).*phi4-B10);
    
    % plasma cells
    dPs = r_Ps.*Ba -mPs.*Ps;

	dP= (1-r).*(1-s).*(r_Ba.*(He./(He+b_n.*Bsum)).*Ba ...
        + r_Bm.*(He./(He+b_m.*Bsum)).*Bma)-mP*P;

   % NK cells
    
    dNK0 = mNK.*(K_NK-NK0)-scenario(3).*phi3(1:jj);
    
    dNK1= alpha_2.*(1+TNFa./(TNFa+k_F)).*(NK0-NK1)-mNK.*NK1...
        -scenario(3).*phi3(jj+1:end) ;
    
    % TNfa
    
    dTNFa = kappa0+kappa1.*(D0+D1)+kappa3.*(Ha+Hma)+kappa5.*NK1-d_a.*TNFa-k_100.*TNFa.*U;
    
% adalimumab [anti-TNFa]
    U0 =dose_U*(1-exp(-kabsorpU*ceil(t/tau_U)*tau_U))/(1-exp(-kabsorpU*tau_U))*exp(-kabsorpU.*rem(t,tau_U));

    dU = kabsorpU*U0./V_U-k_101.*A.*U-d_U.*U;

%  anti-drug antibody
    dA = p_A.*Ps+ p_A.*P-d_A.*A;

% methotrexate

    X0 =doseX*(1-exp(-kabsorpX*ceil(t/tau_X)*tau_X))/(1-exp(-kabsorpX*tau_X))*exp(-kabsorpX.*rem(t,tau_X));

    dMTX = kabsorpX*X0./(V_U*convX)-d_X.*MTX;

   
out=[dD0; dDast; dD1; dH0; dHa; dHm; dHma; dHe;...
     dB0; dBa; dBm; dBma; dNK0; dNK1; dTNFa; dU; ...
     dA; dP; dPs; dB10; dMTX];

end

% 
end
