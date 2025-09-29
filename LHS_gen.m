## generate Latin hypercube samples
## for the virtual patient cohorts
## run in OCTAVE

pkg load stk
pkg load io
#### read values (mean) from spreadsheet DataParameters.xlsx
## export the values to csv files 

##
### Transient ADA
##
M=xlsread('DataParameters.xlsx','TransientADA','B2:AH2');
##
r_Ps_read=M(:,end-1);
q_m_read=M(:,end-6);
r_Ba_read=M(:,end-3);
r_Ha_read=M(:,end-2);
V_U_read=M(:,end);
##
lb=0.85;ub=1.15;
##
BOX=[r_Ps_read*lb,r_Ps_read*ub;...
 q_m_read*lb,0.99;...
 r_Ba_read*lb,r_Ba_read*ub;...
 r_Ha_read*lb,r_Ha_read*ub; 
 V_U_read*lb, V_U_read*ub];

BOX=BOX';
DIM=size(BOX,2);
N2=30;
XADA0 = stk_sampling_randomlhs (N2, DIM, BOX);
out=struct(XADA0).data;
##
csvwrite('LH15_transientADA.csv',out);

##
## low ADA, low subcohort, titre <30AU/ml
##

M=xlsread('DataParameters.xlsx','ADA-low','B2:AM2');
##
aAPC_U_read=M(:,1);
aB_T_read=M(:,2);
aT_APC_read=M(:,6);
hn_read=M(:,11);
bn_read=M(:,12);
p_m_read=M(:,end-9);
q_m_read=M(:,end-8);
r_Ha_read=M(:,end-3);
V_U_read=M(:,end);
##
lb=0.85;ub=1.15;
##
BOX=[aAPC_U_read*lb,aAPC_U_read*ub;...
 aB_T_read*lb,aB_T_read*ub;...
 aT_APC_read*lb,aT_APC_read*ub;...
 hn_read*lb,hn_read*ub; ...
 bn_read*lb, bn_read*ub; ...
  p_m_read*lb,p_m_read*ub;...
 q_m_read*lb,0.99;...
 r_Ha_read*lb,r_Ha_read*ub; 
 V_U_read*lb, V_U_read*ub];

BOX=BOX';
DIM=size(BOX,2);
N1=250;
XADAlow1 = stk_sampling_randomlhs (N1, DIM, BOX);
out=out=struct(XADAlow1 ).data;
##
csvwrite('LH15_lowADA1n.csv',out);


##
#### low ADA, high subcohort, titre >30AU/ml,<100AU/ml
##

M=xlsread('DataParameters.xlsx','ADA-low','B3:AM3');
##
bn_read=M(:,12);
aB_T_read=M(:,2);
aT_APC_read=M(:,6);
p_m_read=M(:,end-9);
r_Ha_read=M(:,end-3);
V_U_read=M(:,end);
##
##
##
BOX=[aB_T_read*lb,aB_T_read*ub; ...
 aT_APC_read*lb,aT_APC_read*ub;...
 bn_read*lb, bn_read*ub; 
 p_m_read*lb,p_m_read*ub;...
 r_Ha_read*lb,r_Ha_read*ub; V_U_read*lb, V_U_read*ub];

BOX=BOX';
DIM=size(BOX,2);
N2=30;
XADAlow2 = stk_sampling_randomlhs (N2, DIM, BOX);
out=struct(XADAlow2 ).data;

csvwrite('LH15_lowADA2.csv',out);


##
##### high ADA, <1000 AU
##
M=xlsread('DataParameters.xlsx','ADA-high','B2:AM2');

N3=70;

aT_APC_read=M(:,6);

bm_read=M(:,14);

p_m_read=M(:,end-9);
q_m_read=M(:,end-8);
r_m_read=M(:,end-7);
s_m_read=M(:,end-6);

r_Ha_read=M(:,end-3);
V_U_read=M(:,end);
##
BOX=[ aT_APC_read*lb,aT_APC_read*ub;...
 bm_read*lb, bm_read*ub;...
 p_m_read*lb,p_m_read*ub;...
 q_m_read*lb,q_m_read*ub; ...
 r_m_read*lb,r_m_read*ub; ...
 s_m_read*lb,s_m_read*ub; 
 r_Ha_read*lb,r_Ha_read*ub; 
 V_U_read*lb, V_U_read*ub];
 
 BOX=BOX';
DIM=size(BOX,2);

XADAhigh1 = stk_sampling_randomlhs (N3, DIM, BOX);
out=struct(XADAhigh1 ).data;
##

csvwrite('LH15_highADA1.csv',out);

##
##### high ADA, >1000 AU
##

M=xlsread('DataParameters.xlsx','ADA-high','B3:AM3');

N3=30;

aT_APC_read=M(:,6);

bm_read=M(:,14);

p_m_read=M(:,end-9);
q_m_read=M(:,end-8);
r_m_read=M(:,end-7);
s_m_read=M(:,end-6);

r_Ha_read=M(:,end-3);
V_U_read=M(:,end);
##
BOX=[ aT_APC_read*lb,aT_APC_read*ub;...
 bm_read*lb, bm_read*ub;...
 p_m_read*lb,p_m_read*ub;...
 q_m_read*lb,q_m_read*ub; ...
 r_m_read*lb,r_m_read*ub; ...
 s_m_read*lb,s_m_read*ub; 
 r_Ha_read*lb,r_Ha_read*ub; 
 V_U_read*lb, V_U_read*ub];
 
 BOX=BOX';
DIM=size(BOX,2);

XADAhigh2 = stk_sampling_randomlhs (N3, DIM, BOX);
out=struct(XADAhigh2 ).data;
##

csvwrite('LH15_highADA2.csv',out);


