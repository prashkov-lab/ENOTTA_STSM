%%----------------------------------------------------------------------------------------
%% uses MATLAB
%% Code performs statistical analysis of the virtual patient cohorts
%% according to their status: ADA-negative, ADA-low, ADA-high
%% exports results to txt file
%%----------------------------------------------------------------------------------------

% D_X=10 mg/mL
fileID1=fopen("statistics_virtualpatientsLH15a.txt",'w');

fileID2=fopen("statistics_virtualpatientsLH15cumula.txt",'w');
fprintf(fileID1, "# scenario, cohort size, perc. ADA-positive [> 12] at W26; perc. low ADA [>12 and <150]; perc. high ADA [>150] \n");
fprintf(fileID2, "# scenario, cohort size, perc. ADA-positive [> 12] at W26; perc. low ADA [>12 and <150]; perc. high ADA [>150] \n");

formatSpec="%g %g %g %g %g \n";
formatSpec2="# scenario no. %g \n";

Ntot=550;
% % read data from noMTX cohorts, ADA-low, ADA-transient
M=readmatrix("simul_results/LHsimul_lowADAtitre15.xlsx",'Sheet','ADAtitre');
M1=readmatrix("simul_results/LHsimul_transADA15.xlsx",'Sheet','ADAtitre');

M=[M;M1];
Ncohort1=size(M,1);
% % number ADA-positive, > 12 AU/ml at W26
titreW26=M(:,end);
NposADA0=length(titreW26(titreW26>12));

% % find the maximum titre
maxtitre0=max(M,[],2);

% % number of ADA-neg, <=12 AU/ml
NADAneg0=length(maxtitre0(maxtitre0<=12));

% % number of low ADA, >12 AU/ml and <100 AU/ml
NlowADA0=length(maxtitre0(maxtitre0<=100 & maxtitre0>12));


% % number of high ADA, >100 AU/ml
NhighADA0=length(maxtitre0(maxtitre0>100));


% % extract indices of ADA-neg
indADAneg=find(maxtitre0<=12);

% % extract indices of high ADA
% % high ADA
indADAhigh=find(maxtitre0>100);


NposADA1=zeros(3,1);
NnoADA1=NposADA1;
NhighADA1=NposADA1;

% % scenario number
N_scen=31;

out1=zeros(N_scen,5);
fprintf(fileID1, "#  [D_X=10mg/mL] \n");
% 
for ii= 1:N_scen
filename1="simul_results/simul"+ii+"MTX10_lowADAtitre.xlsx";
% % originating from low ADA pre-MTX
A1=readmatrix(filename1,'Sheet','ADAtitre');

Ncohort1=size(A1,1)-(NhighADA0 +NADAneg0);
% %number ADA-positive, > 12 AU/ml at W26
titreW26=A1(:,end);
% number of ADA-pos coming from ADA-low
NposADA1(1)=length(titreW26(titreW26>12 & (maxtitre0<=100 & maxtitre0>12)));

% number of ADA-pos coming from ADA-high
NposADA1(2)=length(titreW26(titreW26>12 & (maxtitre0>100)));

% number of ADA-pos coming from ADA-neg
NposADA1(3)=length(titreW26(titreW26>12 & (maxtitre0<=12)));
% find the maximum titre
maxtitre=max(A1,[],2);

% number of ADA-neg, <=12 AU/ml
% coming from ADA-low
% NnoADA1(1)=length(maxtitre(maxtitre<=12 & (maxtitre0<=100 & maxtitre0>12)));
% %  coming from ADA-high
% NnoADA1(2)=length(maxtitre(maxtitre<=12 & (maxtitre0>100)));
% % coming from ADA-neg
% NnoADA1(3)=length(maxtitre(maxtitre<=12 & (maxtitre0<=12)));

% number of low ADA, >12 AU/ml and <100 AU/ml
NlowADA1(1)=length(maxtitre(maxtitre<=100 & maxtitre>12 & (maxtitre0<=100 & maxtitre0>12)));
NlowADA1(2)=length(maxtitre(maxtitre<=100 & maxtitre>12 & (maxtitre0>100 )));
NlowADA1(3)=length(maxtitre(maxtitre<=100 & maxtitre>12 & (maxtitre0<=12)));


% number of high ADA, >100 
NhighADA1(1)=length(maxtitre(maxtitre>100 & (maxtitre0<=100 & maxtitre0>12)));
NhighADA1(2)=length(maxtitre(maxtitre>100 & (maxtitre0>100)));
NhighADA1(3)=length(maxtitre(maxtitre>100 & (maxtitre0<=12)));


out=[ii, Ncohort1, NposADA1(1)/Ncohort1, ...%NnoADA1(1)/Ncohort1, 
    NlowADA1(1)/Ncohort1, NhighADA1(1)/Ncohort1];

fprintf(fileID1,formatSpec2,ii);
fprintf(fileID1,'# from low antibody (12,100) AU/ml \n');
% 
fprintf(fileID1, formatSpec,out);
% 

% 
% % originating from high ADA pre-MTX
% 
fprintf(fileID1,'# from high antibody (>100) AU/ml \n');

filename2="simul_results/simul"+ii+"MTX10_highADAtitre.xlsx";
A2=readmatrix(filename2,'Sheet','ADAtitre');
% 
Ncohort2=size(A2,1) + NhighADA0;
% 
% % number ADA-positive, > 12 AU at W26
titreW26=A2(:,end);
NposADA2=length(titreW26(titreW26>12));
% % find the maximum titre
maxtitre=max(A2,[],2);
% 
% % number of no ADA, <=12 
% NnoADA2=length(maxtitre(maxtitre<=12)) + NnoADA1(2);
% 
% % number of low ADA, >12 and <100
NlowADA2=length(maxtitre(maxtitre<=100 & maxtitre>12)) + NlowADA1(2);

% 
% % number of high ADA, >100 AU/ml
NhighADA2=length(maxtitre(maxtitre>100)) + NhighADA1(2);

% 
% 
out=[ii, Ncohort2, NposADA2/Ncohort2, ... % NnoADA2/Ncohort2, 
    NlowADA2/Ncohort2, NhighADA2/Ncohort2];

% 
fprintf(fileID1, formatSpec,out);

% overall distribution
fprintf(fileID1, "# overall distribution \n");
out1(ii,1:5)=[ii, Ntot, (NposADA2 + NposADA1(1))/Ntot, ... %(NnoADA1(1)+NnoADA2+Ntot-(Ncohort1+Ncohort2))/Ntot,
    (NlowADA1(1)+NlowADA2)/Ntot, (NhighADA1(1)+ NhighADA2)/Ntot];
fprintf(fileID2, formatSpec,out1(ii,:));
end
% 
fclose(fileID1);
fclose(fileID2);
% probabilities from study:
p= [0.25,...% 0.673, 
    0.269, 0.058];

ee=ones(N_scen,3);
outf=out1;

outf(:,2:4)=abs(p.*ee-out1(:,3:5));
outf(:,5)=sum(outf(:,2:4),2);
[minval,indexmin] = min(outf(:,5))
sortval=sortrows(outf,5);
sortval1=sortval(:,[1,5])

out_sort=zeros(16,4);
out_sort(1:16,1:2)=sortval1(1:16,:);
out_sort(1:15,3:4)=sortval1(17:end,:);
dlmwrite('sorted_a.txt',out_sort,'delimiter',';','precision',3)
