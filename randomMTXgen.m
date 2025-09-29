

% MTX+ prameters
% gamma1-4, c4, a10, c40, c41
% MTX=10
avgs= [ 209.4, 0.323, 0.423, 0.704, 0.015, 9.393,0.133,0.893];
sdevs=[ 60, 0.026, 0.138, 0.102, 0.006, 0.88,0,0];

intervals=[0.15,0.15];
N_cohort=290;

out=zeros(N_cohort,8);

out(:,1:6)=avgs(1:6)+sdevs(1:6).* randn(N_cohort,6);

out(:,7:8)=avgs(7:8).*(1-intervals)+2*intervals.*rand(N_cohort,2);

writematrix(out,'param_MoA_MTX.csv','Delimiter',',')