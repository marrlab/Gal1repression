clearvars;
clc;

%no degradation
simData(1,:) = getSimulation(100, 0); 

%super slow degradation
simData(2,:) = getSimulation(100, 1/24); 

%slow degradation
simData(3,:) = getSimulation(100, 0.2); 

%fast degradation
simData(4,:) = getSimulation(100, 0.5);

save('./Data/simData','simData')

getParametersSim(1,1,0,0)

getParametersSim(2,1,0,0)

load('scR1_model1_sim');
scR1 = scR;
BIC_m1 = [scR1(1).sol.BIC, scR1(2).sol.BIC, scR1(3).sol.BIC, scR1(4).sol.BIC];

load('scR1_model2_sim');
scR2 = scR;
BIC_m2 = [scR2(1).sol.BIC, scR2(2).sol.BIC, scR2(3).sol.BIC, scR2(4).sol.BIC];

est_tdelay = 10.^([scR2(1).sol.MS.par(2,1), scR2(2).sol.MS.par(2,1), scR2(3).sol.MS.par(2,1), scR2(4).sol.MS.par(2,1)]);
est_tdelay

BIC_diff = BIC_m2-BIC_m1;
BIC_diff
