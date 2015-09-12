close all; clear
%% Simulate the DF
MDF = DFSimulator();
simNF(MDF);
simDF(MDF);
simAL(MDF);
simAE(MDF);
showStiffness(MDF);
%% Simulate the VF
MVF = VFSimulator();
simNF(MVF);
simVF(MVF);
simAL(MVF);
simAE(MVF);
showStiffness(MVF);
%% Valide Fitts Law
