close all; clear
%% Simulate the DF
M = DFSimulator();
simNF(M);
simDF(M);
simAL(M);
simAE(M);
showStiffness(M);
%% Simulate the VF
M = VFSimulator();
simNF(M);
simVF(M);
simAL(M);
simAE(M);
showStiffness(M);