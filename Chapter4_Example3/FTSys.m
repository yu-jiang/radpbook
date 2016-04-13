function dx = FTSys(x, Params, K)
%% --------------------------------------------------------------
% FTSys - The actual system dynamics
% ------------------------------------------------------------------------
Fc = Params.F + Params.G*K;
dx = Fc*BasisPlanMono(x(1),x(2));
end