function dxx = FTSystemWrapper(t,x,Params,K)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FTSystemWrapper: System Dynamics with external integrators for learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create local copies of the parameters
[F,G,Q,R] = deal(Params.F,Params.G,Params.Q,Params.R); 

[e1,e2] = ExplNoise(t);
xm = BasisPlanMono(x(1),x(2));
u = K*xm + [e1;e2];

dx = F*xm + G*u;
dphi = [BasisQuadPlantMono(x(1),x(2));2*e1*R*xm;2*e2*R*xm];
dQ = xm'*(Q+K'*R*K)*xm;

dxx=[dx;dphi;dQ]; %size 2+ 34 +1 = 37
end