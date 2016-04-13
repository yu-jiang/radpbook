function dX = SystemWrapper(t,x,SysParams)
%% SystemWrapper
% Local function to Wrap the system with externally
% specified integators for learning purpose

% Get local copies for parameters
K = SysParams.K(:)';
Q = SysParams.Q;

x1 = x(1);
sgm  =  x1.^[1 2 3]';

if SysParams.noiseFlag
 e  = ExplorationNoise(t);
else
 e = 0;
end

 u  = -1/2*K*sgm + e; 
 dx = LocalSystemKernel(x1,u,SysParams.F); 
 dZ  =  x1.^[2 3 4 5 6]';
 deZ = sgm*e;
 dQ  = sgm'*(Q + 1/4*K'*K)*sgm; 
 dX  = [dx;dZ;deZ;dQ]; % length 1 + 5 + 3 + 1 = 10
end

%% LocalSystemKernel
% Local function to implement polynomial system dynamics
function dx = LocalSystemKernel(x,u, F)
 sgm  =  x.^[1 2 3]';
 dx = F * sgm + u;
end