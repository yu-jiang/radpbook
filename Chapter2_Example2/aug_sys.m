function dX = aug_sys(t,X,K,ifLearned,expl_noise_freq)
%% Augmented system, a wrapper for the actual dynamic system 

x = X(1:6);

if ~ifLearned;   % See if learning is stopped
	u = sum(sin(expl_noise_freq*t),2);
else
	u = -K*x;    % Exploitation
end

dx = actual_sys(x,u);
dxx = kron(x',x')';
dux = kron(x',u')';
dX  = [dx;dxx;dux];
end


function dx = actual_sys(x,u)
%% Actual dynamics of the desel engine. 
%  This is the system you can customize.

A = [-0.4125    -0.0248        0.0741     0.0089   0          0;
    101.5873    -7.2651        2.7608     2.8068   0          0;
    0.0704       0.0085       -0.0741    -0.0089   0          0.0200;
    0.0878       0.2672        0         -0.3674   0.0044     0.3962;
    -1.8414      0.0990        0          0       -0.0343    -0.0330;
    0            0             0       -359      187.5364   -87.0316];

B = [-0.0042  0.0064
    -1.0360  1.5849
    0.0042  0;
    0.1261  0;
    0      -0.0168;
    0       0];

dx = A*x+B*u;
end