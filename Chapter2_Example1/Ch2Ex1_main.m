%% A Simulink model that implements the linear ADP algorithm
%
% This example shows how to implement the linear ADP algorithm in Simulink.
% The plant dynamics are not known to the controller, yet the controller 
% learns the optimal performance via real-time data.
%
% Copyright 2016 Yu Jiang

%% Check your MATLAB version
mlver = ver('MATLAB');
rlStr = mlver.Release(5:7); % Get release number
if ismember(rlStr, {'14b', '15a', '15b'});
   mdl = ['LinearADPSimullink_R20', rlStr];
else
   error('The example is only supported in MATLAB R2014b--R2015b')
end

%% Open the model
open_system(mdl);

%% Simulate the system
%
sim(mdl);
%% Clean up
close_system(mdl,0);