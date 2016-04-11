%% A Simulink model that implements the linear ADP algorithm
%
% This example shows how to implement the linear ADP algorithm in Simulink.
% The plant dynamics are not known to the controller, yet the controller 
% learns the optimal performance via real-time data.
%
% Copyright 2016 Yu Jiang
mdl = 'LinearADPSimullink_R2015a';
%% Open the model
%
open_system(mdl);
%% Simulate the system
%
sim(mdl);
%% Clean up
close_system(mdl,0);