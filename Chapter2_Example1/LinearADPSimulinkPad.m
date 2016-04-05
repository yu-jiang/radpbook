%% A Simulink model that implements the linear ADP algorithm
%
% This example shows how to implement the linear ADP algorithm in Simulink.
% The plant dynamics are not known to the controller, yet the controller 
% learns the optimal performance via real-time data.
%
% Copyright 2016 Yu Jiang
%% Open the model
%
open_system('LinearADPSimullink_R2015a');
%% Simulate the system
%
sim('LinearADPSimullink_R2015a');