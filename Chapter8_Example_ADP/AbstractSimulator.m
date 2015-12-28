classdef (Abstract) AbstractSimulator < handle
% Abstract class for movement simulation. 
% Author: Yu Jiang
% Contact: yu.jiang@nyu.edu
% Copyright 2015 Yu Jiang

    properties (Constant)      
        % Constant simulation parameters
        tau = 0.05;   % Time constant
        m1 = 2;       % mass on x direction
        m2 = 2;       % mass on y direction
        
        c1 = 0.15/2;  % noise scale
        c2 = 0.05/2;  % noise scale
        dt_ = 0.005;   % sample time for learning
        
        Q0 = [500 0; 0 1000];   % Initial weighting matrices
        R = diag([0.01,0.01]);  % Initial weighting matrices
        
        % The null filed dynamics
        A0 = [zeros(2) eye(2) zeros(2);
              zeros(2,4) diag([1/AbstractSimulator.m1 1/AbstractSimulator.m2])
              zeros(2,4) -diag([1/AbstractSimulator.tau 1/AbstractSimulator.tau])];
        B =  [zeros(4,2); diag([1/AbstractSimulator.tau 1/AbstractSimulator.tau])];
    end
    
    properties
        % Variables 
        dt % Sample time for simulation               
        A  % Dynamics with Force field
    end
    
    
    properties        
        % Others
        fig1;
        fig2;
        fig3;
    end
end