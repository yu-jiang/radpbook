classdef MovementSimulator < handle
    
    properties
        tau = 0.05; % (s)
        m1 = 2;        % kg
        m2 = 2;        % kg
        d11 = 0;
        d12 = 0;
        d21 = 0;
        d22 = 0;
        c1 = 0.15/2;
        c2 = 0.05/2;
        dt = 0.005;
%         
%         in = 0;
%         
        A0 = [0 0 1        0      0      0;
              0 0 0        1      0      0;
              0 0 0        0     1/m1    0;
              0 0 0        0      0      1/m2;
              0 0 0        0      -1/tau 0;
              0 0 0        0      0      -1/tau];
        
        A = A0+[0         0 0     0      0 0;
            0         0 0     0      0 0;
            150/m1    0 0     0      0 0;
            0         0 0     0      0 0;
            0         0 0     0      0 0;
            0         0 0     0      0 0];

        B = [0 0;
            0 0;
            0 0;
            0 0;
            1/tau 0;
            0 1/tau];
        
        Q0 = diag([500 1000]);  % weighting matrix
        R = diag([0.01,0.01]);  % weighting matrix 
        theta = 0*15/180*pi;    %0.9453;
        Qc = Q0;
        Q = blkdiag(Qc,0.01*Qc,0.00005*Qc);

        theta1 = 0/180*pi;%0.9453;
        TM1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
%         %Qc1 = Q0+ 250*[log(1+150) 0; 0 0];
%         Qc1 = Q0+ 1e4*[.03 0; 0 0];
%         Q1 = blkdiag(Qc1,0.01*Qc1,0.00005*Qc1);
       R = TM1'*R*TM1;
%         
%         %Q=10*diag([50,100,2,2,0.001,0.001]);R=diag([.01,.01]);
%         %Q1=10*diag([150,100,2,2,0.001,0.001]);
%         
%         K=lqr(A0,B,Q,R);
%         
%         Ko=lqr(A,B,Q1,R);
%         K0=K;
%         K0(1,1)=K(1,1)+350;%230;
%         K0(2,2)=K(1,1);
%         
%         
%         %K=K1
%         %eig(A-B*K);
        
    end
    
    methods
        function this = MovementSimulator(fieldDynamics)
            Initialize(this);
            
        end
        
        function this = Initialize(this)
        end
        
        function make_a_movement(this)
        end
        
    end
end