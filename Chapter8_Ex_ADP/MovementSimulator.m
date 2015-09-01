classdef MovementSimulator < handle
    
    properties
        tau = 0.05;
        m1 = 2;
        m2 = 2;
        d11 = 0;
        d12 = 13;
        d21 = -13;
        d22 = 0;
        
        c1 = 0.15/2;
        c2 = 0.05/2;
        dt = 0.005;
        
        in = 0;
        A0
        A
        B
        Q0 = [500 0; 0 1000];
        Q;
        R = diag([0.01,0.01]);
        
        theta = 0*15/180*pi;%0.9453;
        theta1 = 0/180*pi;%0.9453;
        
        K
        Ko
        K0
        K_
        % Fig_xy = figure;
        % Fig_tx = figure;
    end
    
    % Public methods for interfacing
    methods
        function this = MovementSimulator(this)
            Initialize(this);
        end
        
        function delete(this)
            delete(figure(1));
            delete(figure(2));
        end
        
        % Simulation for the NF
        function SimNF(this)
            [t1,y1] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t2,y2] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t3,y3] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t4,y4] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t5,y5] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            updateFigure(t1, t2, t3, t4, t5, y1, y2, y3, y4, y5, 'NF');
        end
        
        % Simulation for learning from unstable to stable during the
        %initial exposure in the DF
        function SimDF(this)
            reset(this);
            [t1,y1] = LocalSDESolver(0,2,[0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            this.K = this.K+[30,0,0,0,0,0;0,0,0,0,0,0];
            [t2,y2] = LocalSDESolver(0,2,[-0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            this.K = this.K+[30,0,0,0,0,0;0,0,0,0,0,0];
            [t3,y3] = LocalSDESolver(0,2,[0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            this.K = this.K+[30,0,0,0,0,0;0,0,0,0,0,0];
            [t4,y4] = LocalSDESolver(0,2,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            this.K = this.K+[30,0,0,0,0,0;0,0,0,0,0,0];
            [t5,y5] = LocalSDESolver(0,2,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            
            % Trim trajectories
            [t1,y1] = trimTraj(t1,y1);
            [t2,y2] = trimTraj(t2,y2);
            [t3,y3] = trimTraj(t3,y3);
            [t4,y4] = trimTraj(t4,y4);
            [t5,y5] = trimTraj(t5,y5);
            
            updateFigure(t1,t2,t3,t4,t5, ...
                y1, y2, y3, y4, y5, 'DF');
        end
        
        function SimMoveNLearn(this)
            x_save=[];t_save=[];
            [xn,un] = size(B);%size of x
            N=400; %length of the window, should be at least greater than xn^2
            NN=4;  %max iteration times
            T = 0.001;
            dt = 0.0005;%005;%0.0001; % period to data recording
            
            X=[0,-0.25,0,0,0,0,zeros(1,12+1+4)];
            Dxx=[];
            Iq=[];
            Ixu=[];
            Iuu=[];
            
            for i=1:N
                [t,X]=ode_yu_lrn((i-1)*T,i*T,X(end,:)',dt);
                Dxx=[Dxx; kron(X(end,1:6),X(end,1:6))-kron(X(1,1:6),X(1,1:6))];
                Iq=[Iq; X(end,6+1)-X(1,6+1)];
                Ixu=[Ixu; X(end,6+1+1:6+1+12)-X(1,6+1+1:6+1+12)];
                Iuu=[Iuu; X(end,end-3:end)-X(1,end-3:end)];
                x_save=[x_save;X];
                t_save=[t_save;t'];
            end           
            
            Dxx = Dxx(:,[1:6,8:12,15:18,22:24,29:30,36]); 
            Iuu = Iuu(:,[1,2,4]);
            
            % Learning
            Y = -Iq;
            X1 = [Dxx,-2*Ixu, Iuu];
            pp = inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
            P = [pp(1)  pp(2)     pp(3)     pp(4)   pp(5)    pp(6) ;
                0     pp(7)     pp(8)     pp(9)   pp(10)   pp(11) ;
                0         0     pp(12)    pp(13)  pp(14)   pp(15) ;
                0         0        0      pp(16)  pp(17)   pp(18) ;
                0         0        0          0   pp(19)   pp(20) ;
                0         0        0          0   0        pp(21)];
            P = (P+P')/2
            
            [t,X]=ode_yu_lrn(t(end),1,X(end,:)',dt);
            x_save=[x_save;X]; t_save=[t_save;t'];
            Kadp=[pp(22:2:32)';pp(23:2:33)']         
            
        end
        
        function SimAL(this)
        end
        
        function SimAE(this)
            
            function reset(this)
                this.K = this.K_;
            end
            
        end
        
        % Private methods for internal use
        methods (Access = private)
            function Initialize(this)
                
                this.A0 = [0 0 1        0      0      0;
                    0 0 0        1      0      0;
                    0 0 0        0     1/this.m1    0;
                    0 0 0        0      0      1/this.m2;
                    0 0 0        0      -1/this.tau 0;
                    0 0 0        0      0      -1/this.tau];
                
                this.A = this.A0 + [zeros(2,6);
                    150/this.m1  zeros(1,5);
                    zeros(3,6)];
                
                this.B=[0 0;
                    0 0;
                    0 0;
                    0 0;
                    1/this.tau 0;
                    0 1/this.tau];
                
                Qc = this.Q0;
                this.Q = blkdiag(Qc,0.01*Qc,0.00005*Qc);
                
                TM1 = [cos(this.theta1) -sin(this.theta1);
                    sin(this.theta1) cos(this.theta1)];
                
                Qc1 = this.Q0+ 1e4*[.03 0; 0 0];
                Q1 = blkdiag(Qc1,0.01*Qc1,0.00005*Qc1);
                this.R = TM1'*this.R*TM1;
                
                this.K_ = lqr(this.A0, this.B, this.Q, this.R);
                
                this.Ko = lqr(this.A,this.B,Q1,this.R);
                this.K = this.K_;
                this.K0(1,1) = this.K(1,1) + 350;%230;
                this.K0(2,2) = this.K(1,1);
            end
            
            
        end
    end
    
    %%  ================= Local Functions ========================== %
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LocalFunctions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [t,Y] = LocalSDESolver(t0,tf,x0,this, disturb)
    if disturb == true
        A = this.A;
    else
        A = this.A0;
    end
    dt = this.dt;
    t = t0:dt:tf;
    y = x0;
    Y = [];
    for clock = t
        y = y + motor2D(y,A,this.B,this.Q,this.R, ...
            this.c1,this.c2,this.K,dt)*dt;
        Y = [Y y];
    end
    Y = Y';
    end
    
    function dX = motor2D(X,A,B,Q,R,c1,c2,K,dt)
    x = X(1:6);
    u = -K*x;
    w = randn(2,1)*sqrt(dt);
    M = [c1*u(1) c2*u(2);
        -c2*u(1) c1*u(2)];
    v = M*w;                     % control dependent noise
    dx = A*x + B*u + B*v./dt;    %6
    dIq = x'*(Q+K'*R*K)*x;       %1
    dxxi = kron(x',v'*R)'./dt;   %12;
    duu = kron(v,v)./dt;         % 4
    dX =[dx;                     % \dot x
        dIq;                     % reward
        dxxi;
        duu];
    end
    
    %% ======== Utility Functions ===================================
    function h = fCirc(center,r,N,color)
    theta = linspace(0,2*pi,N);
    rho = ones(1,N)*r;
    [X,Y] = pol2cart(theta,rho);
    X = X+center(1);
    Y = Y+center(2);
    h = fill(X,Y,color);
    axis square;
    end
    
    % Trim Trajectories
    function [t,y] = trimTraj(t,y)
    ct = 1;
    while abs(y(ct,1))<0.03 && norm(y(ct,1:2))>0.007 && (ct<length(y))
        ct = ct +1;
    end
    y(ct+1:end,:) = []; t(ct+1:end) = [];
    end
    
    function updateFigure(t1, t2, t3, t4, t5, y1, y2, y3, y4 ,y5, expStage)
    figure(1)
    switch expStage
        case 'NF'
            subplot(1,4,1)
            fCirc([0,0],0.0075,1000,'r')
            hold on
            plot(y1(:,1),y1(:,2), 'b-', ...
                y2(:,1),y2(:,2), 'b-', ...
                y3(:,1),y3(:,2), 'b-', ...
                y4(:,1),y4(:,2), 'b-', ...
                y5(:,1),y5(:,2), 'b-')
            hold off
            axis equal
            axis([-0.05 0.05 -0.3 .05])
            xlabel('x-position (m)')
            ylabel('y-position (m)')
            title('A')
        case 'DF'
            subplot(1,4,2)
            fCirc([0,0],0.0075,1000,'r')
            hold on
            line([-0.03 -0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0]);
            line([0.03 0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0]);
            plot(y1(:,1),y1(:,2), 'b-', ...
                y2(:,1),y2(:,2), 'b-', ...
                y3(:,1),y3(:,2), 'b-', ...
                y4(:,1),y4(:,2), 'b-', ...
                y5(:,1),y5(:,2), 'b-')
            hold off
            axis equal
            axis([-0.05 0.05 -0.3 .05])
            xlabel('x-position (m)')
            ylabel('y-position (m)')
            title('B')
        case 'AL'
            title('C')
        case 'AE'
            title('D')
    end
    
    figure(2)
    switch expStage
        case 'NF'
            subplot(4, 4, 1)
            title('A')
            ylabel('y-endpoint force (N)')
            plot(t1,y1(:,3),'b-', ...
                t2,y2(:,3),'b-', ...
                t3,y3(:,3),'b-', ...
                t4,y4(:,3),'b-', ...
                t5,y5(:,3),'b-')
            xlabel('time (s)')
            ylabel('x-velocity (m/s)')
            axis([0 0.6 -0.5 0.5])
            subplot(4, 4, 4 + 1)
            plot(t1,y1(:,4),'b-', ...
                t2,y2(:,4),'b-', ...
                t3,y3(:,4),'b-', ...
                t4,y4(:,4),'b-', ...
                t5,y5(:,4),'b-')
            axis([0 0.6 -.5 1.5])
            xlabel('time (s)')
            ylabel('y-velocity (m/s)')
            subplot(4,4,8 + 1)
            plot(t1,y1(:,5),'b-', ...
                t2,y2(:,5),'b-', ...
                t3,y3(:,5),'b-', ...
                t4,y4(:,5),'b-', ...
                t5,y5(:,5),'b-')
            xlabel('time (s)')
            ylabel('x-endpoint force (N)')
            axis([0 0.6 -15 15])
            subplot(4,4,12 + 1)
            plot(t1,y1(:,6),'b-', ...
                t2,y2(:,6),'b-', ...
                t3,y3(:,6),'b-', ...
                t4,y4(:,6),'b-', ...
                t5,y5(:,6),'b-')
            axis([0 0.6 -30 30])
            xlabel('time (s)')
        case 'DF'
            subplot(4, 4, 2)
            title('B')
            plot(t1,y1(:,3),'b-', ...
                t2,y2(:,3),'b-', ...
                t3,y3(:,3),'b-', ...
                t4,y4(:,3),'b-', ...
                t5,y5(:,3),'b-')
            xlabel('time (s)')
            axis([0 0.6 -0.5 0.5])
            subplot(4, 4, 4 + 2)
            plot(t1,y1(:,4),'b-', ...
                t2,y2(:,4),'b-', ...
                t3,y3(:,4),'b-', ...
                t4,y4(:,4),'b-', ...
                t5,y5(:,4),'b-')
            axis([0 0.6 -.5 1.5])
            xlabel('time (s)')
            subplot(4,4,8 + 2)
            plot(t1,y1(:,5),'b-', ...
                t2,y2(:,5),'b-', ...
                t3,y3(:,5),'b-', ...
                t4,y4(:,5),'b-', ...
                t5,y5(:,5),'b-')
            xlabel('time (s)')
            axis([0 0.6 -15 15])
            subplot(4,4,12 + 2)
            plot(t1,y1(:,6),'b-', ...
                t2,y2(:,6),'b-', ...
                t3,y3(:,6),'b-', ...
                t4,y4(:,6),'b-', ...
                t5,y5(:,6),'b-')
            axis([0 0.6 -30 30])
            xlabel('time (s)')
        case 3
            title('C')
            figOffset = 3;
        case 4
            title('D')
            figOffset = 4;
    end
    
    
    end
