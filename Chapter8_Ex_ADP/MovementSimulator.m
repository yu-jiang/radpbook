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
        
        
        
        theta=0*15/180*pi;%0.9453;
        
        
        theta1=0/180*pi;%0.9453;
        
        K
        Ko
        K0
    end
    
    % Public methods for interfacing
    methods
        function this = MovementSimulator(this)
            Initialize(this);
        end
        
        % Simulation for the NF
        function SimNF(this)
            [t,y1]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this.dt, this);
            [t,y2]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this.dt, this);
            [t,y3]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this.dt, this);
            [t,y4]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this.dt, this);
            [t,y5]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this.dt, this);
            
            figure(1)
            subplot(141)
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
            
            
            figure(4)
            subplot(441)
            plot(t,y1(:,3),'b-', ...
                t,y2(:,3),'b-', ...
                t,y3(:,3),'b-', ...
                t,y4(:,3),'b-', ...
                t,y5(:,3),'b-')
            xlabel('time (s)')
            ylabel('x-velocity (m/s)')
            title('A')
            axis([0 0.6 -0.5 0.5])
            subplot(445)
            plot(t,y1(:,4),'b-', ...
                t,y2(:,4),'b-', ...
                t,y3(:,4),'b-', ...
                t,y4(:,4),'b-', ...
                t,y5(:,4),'b-')
            
            axis([0 0.6 -.5 1.5])
            xlabel('time (s)')
            ylabel('y-velocity (m/s)')
            
            subplot(4,4,9)
            plot(t,y1(:,5),'b-', ...
                t,y2(:,5),'b-', ...
                t,y3(:,5),'b-', ...
                t,y4(:,5),'b-', ...
                t,y5(:,5),'b-')
            xlabel('time (s)')
            ylabel('x-endpoint force (N)')
            axis([0 0.6 -15 15])
            
            
            subplot(4,4,13)
            plot(t,y1(:,6),'b-', ...
                t,y2(:,6),'b-', ...
                t,y3(:,6),'b-', ...
                t,y4(:,6),'b-', ...
                t,y5(:,6),'b-')
            
            axis([0 0.6 -30 30])
            xlabel('time (s)')
            ylabel('y-endpoint force (N)')
            
            toc
        end
        
        %         function SimDF
        %         end
        %
        %         function SimAE
        %         end
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
            
            this.K = lqr(this.A0, this.B, this.Q, this.R);
            
            this.Ko = lqr(this.A,this.B,Q1,this.R);
            this.K0 = this.K;
            this.K0(1,1) = this.K(1,1) + 350;%230;
            this.K0(2,2) = this.K(1,1);
        end
    end
end

%%  ================= Local Functions ========================== %


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalFunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,Y] = ode_yu_NF(t0,tf,x0,dt,this)
h = dt;
t = t0:dt:tf;
y = x0;
Y = [];
for clock = t
    y = y + motor2DNF(y,this.A0,this.B,this.Q,this.R, ...
        this.c1,this.c2,this.K,h)*h;
    Y = [Y y];
end
Y = Y';
end

function dX = motor2DNF(X,A0,B,Q,R,c1,c2,K,dt)
x = X(1:6);
u = -K*x;
w = randn(2,1)*sqrt(dt);
M = [c1*u(1) c2*u(2); 
     -c2*u(1) c1*u(2)];
v = M*w; % control dependent noise
dx = A0*x + B*u + B*v./dt; %6
dIq = x'*(Q+K'*R*K)*x; %1
dxxi = kron(x',v'*R)'./dt; %12;
duu = kron(v,v)./dt;         % 4
dX =[dx;  % \dot x
    dIq; % reward
    dxxi;
    duu];
end

% ========
function h = fCirc(center,r,N,color)
theta = linspace(0,2*pi,N);
rho = ones(1,N)*r;
[X,Y] = pol2cart(theta,rho);
X = X+center(1);
Y = Y+center(2);
h = fill(X,Y,color);
axis square;
end