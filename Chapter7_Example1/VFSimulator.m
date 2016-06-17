classdef VFSimulator < AbstractSimulator
% Class for movement simulation in Velocity-depedent field. 
% Author: Yu Jiang
% Contact: yu.jiang@nyu.edu
% Copyright 2015 Yu Jiang
    
    properties      
        Q;
        Q1        
        K
        Ko
        K_
        status = 0; % 0: unlearned. 1: learned
    end
    
    % Public methods for interfacing
    methods
        function this = VFSimulator()
            Initialize(this);
        end
        
        % Simulation for the NF
        function simNF(this)
            % if this.status == 1;
            this.reset();
            % end
            [t1,y1] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t2,y2] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t3,y3] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t4,y4] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            [t5,y5] = LocalSDESolver(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            updateFigure(this, t1, t2, t3, t4, t5, y1, y2, y3, y4, y5, 'NF');
        end
        
        % Simulation for learning from unstable to stable during the
        % initial exposure in the VF
        function simVF(this)
            this.reset();
            this.status = 0;
            [t1,y1] = LocalSDESolver(0,2,[0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',this,1);
            % after the first trial, the CNS directly increase the stiffness
            % to formulate a new control policy which is also the initial
            % control policy to perform ADP
            
            this.K = 3*this.K;
            disp('Simulating the 2nd trial...')
            [~,x_save,t_save] = simMoveNLearn(this);
            t2 = t_save(:); y2 = x_save(:,1:6);
            
            disp('Simulating the 3rd trial...')
            [~,x_save,t_save] = simMoveNLearn(this);
            t3 = t_save(:); y3 = x_save(:,1:6);
            
            disp('Simulating the 4th trial...')
            [~,x_save,t_save] = simMoveNLearn(this);
            t4 = t_save(:); y4 = x_save(:,1:6);
            
            disp('Simulating the 5th trial...')
            [~,x_save,t_save] = simMoveNLearn(this);
            t5 = t_save(:); y5 = x_save(:,1:6);
            
            updateFigure(this, t1,t2,t3,t4,t5, ...
                y1, y2, y3, y4, y5, 'VF');
        end
        
        % This simulates one movement and
        function [K, x_save, t_save] = simMoveNLearn(this)
            N = 400; %length of the window, should be at least greater than xn^2
            NN = 4;  %max iteration times
            T = 0.001;
            
            this.dt = 0.0005; % period to data recording
            x_save = [];t_save=[];
            X = [0,-0.25,0,0,0,0,zeros(1,12+1+4)];
            Dxx = [];
            Iq = [];
            Ixu = [];
            Iuu = [];
            
            for i=1:N
                [t,X] = LocalSDESolver((i-1)*T,i*T, X(end,:)',this, true);
                %= ode_yu_lrn((i-1)*T,i*T,,dt);
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
            X1 = [Dxx,-2*Ixu,-Iuu];
            pp = X1\Y;% inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
            P = [pp(1)  pp(2)     pp(3)     pp(4)   pp(5)    pp(6) ;
                0     pp(7)     pp(8)     pp(9)   pp(10)   pp(11) ;
                0         0     pp(12)    pp(13)  pp(14)   pp(15) ;
                0         0        0      pp(16)  pp(17)   pp(18) ;
                0         0        0          0   pp(19)   pp(20) ;
                0         0        0          0   0        pp(21)];
            P = (P+P')/2;
            
            this.K = [pp(22:2:32)';pp(23:2:33)'];
            this.dt = this.dt_;
            K = this.K;
            %% simulate the rest part of the movement
            [t,X] = LocalSDESolver(t(end),1, X(end,:)',this, true);
            x_save = [x_save;X];
            t_save = [t_save;t'];
        end
        
        function simAL(this)
            if this.status == 0;
                % Make 25 trials to learn
                for ct = 1:25
                    simMoveNLearn(this);
                end
                this.status = 1;
            end
            
            [t1,y1] = LocalSDESolver(0,1,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,true);
            [t2,y2] = LocalSDESolver(0,1,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,true);
            [t3,y3] = LocalSDESolver(0,1,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,true);
            [t4,y4] = LocalSDESolver(0,1,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,true);
            [t5,y5] = LocalSDESolver(0,1,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,true);
            updateFigure(this, t1, t2, t3, t4, t5, y1, y2, y3, y4, y5, 'AL');
        end
        
        function simAE(this)
            if this.status == 0;
                % Make 30 trials to learn
                for ct = 1:30
                    SimMoveNLearn(this);
                end
                this.status = 1;
            end
            
            [t1,y1] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,false);
            [t2,y2] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,false);
            [t3,y3] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,false);
            [t4,y4] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,false);
            [t5,y5] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,false);
            updateFigure(this, t1, t2, t3, t4, t5, y1, y2, y3, y4, y5, 'AE');
        end
        
        function reset(this)
            this.K = this.K_;
            this.status = 0;
        end
        
        function showStiffness(this)
            %main_show_stiffness.m
            
            tt = linspace(0, 2*pi, 200);
            x = cos(tt);
            y = sin(tt);
            
            Kxy1 = this.K_(1:2,1:2);
            Kxy2 = this.K(1:2,1:2);
            xy1 = (Kxy1)*[x;y];
            xy2 = (Kxy2)*[x;y];
            
            this.fig3.Visible = 'on';
            set(0, 'CurrentFigure', this.fig3);
            plot(xy1(1,:),xy1(2,:),'g',xy2(1,:),xy2(2,:),'r-.','Linewidth',2);
            legend('NF', 'VF')
            axis equal
            axis([-1200 1200 -1000 1000])
            xlabel('x component of stiffness (N m^{-1})')
            ylabel('y component of stiffness (N m^{-1})')
            
            figure1 = this.fig3;
            annotation(figure1,'textbox',...
                [0.312564059900167 0.709265175718848 0.166637271214642 0.135463258785942],...
                'String',{'Divergent force filed'},...
                'FitBoxToText','off',...
                'LineStyle','none', ...
                'FontSize', 12);
            
            % Create textbox
            annotation(figure1,'textbox',...
                [0.56838768718802 0.747603833865815 0.115472545757071 0.092651757188498],...
                'String',{'Movement','direction'},...
                'FitBoxToText','off',...
                'LineStyle','none', ...
                'FontSize', 12);
            
            % Create arrow
            annotation(figure1,'arrow',[0.365622119815668 0.410138248847926],...
                [0.734905511811023 0.661417322834646]);
            
            % Create arrow
            annotation(figure1,'arrow',[0.514592933947773 0.514592933947773],...
                [0.508186351706037 0.700787401574803]);
            
            % Create arrow
            annotation(figure1,'arrow',[0.513087557603687 0.61136712749616],...
                [0.509183727034121 0.509186351706037],'LineStyle',':');
            
            % Create arrow
            annotation(figure1,'arrow',[0.621197523239425 0.719477093131899],...
                [0.509789787640181 0.509792412312097],'LineStyle',':');
            
            % Create arrow
            annotation(figure1,'arrow',[0.515286870318463 0.396907216494845],...
                [0.50881076666116 0.51048951048951],'LineStyle',':');
            
            % Create arrow
            annotation(figure1,'arrow',[0.38703944763805 0.297250859106529],...
                [0.512493750344143 0.512820512820513],'LineStyle',':');
            
            % Create textbox
            annotation(figure1,'textbox',...
                [0.334234941041751 0.416995144424391 0.115472545757071 0.092651757188498],...
                'String',{'Force','direction'},...
                'FitBoxToText','off',...
                'LineStyle','none', ...
                'FontSize', 12);
            
            % Create arrow
            annotation(figure1,'arrow',[0.64952380952381 0.577235023041475],...
                [0.335839598997494 0.433338069583407]);
            
            % Create textbox
            annotation(figure1,'textbox',...
                [0.640959156653913 0.24812030075188 0.157136081441325 0.0803052901555278],...
                'String',{'Null filed'},...
                'FitBoxToText','off',...
                'LineStyle','none', ...
                'FontSize', 12);
            
        end
        
        function Tend = getPostLearningMovementDuration(this, r)
            assert(this.status == 1, 'The method can only be called when learning is finished');
            [t,y] = LocalSDESolver(0,2,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',this,0);
            ct = length(t);
            while norm(y(ct,1:2))<r
                ct = ct-1;
            end
            Tend = t(ct);
        end
        
    end
    
    
    % Private methods for internal use
    methods (Access = private)
        function Initialize(this)
            
            this.dt = this.dt_;            
            kai = 0.7;
            d11 = 13*kai;
            d12 = -18*kai;
            d21 = 18*kai;
            d22 = 13*kai;
            
            this.A = this.A0 + blkdiag(zeros(2), ...
                [d11/this.m1 d12/this.m1; d21/this.m2 d22/this.m2], ...
                zeros(2));
            theta1 = 144/180*pi;%0.9453;
            TM1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
            Qc1 = 1500*[cos(theta1)*cos(theta1) cos(theta1)*sin(theta1); 
				cos(theta1)*sin(theta1) sin(theta1)*sin(theta1)]+this.Q0;
            this.Q1 = blkdiag(Qc1,0.01*Qc1,0.00005*Qc1);
           
            Qc = this.Q0;
            this.Q = blkdiag(Qc,0.01*Qc,0.00005*Qc);
            
            this.K_ = lqr(this.A0, this.B, this.Q, this.R);
            this.Ko = lqr(this.A,this.B,this.Q1,this.R);
            this.K = this.K_;

            this.fig1 = figure('Visible', 'off');
            this.fig2 = figure('Visible', 'off');
            this.fig3 = figure('Visible', 'off');
        end
        
        
        function updateFigure(this, t1, t2, t3, t4, t5, y1, y2, y3, y4 ,y5, expStage)
            this.fig1.Visible = 'on';
            set(0,'CurrentFigure',this.fig1);
            switch expStage
                case 'NF'
                    subplot(1,4,1)
                    fCirc([0,0],0.01,1000,'r');
                    hold on
                    plot(y1(:,1),y1(:,2), 'b-', ...
                        y2(:,1),y2(:,2), 'b-', ...
                        y3(:,1),y3(:,2), 'b-', ...
                        y4(:,1),y4(:,2), 'b-', ...
                        y5(:,1),y5(:,2), 'b-')
                    hold off
                    axis equal
                    axis([-0.1 0.1 -0.3 0.1])
                    xlabel('x-position (m)')
                    ylabel('y-position (m)')
                    title('A')
                case 'VF'
                    subplot(1,4,2)
                    fCirc([0,0],0.01,1000,'r');
                    hold on
                    plot(y1(:,1),y1(:,2), 'r-', ...
                        y2(:,1),y2(:,2), 'g-', ...
                        y3(:,1),y3(:,2), 'm-', ...
                        y4(:,1),y4(:,2), 'k-', ...
                        y5(:,1),y5(:,2), 'b-')
                    hold off
                    axis equal
                    axis([-0.15 0.05 -0.3 .1])
                    xlabel('x-position (m)')
                    ylabel('y-position (m)')
                    title('B')
                case 'AL'
                    subplot(143)
                    fCirc([0,0],0.01,1000,'r');
                    hold on
                    plot(y1(:,1),y1(:,2), 'b-', ...
                        y2(:,1),y2(:,2), 'b-', ...
                        y3(:,1),y3(:,2), 'b-', ...
                        y4(:,1),y4(:,2), 'b-', ...
                        y5(:,1),y5(:,2), 'b-')
                    hold off
                    axis equal
                    axis([-0.1 0.1 -0.3 0.1])
                    title('C')
                    xlabel('x-position (m)')
                case 'AE'
                    this.fig1
                    subplot(144)
                    fCirc([0,0],0.01,1000,'r');
                    hold on
                    plot(y1(:,1),y1(:,2), 'b-', ...
                        y2(:,1),y2(:,2), 'b-', ...
                        y3(:,1),y3(:,2), 'b-', ...
                        y4(:,1),y4(:,2), 'b-', ...
                        y5(:,1),y5(:,2), 'b-')
                    hold off
                    axis equal
                    axis([-0.1 0.1 -0.3 0.1])
                    title('D')
                    xlabel('x-position (m)')
            end
            
            this.fig2.Visible = 'on';
            set(0,'CurrentFigure',this.fig2);
            switch expStage
                case 'NF'
                    subplot(4, 4, 1)
                    ylabel('y-endpoint force (N)')
                    plot(t1,y1(:,3),'b-', ...
                        t2,y2(:,3),'b-', ...
                        t3,y3(:,3),'b-', ...
                        t4,y4(:,3),'b-', ...
                        t5,y5(:,3),'b-')
                    xlabel('time (s)')
                    ylabel('x-velocity (m/s)')
                    title('A')
                    % axis([0 0.6 -0.5 0.5])
                    xlim([0 0.6])
                    subplot(4, 4, 4 + 1)
                    plot(t1,y1(:,4),'b-', ...
                        t2,y2(:,4),'b-', ...
                        t3,y3(:,4),'b-', ...
                        t4,y4(:,4),'b-', ...
                        t5,y5(:,4),'b-')
                    % axis([0 0.6 -.5 1.5])
                    xlim([0 0.6])
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
                    % axis([0 0.6 -15 15])
                    xlim([0 0.6])
                    subplot(4,4,12 + 1)
                    plot(t1,y1(:,6),'b-', ...
                        t2,y2(:,6),'b-', ...
                        t3,y3(:,6),'b-', ...
                        t4,y4(:,6),'b-', ...
                        t5,y5(:,6),'b-')
                    %axis([0 0.6 -30 30])
                    xlim([0 0.6])
                    xlabel('time (s)')
                    ylabel('y-endpoint force (N)')
                case 'VF'
                    subplot(4, 4, 2)
                    plot(t1,y1(:,3),'b-', ...
                        t2,y2(:,3),'b-', ...
                        t3,y3(:,3),'b-', ...
                        t4,y4(:,3),'b-', ...
                        t5,y5(:,3),'b-')
                    xlabel('time (s)')
                    title('B')
                    % axis([0 0.6 -0.5 0.5])
                    xlim([0 0.6]);
                    subplot(4, 4, 4 + 2)
                    plot(t1,y1(:,4),'b-', ...
                        t2,y2(:,4),'b-', ...
                        t3,y3(:,4),'b-', ...
                        t4,y4(:,4),'b-', ...
                        t5,y5(:,4),'b-')
                    %axis([0 0.6 -.5 1.5])
                    xlim([0 0.6]);
                    xlabel('time (s)')
                    subplot(4,4,8 + 2)
                    plot(t1,y1(:,5),'b-', ...
                        t2,y2(:,5),'b-', ...
                        t3,y3(:,5),'b-', ...
                        t4,y4(:,5),'b-', ...
                        t5,y5(:,5),'b-')
                    xlabel('time (s)')
                    %axis([0 0.6 -15 15])
                    xlim([0 0.6]);
                    subplot(4,4,12 + 2)
                    plot(t1,y1(:,6),'b-', ...
                        t2,y2(:,6),'b-', ...
                        t3,y3(:,6),'b-', ...
                        t4,y4(:,6),'b-', ...
                        t5,y5(:,6),'b-')
                    % axis([0 0.6 -30 30])
                    xlim([0 0.6]);
                    xlabel('time (s)')
                case 'AL'
                    subplot(4,4,3)
                    plot(t1,y1(:,3),'b-', ...
                        t2,y2(:,3),'b-', ...
                        t3,y3(:,3),'b-', ...
                        t4,y4(:,3),'b-', ...
                        t5,y5(:,3),'b-')
                    xlabel('time (s)')
                    title('C')
                    % axis([0 0.6 -0.5 .5])
                    xlim([0 0.6]);
                    subplot(4,4,3+4)
                    plot(t1,y1(:,4),'b-', ...
                        t2,y2(:,4),'b-', ...
                        t3,y3(:,4),'b-', ...
                        t4,y4(:,4),'b-', ...
                        t5,y5(:,4),'b-')
                    xlabel('time (s)')
                    % axis([0 0.6 -0.5 1.5])
                    xlim([0 0.6])
                    subplot(4,4,3+4+4)
                    plot(t1,y1(:,5),'b-', ...
                        t2,y2(:,5),'b-', ...
                        t3,y3(:,5),'b-', ...
                        t4,y4(:,5),'b-', ...
                        t5,y5(:,5),'b-')
                    xlabel('time (s)')
                    % axis([0 0.6 -15 15])
                    xlim([0 0.6])
                    subplot(4,4,3+4*3)
                    plot(t1,y1(:,6),'b-', ...
                        t2,y2(:,6),'b-', ...
                        t3,y3(:,6),'b-', ...
                        t4,y4(:,6),'b-', ...
                        t5,y5(:,6),'b-')
                    % axis([0 0.6 -30 30])
                    xlim([0 0.6])
                    xlabel('time (s)')
                    
                case 'AE'
                    subplot(4,4,4)
                    plot(t1,y1(:,3),'b-', ...
                        t2,y2(:,3),'b-', ...
                        t3,y3(:,3),'b-', ...
                        t4,y4(:,3),'b-', ...
                        t5,y5(:,3),'b-')
                    xlabel('time (s)')
                    title('D')
                    xlim([0 0.6])
                    % axis([0 0.6 -0.5 .5])
                    subplot(4,4,4+4)
                    plot(t1,y1(:,4),'b-', ...
                        t2,y2(:,4),'b-', ...
                        t3,y3(:,4),'b-', ...
                        t4,y4(:,4),'b-', ...
                        t5,y5(:,4),'b-')
                    xlabel('time (s)')
                    % axis([0 0.6 -0.5 1.5])
                    xlim([0 0.6])
                    subplot(4,4,4+4+4)
                    plot(t1,y1(:,5),'b-', ...
                        t2,y2(:,5),'b-', ...
                        t3,y3(:,5),'b-', ...
                        t4,y4(:,5),'b-', ...
                        t5,y5(:,5),'b-')
                    xlabel('time (s)')
                    % axis([0 0.6 -15 15])
                    xlim([0 0.6])
                    subplot(4,4,4+4*3)
                    plot(t1,y1(:,6),'b-', ...
                        t2,y2(:,6),'b-', ...
                        t3,y3(:,6),'b-', ...
                        t4,y4(:,6),'b-', ...
                        t5,y5(:,6),'b-')
                    % axis([0 0.6 -30 30])
                    xlim([0 0.6])
                    xlabel('time (s)')
            end
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
    y = y + motor2D(y,A,this.B,this.Q1,this.R, ...
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
dIq = x'*(Q +K'*R*K)*x;       %1
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

