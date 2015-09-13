%function Kadp=ADPlearning(Kadp,iter_max)
clc
disp('[Simulating the online learning process]')
%global Kadp
para;
Kadp=Ko;
[K1,P1,E1]=lqr(A,B,Q,R);
P_save=[norm(Po-P1)];P1_save=[norm(Po-P1)];
K_save=[norm(Ko-K1)];K1_save=[norm(Ko-K1)];
iter_max=10;


disp(['[0-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
for num_iter=1:iter_max
    
    x_save=[];t_save=[];
    N=400; %length of the window, should be at least greater than xn^2
    NN=4;  %max iteration times
    T=.001;
    dt=0.0001;
    
    %X=[0,-.25,0,0,0,0,zeros(1,12+1+4)]; % initial conditions
    X=[-.25,0,0,zeros(1,5)];
    
    Dxx=[];
    Iq=[];
    Ixu=[];
    Iuu=[];
    
    for i=1:N
        [t,X]=ode_yu_lrn((i-1)*T,i*T,X(end,:)',dt,Kadp);
        Dxx=[Dxx; kron(X(end,1:3),X(end,1:3))-kron(X(1,1:3),X(1,1:3))];
        Iq=[Iq; X(end,3+1)-X(1,3+1)];
        Ixu=[Ixu; X(end,3+1+1:3+1+3)-X(1,3+1+1:3+1+3)];
        Iuu=[Iuu; X(end,end)-X(1,end)];
        x_save=[x_save;X];
        t_save=[t_save;t'];
    end
    
    %%
    Dxx=Dxx(:,[1,2,3,5,6,9]);%[1:6,8:12,15:18,22:24,29:30,36]); % Only the distinct columns left % ADP learning
    
    %% Learning
    Y=-Iq;
    X1=[Dxx,-2*Ixu, Iuu];
    pp=inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
    P=[pp(1)  pp(2)     pp(3);
        0     pp(4)     pp(5);
        0         0     pp(6)];
    
    P=(P+P')/2;
    Kadp=[pp(7:9)'];
    disp(['[',num2str(num_iter), '-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
    [K1,P1,E1]=lqr(A,B,Q,R);
    P_save=[P_save norm(P-P1)]
    K_save=[K_save norm(Kadp-K1)]
end
%%



%% with noise
Kadp=Ko;
P1_save=[norm(Po-P1)];
K1_save=[norm(Ko-K1)];
disp(['[0-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
for num_iter=1:iter_max
    
    x_save=[];t_save=[];
    N=400; %length of the window, should be at least greater than xn^2
    NN=4;  %max iteration times
    T=.001;
    dt=0.0001;
    
    %X=[0,-.25,0,0,0,0,zeros(1,12+1+4)]; % initial conditions
    X=[-.25,0,0,zeros(1,5)];
    
    Dxx=[];
    Iq=[];
    Ixu=[];
    Iuu=[];
    
    for i=1:N
        [t,X]=ode_yu_lrn_w((i-1)*T,i*T,X(end,:)',dt,Kadp);
        Dxx=[Dxx; kron(X(end,1:3),X(end,1:3))-kron(X(1,1:3),X(1,1:3))];
        Iq=[Iq; X(end,3+1)-X(1,3+1)];
        Ixu=[Ixu; X(end,3+1+1:3+1+3)-X(1,3+1+1:3+1+3)];
        Iuu=[Iuu; X(end,end)-X(1,end)];
        x_save=[x_save;X];
        t_save=[t_save;t'];
    end
    
    
    Dxx=Dxx(:,[1,2,3,5,6,9]);%[1:6,8:12,15:18,22:24,29:30,36]); % Only the distinct columns left % ADP learning
    
    Y=-Iq;
    X1=[Dxx,-2*Ixu, Iuu];
    pp=inv(X1'*X1)*X1'*Y;           % solve the equations in the least-squares sense
    P=[pp(1)  pp(2)     pp(3);
        0     pp(4)     pp(5);
        0         0     pp(6)];
    
    P=(P+P')/2;
    Kadp=[pp(7:9)'];
    disp(['[',num2str(num_iter), '-th step]', 'the feedback gain is K=[', num2str(Kadp), ']'])
    [K1,P1,E1]=lqr(A,B,Q,R);
    P1_save=[P1_save norm(P-P1)]
    K1_save=[K1_save norm(Kadp-K1)]
end
%%


figure(2)
subplot(211)
plot(0:iter_max, P_save,'-o',0:iter_max, P1_save,'-*')
axis([-0.5 iter_max+0.5 -5 35])
legend('Without measurement noise', 'With measurement noise')
xlabel('Number of iterations', 'Fontsize', 12)
ylabel('A', 'Fontsize', 12)
text(0,0,'||P_{1,k}-P_1||', 'FontSize', 12)

subplot(212)
plot(0:iter_max, K_save,'-^',0:iter_max, K1_save,'->')
axis([-0.5 iter_max+0.5 -5 25])
text(0,0,'||K_{1,k}-K_1||', 'FontSize', 12)
xlabel('Number of iterations', 'Fontsize', 12)
ylabel('B', 'Fontsize', 12)
legend('Without measurement noise', 'With measurement noise')
