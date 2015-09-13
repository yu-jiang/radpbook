% 1.1 simulation without learning or uncertain static dynamics
clc
para;
tf=0.65;
ymax=0.1;
[t,y1]=ode1(0,tf,[-.25,0,0]');
[t,y2]=ode1(0,tf,[-.25,0,0]');
[t,y3]=ode1(0,tf,[-.25,0,0]');
[t,y4]=ode1(0,tf,[-.25,0,0]');
[t,y5]=ode1(0,tf,[-.25,0,0]');

%%

figure(1)
axes('position', [0 0 1 1])
subplot(231)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,y1(:,1), 'b-', ...
     t,y2(:,1), 'b-', ...
     t,y3(:,1), 'b-', ...
     t,y4(:,1), 'b-', ...
     t,y5(:,1), 'b-')
%hold off
%axis equal
axis([0 tf -0.3 ymax])
xlabel('Time (s)', 'FontSize', 12)
ylabel('Distance (m)', 'FontSize', 12)
title('A', 'FontSize', 18)

%figure(2)
subplot(234)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,y1(:,2), 'b-', ...
     t,y2(:,2), 'b-', ...
     t,y3(:,2), 'b-', ...
     t,y4(:,2), 'b-', ...
     t,y5(:,2), 'b-')
%hold off
%axis equal
axis([0 tf -0.3 1.6])
xlabel('Time (s)', 'FontSize', 12)
ylabel('Velocity (m/s)', 'FontSize', 12)




%% 1.2 simulation without learning or uncertain static dynamics
tf1=2;

[t,yuc1]=ode2(0,tf1,[-.25,0,0]');
[t,yuc2]=ode2(0,tf1,[-.25,0,0]');
[t,yuc3]=ode2(0,tf1,[-.25,0,0]');
[t,yuc4]=ode2(0,tf1,[-.25,0,0]');
[t,yuc5]=ode2(0,tf1,[-.25,0,0]');

%% figure(4)
subplot(232)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,yuc1(:,1), 'b-', ...
     t,yuc2(:,1), 'b-', ...
     t,yuc3(:,1), 'b-', ...
     t,yuc4(:,1), 'b-', ...
     t,yuc5(:,1), 'b-')
%hold off
%axis equal
axis([0 tf1 -0.3 ymax])
xlabel('Time (s)', 'FontSize', 12)
%ylabel('Position (m)')
title('B', 'FontSize', 18)

%figure(5)
subplot(235)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,yuc1(:,2), 'b-', ...
     t,yuc2(:,2), 'b-', ...
     t,yuc3(:,2), 'b-', ...
     t,yuc4(:,2), 'b-', ...
     t,yuc5(:,2), 'b-')
%hold off
%axis equal
axis([0 tf1 -0.3 1.15])
xlabel('Time (s)', 'FontSize', 12)
%ylabel('Velocity (m/s)')


% figure(6)
% %subplot(141)
% %filledCircle([0,0],0.0075,1000,'r')
% %hold on
% plot(t,yuc1(:,3), 'b-', ...
%      t,yuc2(:,3), 'b-', ...
%      t,yuc3(:,3), 'b-', ...
%      t,yuc4(:,3), 'b-', ...
%      t,yuc5(:,3), 'b-')
% %hold off
% %axis equal
% axis([0 tf -5 30])
% xlabel('Time (s)')
% ylabel('Endpoint force (N)')


%% 1.3 simulation after learning with uncertain static dynamics
tf2=.8;
global Kadp
m0=1;  % mass
m=3;
b=10; % viscosity
c  = 0.05;
tau = 0.05;
A0=[    0    1          0;  0  -b/m0        1/m0;   0     0     -1/tau];
A=[    0    1          0;  0  -b/m        1/m;   0     0     -1/tau];
B=[    0;       0;   1/tau];
Q=diag([400,1,.001]);
R=0.01;
[Ko,Po,Eo]=lqr(A0,B,Q,R);[Ko1,Po,Eo]=lqr(A,B,Q,R);

Kadp=ADPlearning(Ko,10);
[t,y1_1_3]=ode3(0,tf2,[-.25,0,0]',Kadp);%P=lyap((A-B*Kadp)',Q+Kadp'*R*Kadp);Kadp=1/R*B'*P;
[t,y2_1_3]=ode3(0,tf2,[-.25,0,0]',Kadp);%P=lyap((A-B*Kadp)',Q+Kadp'*R*Kadp);Kadp=1/R*B'*P;
[t,y3_1_3]=ode3(0,tf2,[-.25,0,0]',Kadp);%P=lyap((A-B*Kadp)',Q+Kadp'*R*Kadp);Kadp=1/R*B'*P;
[t,y4_1_3]=ode3(0,tf2,[-.25,0,0]',Kadp);%P=lyap((A-B*Kadp)',Q+Kadp'*R*Kadp);Kadp=1/R*B'*P;
[t,y5_1_3]=ode3(0,tf2,[-.25,0,0]',Kadp);%P=lyap((A-B*Kadp)',Q+Kadp'*R*Kadp);Kadp=1/R*B'*P;

%%
%figure(7)
subplot(233)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,y1_1_3(:,1), 'b-', ...
     t,y2_1_3(:,1), 'b-', ...
     t,y3_1_3(:,1), 'b-', ...
     t,y4_1_3(:,1), 'b-', ...
     t,y5_1_3(:,1), 'b-')
%hold off
%axis equal
axis([0 tf2 -0.3 ymax])
xlabel('Time (s)', 'FontSize', 12)
%ylabel('Position (m)')
title('C', 'FontSize', 18)
% 
%figure(8)
subplot(236)
%filledCircle([0,0],0.0075,1000,'r')
%hold on
plot(t,y1_1_3(:,2), 'b-', ...
     t,y2_1_3(:,2), 'b-', ...
     t,y3_1_3(:,2), 'b-', ...
     t,y4_1_3(:,2), 'b-', ...
     t,y5_1_3(:,2), 'b-')
%hold off
%axis equal
axis([0 tf2 -0.3 1])
xlabel('Time (s)', 'FontSize', 12)
