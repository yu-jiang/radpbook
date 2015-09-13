% main
% generates the figures 
%
%#1 Preparing the Data 
para
tf=1.5;
dt=0.01;
[t1b,y1b]=ode_yu_DF(0,tf,[0,-.25,0,0,0,0,0,0]',dt);
[t2b,y2b]=ode_yu_DF(0,tf,[0,-.25,0,0,0,0,0,0]',dt);
[t3b,y3b]=ode_yu_DF(0,tf,[0,-.25,0,0,0,0,0,0]',dt);
[t4b,y4b]=ode_yu_DF(0,tf,[0,-.25,0,0,0,0,0,0]',dt);
[t5b,y5b]=ode_yu_DF(0,tf,[0,-.25,0,0,0,0,0,0]',dt);
%%
para_nd;
dt=0.01;
[t1a,y1a]=ode_yu_AL(0,1,[0,-.25,0,0,0,0,0,0]',dt);
[t2a,y2a]=ode_yu_AL(0,1,[0,-.25,0,0,0,0,0,0]',dt);
[t3a,y3a]=ode_yu_AL(0,1,[0,-.25,0,0,0,0,0,0]',dt);
[t4a,y4a]=ode_yu_AL(0,1,[0,-.25,0,0,0,0,0,0]',dt);
[t5a,y5a]=ode_yu_AL(0,1,[0,-.25,0,0,0,0,0,0]',dt);

%%
%#2 Plot the position trajectories
figure(1)
subplot(2,4,[1 5])
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(y1b(:,1),y1b(:,2), 'b-', ...
     y2b(:,1),y2b(:,2), 'b-', ...
     y3b(:,1),y3b(:,2), 'b-', ...
     y4b(:,1),y4b(:,2), 'b-', ...
     y5b(:,1),y5b(:,2), 'b-', ...
     'Linewidth', 1.2)
hold off
axis equal
axis([-0.075 0.075 -0.3 .075])
xlabel('x-position (m)')
ylabel('y-position (m)')
title('A', 'Fontsize', 12)
%%
% figure(1)
subplot(2,4,[1 5]+2)
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(y1a(:,1),y1a(:,2), 'b-', ...
     y2a(:,1),y2a(:,2), 'b-', ...
     y3a(:,1),y3a(:,2), 'b-', ...
     y4a(:,1),y4a(:,2), 'b-', ...
     y5a(:,1),y5a(:,2), 'b-', ...
     'Linewidth', 1.2)
hold off
axis equal
axis([-0.075 0.075 -0.3 .075])
xlabel('x-position (m)')
ylabel('y-position (m)')
title('C', 'Fontsize', 12)

%%
% #3 Plot the speed profiles 

subplot(2,4,2)
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(t1b,y1b(:,3), 'b-', ...
     t1b,y2b(:,3), 'b-', ...
     t1b,y3b(:,3), 'b-', ...
     t1b,y4b(:,3), 'b-', ...
     t1b,y5b(:,3), 'b-', ...
     'Linewidth', 1.2)
hold off
% axis equal
axis([0 tf -50 50])
ylabel('x-velocity (m/s)')
xlabel('Time (sec)')
title('B', 'Fontsize', 12)


subplot(2,4,6)
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(t1b,y1b(:,4), 'b-', ...
     t1b,y2b(:,4), 'b-', ...
     t1b,y3b(:,4), 'b-', ...
     t1b,y4b(:,4), 'b-', ...
     t1b,y5b(:,4), 'b-', ...
     'Linewidth', 1.2)
hold off
% axis equal
axis([0 tf -50 50])
ylabel('y-velocity (m/s)')
xlabel('Time (sec)')

%%
% #3 Plot the speed profiles 

subplot(2,4,4,'align')
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(t1a,y1a(:,3), 'b-', ...
     t1a,y2a(:,3), 'b-', ...
     t1a,y3a(:,3), 'b-', ...
     t1a,y4a(:,3), 'b-', ...
     t1a,y5a(:,3), 'b-', ...
     'Linewidth', 1.2)
hold off
% axis equal
axis([0 1 -1 1])
ylabel('x-velocity (m/s)')
xlabel('Time (sec)')
title('D', 'Fontsize', 12)


subplot(2,4,8,'align')
filledCircle([0,0],0.0075,1000,'r')
hold on
plot(t1a,y1a(:,4), 'b-', ...
     t1a,y2a(:,4), 'b-', ...
     t1a,y3a(:,4), 'b-', ...
     t1a,y4a(:,4), 'b-', ...
     t1a,y5a(:,4), 'b-', ...
     'Linewidth', 1.2)
hold off
% axis equal
axis([0 1 -1 1])
ylabel('y-velocity (m/s)')
xlabel('Time (sec)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #4 Calculate and draw the stiffness circle
%%

theta=0:pi/200:2*pi;
x=cos(theta);
y=sin(theta);

para;
Kxy1=K(1:2,1:2);
para_nd;
Kxy2=K(1:2,1:2);
xy1=(Kxy1)*[x;y];
xy2=(Kxy2)*[x;y];
%subplot(4,2,[7 8])
figure(2)
plot(xy1(1,:),xy1(2,:),'g',xy2(1,:),xy2(2,:),'r','Linewidth',2);
axis equal
axis([-1200 1200 -700 700]*1.2)
xlabel('x-stiffness (N m^{-1})')
ylabel('y-stiffness (N m^{-1})')



axis equal
%title('E', 'Fontsize', 12)