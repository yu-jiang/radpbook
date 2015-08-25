%SimuAE.m

para
[t,y1]=ode_yu_AE(0,2,[0*0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y2]=ode_yu_AE(0,2,[0*0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y3]=ode_yu_AE(0,2,[0*0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y4]=ode_yu_AE(0,2,[0*0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y5]=ode_yu_AE(0,2,[0*0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);



figure(1)
subplot(144)
filledCircle([0,0],0.01,1000,'r');

hold on
plot(y1(:,1),y1(:,2), 'b-', ... 
     y2(:,1),y2(:,2), 'b-', ...
     y3(:,1),y3(:,2), 'b-', ...
     y4(:,1),y4(:,2), 'b-', ...
     y5(:,1),y5(:,2), 'b-')
hold off

axis equal
axis([-0.1 0.1 -0.3 .1]) 
title('D')
xlabel('x-position (m)')

%%
% figure(4)
% subplot(244)
% plot(t,y1(:,4),'b-', ...
%      t,y2(:,4),'b-', ...
%      t,y3(:,4),'b-', ...
%      t,y4(:,4),'b-', ...
%      t,y5(:,4),'b-')
% xlabel('time (s)')
% %ylabel('y-velocity (m/s)')
% axis([0 0.6 -0.5 1.5])
% 
% 
% subplot(248)
% plot(t,y1(:,6),'b-', ...
%     t,y2(:,6),'b-', ...
%     t,y3(:,6),'b-', ...
%     t,y4(:,6),'b-', ...
%     t,y5(:,6),'b-')
% 
% axis([0 0.6 -10 20])
% xlabel('time (s)')
% ylabel('y-acceleration (m/s^2)')


figure(4)
subplot(4,4,4)
plot(t,y1(:,3),'b-', ...
    t,y2(:,3),'b-', ...
    t,y3(:,3),'b-', ...
    t,y4(:,3),'b-', ...
    t,y5(:,3),'b-')
xlabel('time (s)')
% ylabel('x-velocity (m/s)')
title('D')
axis([0 0.6 -0.5 .5])
subplot(4,4,4+4)
plot(t,y1(:,4),'b-', ...
    t,y2(:,4),'b-', ...
    t,y3(:,4),'b-', ...
    t,y4(:,4),'b-', ...
    t,y5(:,4),'b-')
xlabel('time (s)')
% ylabel('y-velocity (m/s)')
axis([0 0.6 -0.5 1.5])

subplot(4,4,4+4+4)
plot(t,y1(:,5),'b-', ...
    t,y2(:,5),'b-', ...
    t,y3(:,5),'b-', ...
    t,y4(:,5),'b-', ...
    t,y5(:,5),'b-')

xlabel('time (s)')
% ylabel('x-acceleration (m/s^2)')
axis([0 0.6 -15 15])


subplot(4,4,4+4*3)
plot(t,y1(:,6),'b-', ...
    t,y2(:,6),'b-', ...
    t,y3(:,6),'b-', ...
    t,y4(:,6),'b-', ...
    t,y5(:,6),'b-')

axis([0 0.6 -30 30])
xlabel('time (s)')
% ylabel('y-acceleration (m/s^2)')


