global Kadp
disp('Simulating the 1st trial...')
para;Kadp=K;
[t1,y1]=ode_yu_lrn(0,2,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
disp('Simulating the 2rd trial...')
Kadp=K0;  % after the first trial, the CNS directly increase the stiffness 
          % to formulate a new control policy which is also the initial 
          % control policy to perform ADP

%% learning
disp('Simulating the 2nd trial...')
ADPlearning; t2=t_save(:); y2=x_save(:,1:6);
disp('Simulating the 3rd trial...')
ADPlearning; t3=t_save(:); y3=x_save(:,1:6);
disp('Simulating the 4th trial...')
ADPlearning; t4=t_save(:); y4=x_save(:,1:6);
disp('Simulating the 5th trial...')
ADPlearning; t5=t_save(:); y5=x_save(:,1:6);    
    
%%
figure(4)
subplot(4,4,2)
plot(t1,y1(:,3),'b-', ...
    t2,y2(:,3),'b-', ...
    t3,y3(:,3),'b-', ...
    t4,y4(:,3),'b-', ...
    t5,y5(:,3),'b-')
xlabel('time (s)')
title('B')
% ylabel('x-velocity (m/s)')
axis([0 0.6 -0.5 .5])
subplot(4,4,6)
plot(t1,y1(:,4),'b-', ...
    t2,y2(:,4),'b-', ...
    t3,y3(:,4),'b-', ...
    t4,y4(:,4),'b-', ...
    t5,y5(:,4),'b-')
xlabel('time (s)')
% ylabel('y-velocity (m/s)')
axis([0 0.6 -0.5 1.5])

subplot(4,4,2+4+4)
plot(t1,y1(:,5),'b-', ...
    t2,y2(:,5),'b-', ...
    t3,y3(:,5),'b-', ...
    t4,y4(:,5),'b-', ...
    t5,y5(:,5),'b-')

xlabel('time (s)')
% ylabel('x-acceleration (m/s^2)')
%axis([0 0.6 -15 15])
axis([0 0.6 -20 20])

subplot(4,4,2+4*3)
plot(t1,y1(:,6),'b-', ...
    t2,y2(:,6),'b-', ...
    t3,y3(:,6),'b-', ...
    t4,y4(:,6),'b-', ...
    t5,y5(:,6),'b-')

axis([0 0.6 -15 20])
xlabel('time (s)')
% ylabel('y-acceleration (m/s^2)')



%%
figure(1)
subplot(142)
filledCircle([0,0],0.01,1000,'r')
hold on
% line([-0.03 -0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
% line([0.03 0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
plot(y1(:,1),y1(:,2), 'r-', ...
    y2(:,1),y2(:,2), 'g-', ...
    y3(:,1),y3(:,2), 'm-', ...
    y4(:,1),y4(:,2), 'k-', ...
    y5(:,1),y5(:,2), 'b-')
hold off
axis equal
axis([-0.12 0.08 -0.3 .1])
%axis([-0.14 0.06 -0.3 .1])
title('B')
xlabel('x-position (m)')


