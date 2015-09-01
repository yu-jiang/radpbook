%DFsimu_init.m
global Kadp
para;

Kadp=K;
[t1,y1]=ode_yu_lrn(0,2,[0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

Kadp=K+1*[30,0,0,0,0,0;0,0,0,0,0,0];

[t2,y2]=ode_yu_lrn(0,2,[-0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

Kadp=K+2*[30,0,0,0,0,0;0,0,0,0,0,0];

[t3,y3]=ode_yu_lrn(0,2,[0.001,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

Kdap=K+3*[30,0,0,0,0,0;0,0,0,0,0,0];

[t4,y4]=ode_yu_lrn(0,2,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

Kadp = K+4*[30,0,0,0,0,0;0,0,0,0,0,0];

[t5,y5]=ode_yu_lrn(0,2,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

%y1x=y1(:,4); y1y=y1(:,6);
%y2x=y2(:,4); y2y=y2(:,6);

i=1;
while abs(y1(i,1))<0.03 & norm(y1(i,1:2))>0.007 & (i<length(y1))
    i=i+1;
end
y1(i+1:end,:)=[];t1(i+1:end)=[];

i=1;
while abs(y2(i,1))<0.03 & norm(y2(i,1:2))>0.007 & (i<length(y2))
    i=i+1;
end
y2(i+1:end,:)=[];t2(i+1:end)=[];

i=1;
while abs(y3(i,1))<0.03 & norm(y3(i,1:2))>0.007  & (i<length(y3))
    i=i+1;
end
y3(i+1:end,:)=[];t3(i+1:end)=[];

i=1;
while abs(y4(i,1))<0.03 & norm(y4(i,1:2))>0.007 & (i<length(y4))
    i=i+1;
end
y4(i+1:end,:)=[];t4(i+1:end)=[];

i=1;
while abs(y5(i,1))<0.03 & norm(y5(i,1:2))>0.007 & (i<length(y5))
    i=i+1;
end
y5(i+1:end,:)=[];t5(i+1:end)=[];


% learning

% Kadp=K0;
% 
% [t3,y3]=ode_yu_lrn(0,1,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt)
% 
% Kadp = [662.0996         0   84.7902         0    3.8938         0
%     0  316.2278         0   58.8120         0    2.8878];
% 
% [t4,y4]=ode_yu_lrn(0,1,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt)
% 
% 
% Kadp = [539.5074    0.0004   61.7851    0.0000    2.9824    0.0000
%     0.0003  315.8477    0.0000   58.7838    0.0000    2.8881]
% 
% [t5,y5]=ode_yu_lrn(0,1,[0,-.25,0,0,0,0,zeros(1,12+1+4)]',dt)


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
plot(t1,y1(:,5)*m1,'b-', ...
    t2,y2(:,5)*m1,'b-', ...
    t3,y3(:,5)*m1,'b-', ...
    t4,y4(:,5)*m1,'b-', ...
    t5,y5(:,5)*m1,'b-')

xlabel('time (s)')
% ylabel('x-acceleration (m/s^2)')
axis([0 0.6 -30 30])


subplot(4,4,2+4*3)
plot(t1,y1(:,6)*m2,'b-', ...
    t2,y2(:,6)*m2,'b-', ...
    t3,y3(:,6)*m2,'b-', ...
    t4,y4(:,6)*m2,'b-', ...
    t5,y5(:,6)*m2,'b-')

axis([0 0.6 -50 50])
xlabel('time (s)')
% ylabel('y-acceleration (m/s^2)')



%%
figure(1)
subplot(142)
filledCircle([0,0],0.0075,1000,'r')
hold on
line([-0.03 -0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
line([0.03 0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
plot(y1(:,1),y1(:,2), 'b-', ...
    y2(:,1),y2(:,2), 'b-', ...
    y3(:,1),y3(:,2), 'b-', ...
    y4(:,1),y4(:,2), 'b-', ...
    y5(:,1),y5(:,2), 'b-')
hold off
axis equal
axis([-0.05 0.05 -0.3 .05])
title('B')
xlabel('x-position (m)')
