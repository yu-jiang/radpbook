%SimuNF.m
tic
para
[t,y1]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y2]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y3]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y4]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);
[t,y5]=ode_yu_NF(0,.7,[0.001*(rand-0.5),-.25,0,0,0,0,zeros(1,12+1+4)]',dt);

%%
figure(1)
subplot(141)
filledCircle([0,0],0.0075,1000,'r')
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