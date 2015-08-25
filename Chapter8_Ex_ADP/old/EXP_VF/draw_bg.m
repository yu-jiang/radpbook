%draw_bg.m Provide the background of the figure1


close all
figure(1)
subplot(141)
filledCircle([0,0],0.01,1000,'r')
hold on
NFsimu % plot the NF simulation
hold off
axis equal
axis([-0.05 0.05 -0.3 .05])
xlabel('x-position (m)')
ylabel('y-position (m)')
title('A')

subplot(142)
filledCircle([0,0],0.0075,1000,'r')
axis equal
axis([-0.05 0.05 -0.3 .05])
title('B')
xlabel('x-position (m)')
ylabel('y-position (m)')

subplot(143)
filledCircle([0,0],0.0075,1000,'r')
axis equal
axis([-0.05 0.05 -0.3 .05]) 
title('C')
xlabel('x-position (m)')
ylabel('y-position (m)')


subplot(144)
filledCircle([0,0],0.0075,1000,'r')
axis equal
axis([-0.05 0.05 -0.3 .05]) 
title('D')
xlabel('x-position (m)')
ylabel('y-position (m)')