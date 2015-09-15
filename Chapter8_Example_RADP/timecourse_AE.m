%SimuNF_nl.m
%tic
function []=timecourse_AE(x0)
para
tf=2;
[t1,y1]=ode_yu_AE(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t2,y2]=ode_yu_AE(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t3,y3]=ode_yu_AE(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t4,y4]=ode_yu_AE(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t5,y5]=ode_yu_AE(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);


%%
% figure(1)
% %subplot(142)
% 
% hold on
% draw_bg
% line([-0.03 -0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
% line([ 0.03  0.03],[-0.2 0.01], 'linewidth', 2, 'color', [0,0,0])
plot(y1(:,1)+x0(1),y1(:,2)+x0(2), 'b-', ...
     y2(:,1)+x0(1),y2(:,2)+x0(2), 'b-', ...
     y3(:,1)+x0(1),y3(:,2)+x0(2), 'b-', ...
     y4(:,1)+x0(1),y4(:,2)+x0(2), 'b-', ...
     y5(:,1)+x0(1),y5(:,2)+x0(2), 'b-', ...
     'Linewidth', 1.2)
%hold off
%axis equal
%axis([-0.075 0.075 -0.3 .075])


