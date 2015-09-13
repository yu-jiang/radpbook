%SimuNF_nl.m
%tic
function []=timecourse_VF(x0,i)
para
tf=1.5;
[t1,y1]=ode_yu_VF(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t2,y2]=ode_yu_VF(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t3,y3]=ode_yu_VF(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t4,y4]=ode_yu_VF(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);
[t5,y5]=ode_yu_VF(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt);


%%
figure(1)
subplot(222)
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


figure(2)
subplot(8,3,(i-1)*3+2)
plot(t1,sqrt(y1(:,3).^2+y1(:,4).^2), 'b-', ...
     t2,sqrt(y2(:,3).^2+y2(:,4).^2), 'b-', ...
     t3,sqrt(y3(:,3).^2+y3(:,4).^2), 'b-', ...
     t4,sqrt(y4(:,3).^2+y4(:,4).^2), 'b-', ...
     t5,sqrt(y5(:,3).^2+y5(:,4).^2), 'b-', ...
    'Linewidth', 1.2)

figure(1)
%hold on
%hold off
%axis equal
%axis([-0.075 0.075 -0.3 .075])


end