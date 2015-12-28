%SimuNF_nl.m
%tic
function []=timecourse_AE(x0)


dt = 0.005;
tf = 2;
K = [547.6615   -6.6753  100.8891  -19.1842    3.1623         0;
	8.1756  447.1638  -19.1842   96.0024         0    3.1623];
[t1,y1] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, false);
[t2,y2] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, false);
[t3,y3] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, false);
[t4,y4] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, false);
[t5,y5] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, false);

plot(y1(:,1)+x0(1),y1(:,2)+x0(2), 'b-', ...
	y2(:,1)+x0(1),y2(:,2)+x0(2), 'b-', ...
	y3(:,1)+x0(1),y3(:,2)+x0(2), 'b-', ...
	y4(:,1)+x0(1),y4(:,2)+x0(2), 'b-', ...
	y5(:,1)+x0(1),y5(:,2)+x0(2), 'b-', ...
	'Linewidth', 1.2)



