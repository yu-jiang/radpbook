%SimuNF_nl.m
%tic
function [] = plotTimecourse(x0, str, i)


dt = 0.005;
tf = 2;
switch str
	case 'NF'
		isField = false;
		offset = 1;
		K = [632.4555         0   68.3772         0    3.1623         0
			0  632.4555         0   68.3772         0    3.1623];
	case 'VF'
		isField = true;
		offset = 2;
		K = [632.4555         0   68.3772         0    3.1623         0
			0  632.4555         0   68.3772         0    3.1623];
	case 'AL'
		isField = true;
		offset = 3;
		K = [547.6615   -6.6753  100.8891  -19.1842    3.1623         0;
			8.1756  447.1638  -19.1842   96.0024         0    3.1623];
	case 'AE'
		offset = nan;
		isField = false;
		K = [547.6615   -6.6753  100.8891  -19.1842    3.1623         0;
			8.1756  447.1638  -19.1842   96.0024         0    3.1623];
end

[t1,y1] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, isField);
[t2,y2] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, isField);
[t3,y3] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, isField);
[t4,y4] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, isField);
[t5,y5] = LocalSDESolver(0,tf,[-x0(1),-x0(2),0,0,0,0,0,0]',dt, K, isField);

figure(1)
plot(y1(:,1)+x0(1),y1(:,2)+x0(2), 'b-', ...
	y2(:,1)+x0(1),y2(:,2)+x0(2), 'b-', ...
	y3(:,1)+x0(1),y3(:,2)+x0(2), 'b-', ...
	y4(:,1)+x0(1),y4(:,2)+x0(2), 'b-', ...
	y5(:,1)+x0(1),y5(:,2)+x0(2), 'b-', ...
	'Linewidth', 1.2)

if ~isnan(offset)
	figure(2)
	subplot(8,3,(i-1)*3 + offset)
	plot(t1,sqrt(y1(:,3).^2+y1(:,4).^2), 'b-', ...
		t2,sqrt(y2(:,3).^2+y2(:,4).^2), 'b-', ...
		t3,sqrt(y3(:,3).^2+y3(:,4).^2), 'b-', ...
		t4,sqrt(y4(:,3).^2+y4(:,4).^2), 'b-', ...
		t5,sqrt(y5(:,3).^2+y5(:,4).^2), 'b-', ...
		'Linewidth', 1.2)
	figure(1)
end
end



