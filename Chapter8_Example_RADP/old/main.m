figure(1)
subplot(221)
draw_bg
x0=[0, 0.1];timecourse_NF(x0,1);disp('1')
x0=[0,-0.1];timecourse_NF(x0,5);disp('2')
x0=[0.1, 0];timecourse_NF(x0,3);disp('3')
x0=[-.1, 0];timecourse_NF(x0,7);disp('4')

x0=[1,1]/10*sqrt(2)/2;timecourse_NF(x0,2);  disp('5')
x0=[-1,-1]/10*sqrt(2)/2;timecourse_NF(x0,6);disp('6')
x0=[-1,1]/10*sqrt(2)/2;timecourse_NF(x0,8); disp('7')
x0=[1,-1]/10*sqrt(2)/2;timecourse_NF(x0,4); disp('8')

xlabel('x-position (m)')
ylabel('y-position (m)')
title('A', 'fontsize', 12)
hold off

%%
figure(1)
subplot(222)
draw_bg
hold on
x0=[0,.1];timecourse_VF(x0,1);disp('1')
x0=[0,-.1];timecourse_VF(x0,5);disp('2')
x0=[.1,0];timecourse_VF(x0,3);disp('3')
x0=[-.1,0];timecourse_VF(x0,7);disp('4')
x0=[1,1]/10*sqrt(2)/2;timecourse_VF(x0,2);disp('5')
x0=[-1,-1]/10*sqrt(2)/2;timecourse_VF(x0,6);disp('6')
x0=[-1,1]/10*sqrt(2)/2;timecourse_VF(x0,8);disp('7')
x0=[1,-1]/10*sqrt(2)/2;timecourse_VF(x0,4);disp('8')
xlabel('x-position (m)')
ylabel('y-position (m)')
title('B', 'fontsize', 12)
hold off

%%
figure(1)
subplot(223)
draw_bg
x0=[0,.1];timecourse_AL(x0,1);disp('1')
x0=[0,-.1];timecourse_AL(x0,5);disp('2')
x0=[.1,0];timecourse_AL(x0,3);disp('3')
x0=[-.1,0];timecourse_AL(x0,7);disp('4')
x0=[1,1]/10*sqrt(2)/2;timecourse_AL(x0,2);disp('5')
x0=[-1,-1]/10*sqrt(2)/2;timecourse_AL(x0,6);disp('6')
x0=[-1,1]/10*sqrt(2)/2;timecourse_AL(x0,8);disp('7')
x0=[1,-1]/10*sqrt(2)/2;timecourse_AL(x0,4);disp('8')
xlabel('x-position (m)')
ylabel('y-position (m)')
title('C', 'fontsize', 12)
hold off


%%
figure(1)
subplot(224)
draw_bg
x0=[0,.1];timecourse_AE(x0);disp('1')
x0=[0,-.1];timecourse_AE(x0);disp('2')
x0=[.1,0];timecourse_AE(x0);disp('3')
x0=[-.1,0];timecourse_AE(x0);disp('4')
x0=[1,1]/10*sqrt(2)/2;timecourse_AE(x0);disp('5')
x0=[-1,-1]/10*sqrt(2)/2;timecourse_AE(x0);disp('6')
x0=[-1,1]/10*sqrt(2)/2;timecourse_AE(x0);disp('7')
x0=[1,-1]/10*sqrt(2)/2;timecourse_AE(x0);disp('8')
xlabel('x-position (m)')
ylabel('y-position (m)')
title('D', 'fontsize', 12)
hold off
