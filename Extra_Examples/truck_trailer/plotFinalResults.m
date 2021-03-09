%% Setup
p = get_default_truck_trailer_params();
p.noiseLevel = 0.01;
p.trailerWheelbase = 11;
p.forwardTarget = [20 2 0.3];
xx = [0 -5 0 0];
opt = odeset('Events', @obstacleEvents);
dt = 0.2;
p.forwardGain = [0 0 0];
p.feedbackGain = [0 0 0];

bp1 = [1 210 420 630 837 1043];
bp2 = [1 212 423 629 830 1031];
%% setup viz
figure(101)
ax(1) = subplot(121);
tp1 = truck_trailer_plot(p, ax(1));
tp1.updateFig(xx)
title('No learning: Starting Pose')
ax(2) = subplot(122);
tp2 = truck_trailer_plot(p, ax(2));
tp2.updateFig(xx)
title('Learning via ADP: Starting Pose')
set(gcf, 'Color', 'w');
export_fig("results"+num2str(0)+".png")
%%
for ct = 2:6
figure(101)
ax(1) = subplot(121);
tp1 = truck_trailer_plot(p, ax(1));
tp1.updateFig(ys1(bp1(ct),:));
title(['No learning: Trial #', num2str(ct-1)])
hold on
plot(-ys1(bp1(ct-1):bp1(ct), 2), ys1(bp1(ct-1):bp1(ct), 1), ...
    'k:', ...
    'LineWidth', 1.2);
hold off
ax(2) = subplot(122);
tp2 = truck_trailer_plot(p, ax(2));
tp2.updateFig(ys2(bp2(ct),:));
title(['Learning via ADP: Trial #', num2str(ct-1)])
hold on
plot(-ys2(bp2(ct-1):bp2(ct), 2), ys2(bp2(ct-1):bp2(ct), 1), ...
    'k:', ...
    'LineWidth', 1.2);
hold off
drawnow
set(gcf, 'Color', 'w');
export_fig("results"+num2str(ct-1)+".png")
end

%%
figure(301)
subplot(411)
plot(ts1, ys1(:,1), ts2, ys2(:,1), 'LineWidth', 1.2)
xlim([0 ts1(end)])
ylabel('(m)')
title('x: Longitudinal displacement')
legend('No learning', 'Learning via ADP')
grid on

subplot(412)
plot(ts1, ys1(:,2), ts2, ys2(:,2), 'LineWidth', 1.2)
ylim([-5 5])
xlim([0 ts1(end)])
ylabel('(m)')
title('y: Lateral displacement')
legend('No learning', 'Learning via ADP')
grid on

subplot(413)
plot(ts1, ys1(:,3), ts2, ys2(:,3), 'LineWidth', 1.2)
xlim([0 ts1(end)])
ylim([-5 5]/6)
ylabel('(rad)')
title('\theta: Truck heading')
legend('No learning', 'Learning via ADP')
grid on

subplot(414)
plot(ts1, ys1(:,4), ts2, ys2(:,4), 'LineWidth', 1.2)
xlim([0 ts1(end)])
ylim([-5 5]/6)
grid on
ylabel('(rad)')
xlabel('Time (sec)')
legend('No learning', 'Learning via ADP')
title('\gamma: Relative angle')

%%
set(gcf, 'Color', 'w');
export_fig simResults.png