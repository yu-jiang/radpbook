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
%% setup viz
figure(101)
ax(1) = subplot(131);
tp1 = truck_trailer_plot(p, ax(1));
tp1.updateFig([0 -5 0 0])
title('Initial Pose')
ax(2) = subplot(132);
tp2 = truck_trailer_plot(p, ax(2));
tp2.updateFig([p.forwardTarget 0]);
title('Intermediate Target')
ax(3) = subplot(133);
tp2 = truck_trailer_plot(p, ax(3));
title('Final Target')


%% setup viz
figure(101)
ax(1) = subplot(121);
tp1 = truck_trailer_plot(p, ax(1));
tp1.updateFig([0 -5 0 0])
title('Initial Pose')
ax(2) = subplot(122);
tp2 = truck_trailer_plot(p, ax(2));
title('Desired Target')

set(gcf, 'Color', 'w');
export_fig ttproblem.png



%% setup viz
figure(101)
ax(1) = subplot(121);
tp1 = truck_trailer_plot(p, ax(1));
tp1.updateFig([20 2 0.3 0])
title('Intermediate Target #1')
ax(2) = subplot(122);
tp2 = truck_trailer_plot(p, ax(2));
tp2.updateFig([20 0 0  0])
title('Intermediate Target #2')
%%
set(gcf, 'Color', 'w');
export_fig ttinter.png