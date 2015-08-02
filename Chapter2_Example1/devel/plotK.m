% Plot P

a = getElement(logsout, 'Kdata');
t = squeeze(a.Values.time);
K = squeeze(a.Values.data);

plot(t, K(1,:), '-', ...
	 t, K(2,:), '-.',...
     t, K(3,:), ':', ...
     'LineWidth', 2)
 
legend('K_1', 'K_2', 'K_3')

xlabel('Time (sec)')

%%
print('Ch2_ex1_fig3_K', '-depsc')