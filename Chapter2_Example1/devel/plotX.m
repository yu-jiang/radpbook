% Plot X

t = ScopeData.time;
x = ScopeData.signals.values;

plot(t, x(:,1), '-', ...
	 t, x(:,2), '-.',...
     t, x(:,3), ':', ...
     'LineWidth', 2)
 
legend('x_1', 'x_2', 'x_3')

xlabel('Time (sec)')

annotation(gcf,'textarrow',[0.310714285714285 0.223214285714285],...
	[0.753761904761908 0.719047619047623],'String',{'First Iteration'});

% Create textarrow
annotation(gcf,'textarrow',[0.266071428571429 0.224999999999999],...
	[0.564285714285715 0.51428571428572],'String',{'First Iteration'});

% Create textarrow
annotation(gcf,'textarrow',[0.280357142857143 0.237499999999998],...
	[0.276190476190476 0.326190476190481],'String',{'First Iteration'});

% Create textarrow
annotation(gcf,'textarrow',[0.683928571428571 0.639285714285714],...
	[0.576190476190476 0.473809523809525],'String',{'Last Iteration'});


%%
print('Ch2_ex1_fig2_x', '-depsc')