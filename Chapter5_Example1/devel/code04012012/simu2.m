function createfigure(X1, YMatrix1, YMatrix2)
%CREATEFIGURE(X1,YMATRIX1,YMATRIX2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 03-Apr-2012 13:52:38

% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(2,1,1,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 10]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot1,[45 55]);
box(subplot1,'on');
hold(subplot1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',subplot1,'LineWidth',1.5);
set(plot1(1),'DisplayName','Robust ADP');
set(plot1(2),'LineStyle','-.','Color',[1 0 0],'DisplayName','Unlearned');

% Create xlabel
xlabel('time (sec)');

% Create ylabel
ylabel('Frequency (Hz)','FontSize',12);

% Create title
title('Generator 1','FontSize',12);

% Create subplot
subplot2 = subplot(2,1,2,'Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 10]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot2,[48 52]);
box(subplot2,'on');
hold(subplot2,'all');

% Create multiple lines using matrix input to plot
plot2 = plot(X1,YMatrix2,'Parent',subplot2,'LineWidth',1.5);
set(plot2(1),'DisplayName','Robust ADP');
set(plot2(2),'LineStyle','-.','Color',[1 0 0],'DisplayName','Unlearned');

% Create xlabel
xlabel('time (sec)');

% Create ylabel
ylabel('Frequency (Hz)','FontSize',12);

% Create title
title('Generator 2','FontSize',12);

% Create legend
legend(subplot1,'show');

% Create legend
legend(subplot2,'show');

% Create textbox
annotation(figure1,'textbox',...
    [0.247788316735487 0.833726379733572 0.0999239130434784 0.0401459854014595],...
    'String',{'Learning'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation(figure1,'arrow',[0.223463687150838 0.209497206703911],...
    [0.670658682634731 0.744510978043911]);

% Create textbox
annotation(figure1,'textbox',...
    [0.369197109545786 0.625097452021982 0.206282608695653 0.0401459854014595],...
    'String',{'Controller updated'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.247597643915472 0.382666849309884 0.0999239130434784 0.0401459854014595],...
    'String',{'Learning'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.137999635657032 0.132353867027026 0.219108695652174 0.0401459854014595],...
    'String',{'Oscillation started'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation(figure1,'arrow',[0.201117318435754 0.201564245810056],...
    [0.177644710578842 0.265429141716565]);

% Create arrow
annotation(figure1,'arrow',[0.399441340782123 0.375698324022347],...
    [0.658682634730539 0.692614770459081]);

% Create arrow
annotation(figure1,'arrow',[0.380335195530726 0.35659217877095],...
    [0.201556886227543 0.235489021956085]);

% Create textbox
annotation(figure1,'textbox',...
    [0.359533519553072 0.15741138086533 0.206282608695653 0.0401459854014595],...
    'String',{'Controller updated'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.139586956521739 0.629980057944158 0.219108695652174 0.0401459854014595],...
    'String',{'Oscillation started'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
