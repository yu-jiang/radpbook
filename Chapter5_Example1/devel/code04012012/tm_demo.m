%tm_demo

clc
twopara
global KM

G=A1(3,3);
K0=[2.2361    0.2891];
KM=[inv(G)*inv(S)*(N+K0)+M*K0 M]
%KM=[5*36.1484*([-2.9215 -2.9216]+[1.4142 1.4300])+764.8232*[1.4142 1.4300] 764.8232];

[tt1,XX1]=ode45(@tm_sys_demo,[0,1],[0 0 0 0 0 0 0 0 0 0 0 0]');
[tt2,XX2]=ode45(@tm_sys_demo,[1,10],10*[0 0 .1 0 0  -0.1 0 0 .1 0 0 -0.1]');
tt=[tt1;tt2];
XX=[XX1;XX2];

%%
figure(1)
subplot(211)
plot(tt,(XX(:,1)+angle10)*180/pi,tt,(XX(:,7)+angle10)*180/pi,'r-.','Linewidth',1.5)
legend('Robust ADP','Unlearned')
xlabel('time (sec)')
ylabel('Rotor Angle (degree)')
axis([0 10 -100 300])
title('Generator 1')


% Create textbox
annotation(gcf,'textbox',...
    [0.140013062409289 0.348463536897226 0.205095791001451 0.0559105431309892],...
    'String',{'Oscillation started'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.235484760522496 0.617600917088918 0.111394775036285 0.0559105431309892],...
    'String',{'Learning'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(gcf,'textbox',...
    [0.245644412191581 0.151147242967513 0.108492017416546 0.0559105431309892],...
    'String',{'Learning'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);


subplot(212)
plot(tt,(XX(:,4)+angle20)*180/pi,tt,(XX(:,10)+angle20)*180/pi,'r-.','Linewidth',1.5)
legend('Robust ADP','Unlearned')
xlabel('time (sec)')
ylabel('Rotor angle (degree)')
axis([0 10 -10 160])
title('Generator 2')









% Create arrow
annotation(figure1,'arrow',[0.217387518142235 0.207547169811321],...
    [0.363984025559104 0.28594249201278]);

% Create arrow
annotation(figure1,'arrow',[0.188679245283019 0.21044992743106],...
    [0.863217252396166 0.731629392971246]);

% Create textbox
annotation(figure1,'textbox',...
    [0.373046444121916 0.843607306865276 0.240886792452831 0.0559105431309891],...
    'String',{'Controller updated'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.379171262699564 0.339645645523423 0.236213352685051 0.0559105431309892],...
    'String',{'Controller updated'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation(figure1,'arrow',[0.428156748911466 0.380261248185777],...
    [0.851437699680511 0.742811501597444]);

% Create arrow
annotation(figure1,'arrow',[0.400580551523948 0.36742934051144],...
    [0.349840255591054 0.297939778129952]);

% Create textbox
annotation(figure1,'textbox',...
    [0.135978229317852 0.854086540092115 0.205095791001451 0.0559105431309892],...
    'String',{'Oscillation started'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);


figure(2)
subplot(211)
plot(tt,XX(:,2)/2/pi+50,tt,XX(:,8)/2/pi+50,'r-.','Linewidth',1.5)
legend('Robust ADP','Unlearned')
xlabel('time (sec)')
ylabel('Frequency (Hz)')
axis([0 10 45 55])
title('Generator 1')
subplot(212)
plot(tt,XX(:,5)/2/pi+50,tt,XX(:,11)/2/pi+50,'r-.','Linewidth',1.5)
legend('Robust ADP','Unlearned')
xlabel('time (sec)')
ylabel('Frequency (Hz)')
axis([0 10 48 52])
title('Generator 2')
% 
% %%
% 
% figure(3)
% plot(tt,XX(:,1),tt,XX(:,7),'r:')
% xlabel('time (sec)')
% 
% figure(4)
% plot(tt,XX(:,2),tt,XX(:,8),'r:')