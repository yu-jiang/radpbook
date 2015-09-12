%Fitts_law_main.m


% This file generates the curve fitting using Fitts's Law. It generates
% Figure 8 and Table 2 in the paper.

%NF
disp('Simulating the NF...')
wdd=2.2:0.1:6;
r=0.5./(2.^wdd);
T=r;
for i=1:length(r) % generating each single trial
    disp(['Simulating the NF...', 'Trial: #', num2str(i), ' (', num2str(length(r)-i), 'trials left)'])
    % T(i) = fitts_movement_duration_NF(r(i));
    T(i) = MDF.getPreLearningMovementDuration(r(i));
end
%%
figure(1)
subplot(121)
fit1l=polyfit(wdd,T,1)
plot(wdd,T,'+',wdd,fit1l(2)+fit1l(1).*wdd,'linewidth',1.5)
axis([1,7,0,0.8])
xlabel(['\fontsize{14}log_2(2d/s)'])
ylabel('\fontsize{14}Movement time t_f (sec)')
title(['\fontsize{14}A'])
legend('\fontsize{12}Movement duration of each trial', '\fontsize{12}Least-squares fit using the log law')
subplot(122)
fit1p=polyfit(log(0.25./r),log(T),1)
plot(log(0.25./r),log(T),'+',log(0.25./r),fit1p(2)+fit1p(1).*log(0.25./r),'linewidth',1.5)
legend('\fontsize{12}Movement duration  of each trial', '\fontsize{12}Least-squares fit using the power law')
axis([0,4,-2,0])
xlabel(['\fontsize{14}ln(d/s)'])
ylabel(['\fontsize{14}ln(t_f)'])
title(['\fontsize{14}B'])


%% VF
Tvf = zeros(size(r));
for i=1:length(r)
    disp(['Simulating the VF...', 'Trial: #', num2str(i), ' (', num2str(length(r)-i), 'trials left)'])
    % Tvf(i) = fitts_movement_duration_VF(r(i));
    Tvf(i) = MVF.getPostLearningMovementDuration(r(i));
end
%%
figure(2)
subplot(121)
fit2l=polyfit(wdd,Tvf,1)
plot(wdd,Tvf,'+',wdd,fit2l(2)+fit2l(1).*wdd,'linewidth',1.5)
axis([1,7,0,0.9])
xlabel(['\fontsize{14}log_2(2d/s)'])
ylabel('\fontsize{14}Movement time t_f (sec)')
title(['\fontsize{14}C'])
legend('\fontsize{12}Movement duration of each trial', '\fontsize{12}Least-squares fit using the log law')
subplot(122)
fit2p=polyfit(log(0.25./r),log(Tvf),1)
plot(log(0.25./r),log(Tvf),'+',log(0.25./r),fit2p(2)+fit2p(1).*log(0.25./r),'linewidth',1.5)
legend('\fontsize{12}Movement duration  of each trial', '\fontsize{12}Least-squares fit using the power law')
axis([0,4,-2,0.1])
xlabel(['\fontsize{14}ln(d/s)'])
ylabel(['\fontsize{14}ln(t_f)'])
title(['\fontsize{14}D'])


%% DF
Tdf = zeros(size(r));
for i=1:length(r)
    disp(['Simulating the DF...', 'Trial: #', num2str(i), ' (', num2str(length(r)-i), 'trials left)'])
    Tdf(i) = MDF.getPostLearningMovementDuration(r(i));
end
%%
figure(3)
subplot(121)
fit3l=polyfit(wdd,Tdf,1)
plot(wdd,Tdf,'+',wdd,fit3l(2)+fit3l(1).*wdd,'linewidth',1.5)
axis([1,7,0,0.8])
xlabel(['\fontsize{14}log_2(2d/s)'])
ylabel('\fontsize{14}Movement time t_f (sec)')
title(['\fontsize{14}E'])
legend('\fontsize{12}Movement duration of each trial', '\fontsize{12}Least-squares fit using the log law')
subplot(122)
fit3p=polyfit(log(0.25./r),log(Tdf),1)
plot(log(0.25./r),log(Tdf),'+',log(0.25./r),fit3p(2)+fit3p(1).*log(0.25./r),'linewidth',1.5)
legend('\fontsize{12}Movement duration  of each trial', '\fontsize{12}Least-squares fit using the power law')
axis([0,4,-2,0])
xlabel(['\fontsize{14}ln(d/s)'])
ylabel(['\fontsize{14}ln(t_f)'])
title(['\fontsize{14}F'])


%% Display the fitted parameters

disp('Fitting results:')
disp(['NF Log   Law:', 'a=',num2str(fit1l(1)), ' b=', num2str(fit1l(2))])
disp(['NF Power Law:', 'a=',num2str(fit1p(1)), ' b=', num2str(fit1p(2))])
disp(['VF Log   Law:', 'a=',num2str(fit2l(1)), ' b=', num2str(fit2l(2))])
disp(['VF Power Law:', 'a=',num2str(fit2p(1)), ' b=', num2str(fit2p(2))])
disp(['DF Log   Law:', 'a=',num2str(fit3l(1)), ' b=', num2str(fit3l(2))])
disp(['DF Power Law:', 'a=',num2str(fit3p(1)), ' b=', num2str(fit3p(2))])
