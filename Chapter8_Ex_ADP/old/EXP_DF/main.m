%main.m

% this file plot figures 5 and 6, takes about one minute

%%
tic
disp('Simulating the NF...')
SimuNF % plot the NF simulation

disp('Simulating the DF...')
SimuDF % plot the learning

disp('Simulating the AL...')%%
SimuAL % plot the after learning

disp('Simulating the AE...')
SimuAE % plot the After Effects
toc