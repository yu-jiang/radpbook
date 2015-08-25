%main.m

% this file plot figures 1 and 2, takes about one minute

%%
tic
disp('Simulating the NF...')
SimuNF % plot the NF simulation

disp('Simulating the VF...')
SimuVF % plot the learning

disp('Simulating the AL...')%%
SimuAL % plot the after learning

disp('Simulating the AE...')
SimuAE % plot the After Effects
toc