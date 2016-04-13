function e = ExplorationNoise(t)
%% LocalExploration 
% Generate exploratoion noise
  e = (0.01*sin(10*t)+0.01*sin(3*t)+0.01*sin(100*t));
end