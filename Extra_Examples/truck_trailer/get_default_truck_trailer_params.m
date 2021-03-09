function params = get_default_truck_trailer_params()
params.velocity = -1;
params.truckWheelbase = 3.0;
params.trailerWheelbase = 11.0;
params.feedbackGain = local_get_default_k();
params.forwardGain = local_get_fwd_k();

% for control
params.Q = [1 0 0
     0 1 0
     0 0 1];
params.R = 100;
params.noiseLevel = 100;
params.forwardTarget = [20, 2, 0.3];

% for viz
params.truckWidth         =  3;
params.truckMargin        = 1;
params.trailerWidth       = 2.8;
params.trailerMarginFront = 1;
params.trailerMarginBack  = 3;
params.tireLen   = 1.2;
params.tireWidth = 0.3;
end


%%
function K = local_get_default_k()
% offLinePolicyIteration
D1 = 3.0;     % trailer wheelbase
D2 = 6.0;   % tractor wheelbase
A = [  0    -1     0;
       0     0     0;
       0     0    1/D2];
B = [0;
    -1/D1;
    1/D1];
Q = [1 0 0
     0 1 0
     0 0 1];
R = 1000;
[~, ~, K] = care(A,B,Q,R)
end

%%
function K = local_get_fwd_k()
% offLinePolicyIteration
D1 = 3.0;     % trailer wheelbase
D2 = 6.0;   % tractor wheelbase
A = [  0     1;
       0     0];
B = [0;
     1/D1];
Q = [1 0
     0 1];
R = 1000;
[~, ~, K] = care(A,B,Q,R);
K = [K 0];
end
