function dxx = sys_truck_trailer_wrapper(t, xall, params)
    xx = xall(1:4);
    dxx = sys_truck_trailer(t, xx, params);       %4
    Qk = params.Q + params.feedbackGain'*params.R*params.feedbackGain;
    dqr = xx(2:4)' * Qk * xx(2:4);               %1
    dk  = xx(2:4) * exploration_noise(t, params); %3
    ds  = base_sigma(xx(2), xx(3), xx(4));       %6;
    
    dxx = [dxx; dqr; dk; ds];          % 14
end

