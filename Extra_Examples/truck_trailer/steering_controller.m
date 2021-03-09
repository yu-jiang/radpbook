function [usat, e] = steering_controller(t, x, y, h, g, params)
% steering feedback control
if params.velocity <= 0
    K  = params.feedbackGain;
    u = -K * [y; h ; g];
else
    dh = h - params.forwardTarget(3);
    bx = x - params.forwardTarget(1);
    by = y - params.forwardTarget(2);
    dy = -bx * sin(h) + by * cos(h);
    K  = params.forwardGain;
    u = -K * [dy; dh; g];
end
steeringWheelLimit = tan(deg2rad(25)); % Limtied at 25 deg
e = exploration_noise(t, params);
usat = min(max(u+e, -steeringWheelLimit), steeringWheelLimit);
end