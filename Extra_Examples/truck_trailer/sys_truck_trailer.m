function dxx = sys_truck_trailer(t, xx, params)
   x = xx(1);
   y = xx(2);
   h = xx(3);
   g = xx(4);
   
   v  = params.velocity;
   d1 = params.truckWheelbase;
   d2 = params.trailerWheelbase; 
%  K  = params.feedbackGain;
   
%    % steering feedback control
%    steeringWheelLimit = tan(deg2rad(25)); % Limtied at 25 deg
%    u = -K * [y; h; g];
%    e = exploration_noise(t, params);
%    usat = min(max(u+e, -steeringWheelLimit), steeringWheelLimit);
   [usat, ~] = steering_controller(t, x, y, h, g, params);
   
   dx = v * cos(h);
   dy = v * sin(h);
   dh = v/d1 * usat;
   dg = -v/d2 * sin(g) -  v/d1 * usat;
   
   dxx = [dx; dy; dh; dg];
end