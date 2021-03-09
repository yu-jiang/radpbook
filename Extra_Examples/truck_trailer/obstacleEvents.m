function [value,isterminal,direction] = obstacleEvents(t,y)
% dDSQdt is the derivative of the equation for current distance. Local
% minimum/maximum occurs when this value is zero.
value = [y(1); y(1)>20] ;  
isterminal = [1; 1];         % stop at local minimum
direction  = [-1; 1];         % [local minimum, local maximum]
end