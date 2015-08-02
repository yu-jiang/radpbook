function y = isInRegionOfAttraction(x,p,D)
% Check to see if x is in the Region or Attraction 
y = (p'*Phi_fun(x)' <= D);
end
