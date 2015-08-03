function dx = polySystem(x,u)
 % This function discribes the dynamics of a polynomial system
 % dx/dt = F(1)*x + F(2)*x^2 + F(3)*x^3 +  
 F = [0 0.01 0];
 dx = F * (x.^[1 2 3]') + u;
 end