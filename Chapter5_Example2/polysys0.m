function drx=polysys0(t,x)
global W F Q W1 rho
 sgm  =  x(1).^[1 2 3]';
 r = x(end);
 %e  = (0.01*sin(10*t)+0.01*sin(3*t)+0.01*sin(100*t))*noise_on;
 e  = 3*r*x(1)+3*r; 
 u  = -1/2*W1(:)'*sgm-rho*W1(:)'*sgm+e; 
 dx = F * sgm + u;
 dr = - 0.3*r^2-0.3*r*(2*x(1)+x(1)*x(1));
 drx = [dx;dr];
 
 end