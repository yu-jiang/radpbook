function dx=polysys0(t,x)
global W F Q W1
 sgm  =  x.^[1 2 3]';
 %e  = (0.01*sin(10*t)+0.01*sin(3*t)+0.01*sin(100*t))*noise_on;
 u  = -1/2*W1(:)'*sgm;%+e; 
 dx = F * sgm + u;
 end