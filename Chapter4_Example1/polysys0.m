function dx = polysys0(t,x)
global W F Q W1
 sgm  =  x.^[1 2 3]';
 u  = -1/2*W1(:)'*sgm;%+e; 
 dx = F * sgm + u;
 end