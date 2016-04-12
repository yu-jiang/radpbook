function drx=polysys(t,x)
global W F Q noise_on rho
 %F  = -[1 0.1 1];
 %W  =  [2 0 0];
 x1 = x(1);
 r = x(11);
 sgm  =  x1.^[1 2 3]';
 e  = -rho*W(:)'*sgm+3*r*x(1)+3*r+2*(0.01*sin(10*t)+0.01*sin(1000*t)+0.01*sin(.3*t)+0.01*sin(.1*t))*noise_on;
 u  = -1/2*W(:)'*sgm+e; 
 dx = F * sgm + u;
 
 dZ  =  x1.^[2 3 4 5 6]';
 deZ = sgm*e;
 dZZ = kron(sgm,sgm);
 Qk  = Q+1/4*W(:)*W(:)';
 dQ  = sgm'*Qk*sgm; 
 
 dr = - 0.3*r^2-0.3*r*(2*x(1)+x(1)*x(1));
 
 drx  = [dx; 
        dZ;
        deZ;
        dQ;
        %dZZ
        dr]; % 1 + 3*3+3 =13
end