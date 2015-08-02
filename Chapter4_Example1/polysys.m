function dx=polysys(t,x)
global W F Q noise_on
 %F  = -[1 0.1 1];
 %W  =  [2 0 0];
 x1 = x(1);
 sgm  =  x1.^[1 2 3]';
 e  = (0.01*sin(10*t)+0.01*sin(3*t)+0.01*sin(100*t))*noise_on;
 u  = -1/2*W(:)'*sgm+e; 
 dx = F * sgm + u;
 
 dZ  =  x1.^[2 3 4 5 6]';
 deZ = sgm*e;
 dZZ = kron(sgm,sgm);
 Qk  = Q+1/4*W(:)*W(:)';
 dQ  = sgm'*Qk*sgm; 
 dx  = [dx; 
        dZ;
        deZ;
        dQ;
        %dZZ
        ]; % 1 + 3*3+3 =13
end