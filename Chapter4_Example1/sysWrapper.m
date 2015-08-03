function dx = sysWrapper(t,x,K)
 %  
 x1 = x(1);
 sgm  =  x1.^[1 2 3]';
 e  = (0.01*sin(10*t) + 0.01*sin(3*t) + 0.01*sin(100*t))*noise_on;
 u  = -1/2*K(:)'*sgm + e; 
 dx = polySystem(x,u);
 
 dZ  =  x1.^[2 3 4 5 6]';
 deZ = sgm*e;
 dZZ = kron(sgm,sgm);
 Qk  = Q+1/4*K(:)*K(:)';
 dQ  = sgm'*Qk*sgm; 
 dx  = [dx; 
        dZ;
        deZ;
        dQ;
        ]; % length = 1 + 3*3+3 =13
end