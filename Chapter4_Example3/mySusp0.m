function dx = mySusp0(x,u,r)
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

[x1,x2,x3,x4] = deal(x(1),x(2),x(3),x(4));

% State matrices
A = [ 0 1 0 0;
    [-ks -bs ks bs]/mb ; ...
    0 0 0 1;
    [ks bs -ks-kt -bs]/mw];
B = [ 0; 10000/mb ; 0 ; -10000/mw];
B1 = [ 0; 0 ; 0 ; kt/mw];

if ~isdouble(u)
    u = eval(u);
end
    
if t <= 0.1
    r = 0.1;
else
    r = 0;
end
dx = A*x + B*u + B1*r;
end