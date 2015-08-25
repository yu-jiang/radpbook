function dx=motor2D_nmn(t,x)
global Kadp
para;

%x=X(1:6);
u=-Kadp*x;

dx=A*x+B*u;

end