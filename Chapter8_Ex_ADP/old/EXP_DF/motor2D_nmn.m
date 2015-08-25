function dx=motor2D_nmn(t,x)

para;

%x=X(1:6);
u=-K*x;

dx=A*x+B*u;

end