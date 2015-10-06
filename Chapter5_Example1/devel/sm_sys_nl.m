function dxx=sm_sys_nl(t,xx)
global A1 A2 B1 K1 K2 B22 B2 B12 B11 B21 Ef1 Ef2 angle10 angle20 ts
   x=xx(1:3);
   z=xx(4:6);
   if t>=ts
    dx=((A1-B1*K1)*x+[0;0;1]*(B11*Ef1*cos(x(1)-angle10)*x(2)+B12*Ef1*Ef2*cos(x(1)-z(1)-angle10+angle20)*(x(2)-z(2))));
    dz=((A2-B2*K2)*z+[0;0;1]*(B22*Ef2*cos(z(1)-angle20)*z(2)+B12*Ef1*Ef2*cos(z(1)-x(1)+angle10-angle20)*(z(2)-x(2))));
   else
    dx=(A1-B1*K1)*x;
    dz=zeros(3,1);
   end
    dxx=[dx;dz];  
end