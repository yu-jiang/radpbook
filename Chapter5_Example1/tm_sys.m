function dxx=tm_sys(t,xx)
global KM
xx=xx(:);
twopara
K0=[2.2361    0.2891];
x=xx(1:3);
e=xx(4:6);
B11=0.05;
Pd=B11*Ef1*Ef2*(sin(x(1)+angle10-e(1)-angle20)-sin(angle10-angle20));

if t>=3
    u=-KM*x;
else
    u=K1*x+.1*sin(100*t);
end


%for checking purpose



% Pd=0;
dx=A1*x+B1*u-[0;Pd*w0/2/H1;0];

A=A1(1:2,1:2);
B=A1(1:2,3);
F=A1(3,3);
G=B1(3);

xi=xx(20);dxi=(K0*(A-B*K0)-F*K0)*x(1:2)+(F+K0*B)*xi+G*u;

%dx(3)=-(K0*(A-B*K0)-F*K0)*x(1:2)+F*x(3)+G*u;

%xi=xx(20);%dxi=%F*x(3)+K0*dx(1:2)+G*u;

%
%    %for checking purpose
%    dx(1)=0;
%    dx(2)=0;
%    %for checking purpose
%
de=(A2+B2*K2)*e+[0;Pd*w0/2/H2;0];

Ixx=xx(7:10);    Ixu=xx(11:12);
Ixz=xx(13:14);   Izz=xx(15);
Izu=xx(16);      Idx=xx(17:18);
Idz=xx(19);


dIxx=kron(x(1:2),x(1:2));
dIxu=x(1:2)*(u);
dIzz=x(3)^2;
dIxz=x(1:2)*x(3);%*(x(3)-Pd);
dIzu=x(3)*(u);
dIdx=x(1:2)*(-Pd);
dIdz=x(3)*(-Pd);


Ixiu=xx(21);
Ixixi=xx(22);
Ixix=xx(23:24);



dIxiu=xi*u;
dIxixi=xi*xi;
dIxix=xi*x(1:2);


xc=xx(25:27);
ec=xx(28:30);

Pdc=B11*Ef1*Ef2*(sin(xc(1)+angle10-ec(1)-angle20)-sin(angle10-angle20));
dxc=(A1+B1*K1)*xc-[0;Pdc*w0/2/H1;0];
dec=(A2+B2*K2)*ec+[0;Pdc*w0/2/H2;0];


dxx=[dx;de;dIxx;dIxu;dIxz;dIzz;dIzu;dIdx;dIdz;dxi;dIxiu;dIxixi;dIxix;dxc;dec];

end