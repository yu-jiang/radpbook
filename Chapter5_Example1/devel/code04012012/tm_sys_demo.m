function dxx=tm_sys_demo(t,xx)
twopara
global KM
x=xx(1:3);
e=xx(4:6);
xc=xx(7:9);
ec=xx(10:12);
Pd=B11*Ef1*Ef2*(sin(x(1)+angle10-e(1)-angle20)-sin(angle10-angle20));

if t>=3
    dx=(A1-B1*KM)*x-[0;Pd*w0/2/H1;0];
else
    dx=(A1+B1*K1)*x-[0;Pd*w0/2/H1;0];
end

de=(A2+B2*K2)*e+[0;Pd*w0/2/H2;0];

%control group
Pdc=B11*Ef1*Ef2*(sin(xc(1)+angle10-ec(1)-angle20)-sin(angle10-angle20));
dxc=(A1+B1*K1)*xc-[0;Pdc*w0/2/H1;0];
dec=(A2+B2*K2)*ec+[0;Pdc*w0/2/H2;0];

dxx=[dx;de;dxc;dec];
end