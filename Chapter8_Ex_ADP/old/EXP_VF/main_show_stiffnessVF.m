% This m-file generates the stiffness ellipses in the VF. It requires the
% file "para.m"



theta=0:pi/200:2*pi;
x=cos(theta);
y=sin(theta);
para;
Kxy1=K(1:2,1:2);
Kxy2=Ko(1:2,1:2);
xy1=(Kxy1)*[x;y];
xy2=(Kxy2)*[x;y];
figure
plot(xy1(1,:),xy1(2,:),'g',xy2(1,:),xy2(2,:),'r','Linewidth',2);
axis equal
axis([-1200 1200 -1000 1000]*0.8)
xlabel('x component of stiffness (N m^{-1})')
ylabel('y component of stiffness (N m^{-1})')