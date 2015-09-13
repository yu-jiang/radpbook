%stiffness_geometry.m


theta=0:pi/100:2*pi;

y = 316.2278 .* sin(theta); 

x = 1000 .* cos(theta) ; 

y1 = 316.2278 .* sin(theta); 

x1 = 223.6068 .* cos(theta) ; 


plot(x,y,'-',x1,y1,'--','linewidth',2)

axis([-1300 1300 -900 900])