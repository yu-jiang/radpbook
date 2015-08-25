function Tend=movement_duration_NF()
para;
Kadp=K;
[t,y]=ode_yu_md(0,2,[0,-.25,0,0,0,0]',dt,Kadp)
%plot(y1(:,1),y1(:,2))

i=length(t);
while norm(y(i,1:2))<0.01
 i=i-1;
end
Tend=t(i);

end

function [t,Y]=ode_yu_md(t0,tf,x0,dt,Kadp)
h=dt;                                                         
t=t0:dt:tf;                                                                            
y=x0;       
Y=[];
    for clock=t  
        clock;
       y=y+motor2D(clock,y,Kadp)*h;      
       Y=[Y y];
    end 
Y=Y';
end

function dX=motor2D(t,X,Kadp)
para;
%global Kadp
x=X(1:6);
u=-Kadp*x;
w=randn(2,1)*sqrt(dt);
%M=[c1*u(1) c2*u(2); -c2*u(1) c1*u(2)];
%v=M*w; % control dependent noise
dx=A0*x+B*u;%+B*v./dt; %6
%dIq=x'*(Q+K'*R*K)*x; %1
%dxxi=kron(x',v'*R)'./dt; %12;
%duu=kron(v,v)./dt;         % 4
dX=dx;
end