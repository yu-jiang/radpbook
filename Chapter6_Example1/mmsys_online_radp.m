function dxx=mmsys_online_radp(t,xx)
global Kadp
mm_para %load the parameters

x=[0;0;0];
for i=2:Nm
    id=(i-2)*3+1:(i-2)*3+3;
    x(:,i)=xx(id);
end


K=zeros(1,3,Nm-1);
for i=1:Nm-1
    K(:,:,i)=[10   50   0];
end

% calculate angular differences in matrix form
dlt=zeros(Nm);d=zeros(1,Nm);
for i=1:Nm
    dlt(i,:)=x(1,i)+dlt0(i)-x(1,:)-dlt0';
    d(i)= E(i)*(E.*(x(2,i)-x(2,:)))*(BX(i,:).*cos(dlt(i,:))-GX(i,:).*sin(dlt(i,:)))';
end


u=0;
for i=2:Nm
    if t>=4 & t<=5 % learning stage is between 4s and 5s
        u(i)=-K(:,:,i-1)*x(:,i)+0.001*sin(100*t);
    else
        u(i)=-K(:,:,i-1)*x(:,i);
    end
end

if t>5 % update the control policies and start the post-learning stage
    for i=2:Nm
        u(i)=-Kadp(:,:,i-1)*x(:,i);
    end
end

for i=2:Nm
    dx(:,i-1)=A(:,:,i)*x(:,i)+B(:,:,i)*(u(i)-d(i));
end

dIxxu=zeros(12,Nm-1);
if t>=4 & t<=5
    for i=2:Nm
     dIxxu(:,i-1)=[kron(x(:,i),x(:,i));kron(x(:,i),u(i)-d(i))];
    end
end
dxx=[dx(:);dIxxu(:)];
end