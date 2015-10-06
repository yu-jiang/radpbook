% checker

global K0
K0=[2.2361    0.2891];

[t,X]=ode45(@tm_sys,[0,10],[0.5 0 0 -0.5 0 0 zeros(1,13) 2.2361*0.5 zeros(1,4)]');
%tt=[tt;t];
%XX=[XX;X];
Dxx=[kron(X(end,1:2),X(end,1:2))-kron(X(1,1:2),X(1,1:2))];
Dzz=[X(end,3)^2-X(1,3)^2];
Dxz=[X(end,1:2)*X(end,3)-X(1,1:2)*X(1,3)];
Ixx=[X(end,7:10)-X(1,7:10)];
Ixu=[X(end,11:12)-X(1,11:12)];
Ixz=[X(end,13:14)-X(1,13:14)];
Izz=[X(end,15)-X(1,15)];
Izu=[X(end,16)-X(1,16)];
Idx=[X(end,17:18)-X(1,17:18)];
Idz=[X(end,19)-X(1,19)];
%%
n=40
XX(:,1:2)*K0'+XX(:,3)-XX(:,20)