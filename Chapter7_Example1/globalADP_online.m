%Global ADP
clc
clear all
global W F Q
F=-1*[1 0.1 1];
Q=[1 0 0;0 1 0;0 0 1];
W=[0 0 0];
P=[2 0;0 2]*100;

Psave=[];
Qsave=[];

T=0.1;
x=[1;zeros(9,1)]';
for i=1:10
%CXX=[];
C1=[];
C2=[];
C3=[];
CQ=[];
%CZZ=[];
%x=[1;zeros(9,1)]';
for j=0:20
    [t,x] = ode45(@polysys,[j*T,j*T+T],[x(end,1) zeros(1,9)]');
     %[t,x] = ode_yuri(j*T,j*T+T,[x(end,1) zeros(1,9)]',0.001);
    C1 = [C1; 1/2*(x(end,1)^2-x(1,1)^2) 1/3*(x(end,1)^3-x(1,1)^3) 1/4*(x(end,1)^4-x(1,1)^4)];
    C2 = [C2; x(end,2:6)];
    C3 = [C3; x(end,7:9)];
    CQ = [CQ; x(end,10)];
end

cvx_begin sdp

    variable Wn(3,1)
    variable dQs(3,3) symmetric
    %variable Qv(5,1)
    %variable mu


    %Wn=inv((C1-C3)'*(C1-C3))*(-(C1-C3)'*(CQ+C2*Qv))

    Qv=-inv(C2'*C2)*C2'*(CQ+(C1-C3)*Wn(:));

    dQs(1,1)==Qv(1);
    dQs(1,2)+dQs(2,1)==Qv(2);
    dQs(1,3)+dQs(3,1)+dQs(2,2)==Qv(3); 
    dQs(3,2)+dQs(2,3)==Qv(4);
    dQs(3,3)==Qv(5);
    dQs>=0;
    
    
    Pn = [1/2*(Wn(1)) 1/6*(Wn(2));  1/6*(Wn(2)) 1/4*(Wn(3))];
   % mu*eye(2)-P>= 0; %mu*eye(3)-dQs>= 0;
    Pn>=0
    Pn<=P;%-(1e-10)*eye(2);
    %minimize(mu)
    minimize(trace(dQs))

    cvx_end

    Psave=[Psave;P(:)'];
    Qsave=[Qsave;dQs(:)'];
    
    W=Wn;
    P=Pn;
end
Psave
Qsave;
% %


