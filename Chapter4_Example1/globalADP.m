%Global ADP
global W
F=-[1 0.1 1];
Q=[1 0 0;0 1 0;0 0 1];
W=[2 0 0];
P=[2 0;0 2];

Psave=[P(:)'];
Qsave=[];


% for i=1:10
% cvx_begin sdp
% 
% %variable Pn(2,2) symmetric
% variable dQs(3,3) symmetric
% variable Wn(1,3)
% 
% variable mu % Wn(1,3)
% %Wn=[2*Pn(1,1); 6*Pn(1,2); 4*Pn(2,2)]';
% 
% Pn=[1/2*(Wn(1)-W(1)) 1/6*(Wn(2)-W(2));
%    1/6*(Wn(2)-W(2)) 1/4*(Wn(3)-W(3))];
% dQ=-(1/2*(Wn'*(F-1/2*W)+(F-1/2*W)'*Wn)+Q+1/4*(W'*W));%==0;
% dQs(1,1)==dQ(1,1);
% dQs(1,2)+dQs(2,1)==dQ(1,2)+dQ(2,1);
% dQs(1,3)+dQs(3,1)+dQs(2,2)==dQ(1,3)+dQ(3,1)+dQ(2,2);
% dQs(3,2)+dQs(2,3)==dQ(3,2)+dQ(2,3);
% dQs(3,3)==dQ(3,3);
% dQs>=0;
% mu*eye(3)-dQs>= 0;
% Pn>=0;
% minimize(mu)
% 
% cvx_end
% P=[1/2*(Wn(1)) 1/6*(Wn(2));
%    1/6*(Wn(2)) 1/4*(Wn(3))];
% W=Wn;
% 
% Psave=[Psave;P(:)']
% Qsave=[Qsave;dQs(:)']
% end

T=0.2;
x=[.1;zeros(12,1)]';
C1=[];
C2=[];
C3=[];
W=[2 0 0];

W=[2.0084   -0.1074    0.2114]
W=[2.0117   -0.1492    0.2935]
W=[ 2.0146   -0.1864    0.3666]
W=[ 2.0171   -0.2175    0.4277]
W=[  2.0193   -0.2454    0.4825]
W=[2.0212   -0.2697    0.5302]
W=[   2.0229   -0.2910    0.5721]
W=[2.0244   -0.3102    0.6099]
W=[ 2.0257   -0.3265    0.6420]
W=[2.0268   -0.3412    0.6709]
for i=1:10
%pause

x=[10;zeros(12,1)]';
for j=0:10
    [t,x]=ode45(@polysys,[0,T]+j*T,[x(end,1) zeros(1,12)]);
    C1 = [C1; kron([x(end,1) x(end,1)^2],[x(end,1) x(end,1)^2]) - kron([x(1,1) x(1,1)^2],[x(1,1) x(1,1)^2])];
    C2 = [C2; x(end,2:2+8)];
    C3 = [C3; x(end, 2+8+1:2+8+3)];
end


cvx_begin sdp

variable Pn(2,2) symmetric
variable dQs(3,3) symmetric
variable dQ(3,3)
%variable K(3,1)
variable mu % Wn(1,3)
Wn=[2*Pn(1,1); 6*Pn(1,2); 4*Pn(2,2)]';

%Pn=[1/2*(Wn(1)-W(1)) 1/6*(Wn(2)-W(2));
%    1/6*(Wn(2)-W(2)) 1/4*(Wn(3)-W(3))];
%1/2*(Wn'*(F-1/2*W)+(F-1/2*W)'*Wn)+Q+dQ+1/4*(W'*W)==0;
Qk=Q+1/4*W'*W;

C2rd=C2(:,[1,2,3,5,6,9]);     
%C1(:,[1 2 3])*[Pn(1,1) Pn(1,2)*2 Pn(2,2)]' + C2rd*dQv==C3*(1/2*W)'-C2*Qk(:);
C1*Pn(:) + C2*dQ(:)+C2*Qk(:)-C3*W'<=0
C1*Pn(:) + C2*dQ(:)+C2*Qk(:)-C3*W'>=0

dQs(1,1)==dQ(1,1);
dQs(1,2)+dQs(2,1)==dQ(1,2)+dQ(2,1);
dQs(1,3)+dQs(3,1)+dQs(2,2)==dQ(1,3)+dQ(3,1)+dQ(2,2);
dQs(3,2)+dQs(2,3)==dQ(3,2)+dQ(2,3);
dQs(3,3)==dQ(3,3);
dQs>=0;
Pn>=0;
Pn<=P;
mu*eye(3)-dQs>= 0;
minimize(mu)

cvx_end
W=[2*Pn(1,1); 6*Pn(1,2); 4*Pn(2,2)]';
P=Pn;
Psave=[Psave;P(:)'];
Qsave=[Qsave;dQs(:)']
end 
% W=Wn;
% % 


