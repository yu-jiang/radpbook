%Global ADP
% global W F Q
F = [0 0.01 0];

Q=[1 0 0;0 1 0;0 0 1];
Q =0.01*[1 0 0;0 1 0;0 0 0];
W=[0.1 0 0.01];
P=[10 0;0 10]; 
Psave=[];
Qsave=[];

Qk = Q + 1/4*W'*W;

for i=1:3
    cvx_begin sdp
    
    variable dQs(3,3) symmetric
    variable Wn(1,3)
    variable mu
    
    % Wn(1,3)
    %Wn=[2*Pn(1,1); 6*Pn(1,2); 4*Pn(2,2)]';
    %Pn=[1/2*(Wn(1)-W(1)) 1/6*(Wn(2)-W(2));1/6*(Wn(2)-W(2)) 1/4*(Wn(3)-W(3))];
    %dQ=-(1/2*(Wn'*(F-1/2*W)+(F-1/2*W)'*Wn)+Q+1/4*(W'*W));%==0;
    
    dQ=-(Wn'*(F-1/2*W)+Q+1/4*(W'*W));
    dQs(1,1)==dQ(1,1);
    dQs(1,2)+dQs(2,1)==dQ(1,2)+dQ(2,1);
    dQs(1,3)+dQs(3,1)+dQs(2,2)==dQ(1,3)+dQ(3,1)+dQ(2,2);
    dQs(3,2)+dQs(2,3)==dQ(3,2)+dQ(2,3); 
    dQs(3,3)==dQ(3,3);
    dQs>=0;
    Pn=[1/2*(Wn(1)) 1/6*(Wn(2)); 1/6*(Wn(2)) 1/4*(Wn(3))];
    Pn>=0;
    %P<=(1e-8)*eye(2);
    %mu*eye(3)-dQs>= 0;
    %mu*eye(2)-P>= 0;
    Pn<=P-(1e-10)*eye(2);
    minimize(trace(dQs))
    
    cvx_end;
    W=Wn;
    P=Pn;
    Psave=[Psave;P(:)']
    Qsave=[Qsave;dQs(:)']
end  





