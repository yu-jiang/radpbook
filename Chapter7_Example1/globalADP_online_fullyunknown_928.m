%Global ADP
clc
clear all
close all
global W W1 F Q noise_on
%F=[-.1 0.05 -.1];
F=[0 0.01 0];
Q=.01*[1 0 0;0 1 0;0 0 0];
W=[.1 0 .01];
wsave=W;
W1=W;
P=[2 0;0 2]*100;
P=eye(2)*10;
Pold=-100*eye(2);
noise_on=1;

Psave=[];
Qsave=[];
ua=[];

Trjsave=[];
Tsave=[];

T=0.1;
x=[2;zeros(9,1)]';
for i=0:9
    %CXX=[];
    C1=[];
    C2=[];
    C3=[];
    CQ=[];
    %CZZ=[];
    %x=[1;zeros(9,1)]';
    for j=0:49
        [t,x] = ode45(@polysys,[j*T,j*T+T]+50*i*T,[x(end,1) zeros(1,9)]');
        %[t,x] = ode_yuri(j*T,j*T+T,[x(end,1) zeros(1,9)]',0.001);
        C1 = [C1; 1/2*(x(end,1)^2-x(1,1)^2) 1/3*(x(end,1)^3-x(1,1)^3) 1/4*(x(end,1)^4-x(1,1)^4)];
        C2 = [C2; x(end,2:6)];
        C3 = [C3; x(end,7:9)];
        CQ = [CQ; x(end,10)];
        Tsave=[Tsave;t];
        Trjsave=[Trjsave;x(:,1)];
        for k=1:length(t)
            ua=[ua; W(:)'*x(k,1).^[1 2 3]'+(0.01*sin(10*t(k))+0.01*sin(3*t(k))+0.01*sin(100*t(k)))*noise_on];
        end
    end
    
    if norm(P(:)-Pold(:))>0.001
        cvx_begin sdp       
        variable Wn(3,1)
        variable dQs(3,3) symmetric
        Qv=-inv([C2 -C3]'*[C2 -C3])*[C2 -C3]'*(CQ+C1*Wn(:));
        dQs(1,1)==Qv(1);
        dQs(1,2)+dQs(2,1)==Qv(2);
        dQs(1,3)+dQs(3,1)+dQs(2,2)==Qv(3);
        dQs(3,2)+dQs(2,3)==Qv(4);
        dQs(3,3)==Qv(5);
        dQs>=0;  
        Pn = [1/2*(Wn(1)) 1/6*(Wn(2));  1/6*(Wn(2)) 1/4*(Wn(3))];
        Pn<=P;
        
        %minimize([-1 100]*Pn*[-1 ;100]+[1 100]*Pn*[1 ;100])
        minimize(Pn(1,1)+Pn(2,2))
        W=Qv(6:8);
        wsave=[wsave;W(:)' ];
        cvx_end
        noise_on=1;
        Psave=[Psave;P(:)'];
        Qsave=[Qsave;dQs(:)'];
        Pold=P;
        P=Pn;
    else
        noise_on=0;
        disp(num2str(i))
    end
    

    
end
Psave
Qsave;
%%
figure(1)
[t0,y0]=ode45(@polysys0,[0 50],2)


for i=1:length(t0)
 u0(i)=W1*y0(i).^[1 2 3]';
end


plot(Tsave,Trjsave, 'b-', t0,y0, 'r-.', 'linewidth', 2)
axis([0 50 -.5 2])
legend('With GADP-based controller', 'With initial controller')
xlabel('time (sec)')
ylabel('x')
% %
syms v(y) 
%vsx=dsolve(diff(v)*(F(1)*y+F(2)*y^2+F(3)*y^3)+0.01*(y^2+y^4)-1/4*(diff(v))^2==0, v(0)==0)

%x=-1.5:0.05:1.5;
x=-2.5:.01:2.5;
vn=[];
v1=[];    
vs=[];
us=[];
u1=[];
un=[];
P1=[Psave(2,1) Psave(2,2);Psave(2,3) Psave(2,4)] ;

for y=x 

vn=[vn [y y^2]*P*[y ;y^2]];
v1=[v1 [y y^2]*P1*[y ;y^2]];
%vsx=(y^2 + 2)^(3/2)/15 - (2*2^(1/2))/15 - y^2/10;
 vsx=y^3/150 + (101*y^2 + 100)^(3/2)/15150 - 20/303;
vs=[vs vsx];
u1=[u1 -1/2*W1*[y;y^2;y^3]];
un=[un -1/2*W'*[y;y^2;y^3]];
us=[us -1/2*((y*(y*(101*y^2 + 100)^(1/2) + 101*y^2 + 100))/(50*(101*y^2 + 100)^(1/2)))];
end
state_anno


figure(2)
plot(Tsave,ua,'Linewidth',2)
legend('u')
xlabel('time (sec)')
axis([0 50 -0.5 2])
control_anno

figure(3)
plot(x,v1,'g:',x,vn,'r-.',x,vs,'b','linewidth',2)
legend('V_1 : Initial cost', 'V^5: Improved cost', 'V^o: Optimal cost')

figure(4)
plot(x,u1,'g:',x,un,'r-.',x,us,'b','linewidth',2)
legend('u_1: Initial control policy', 'u^5: Improved control policy', 'u^o: Optimal control policy')