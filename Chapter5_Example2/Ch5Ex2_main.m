%% Global Robust Adaptive Dynamic Programming
global W W1 F Q noise_on rho
rho = 0.5 ; % Robust Redesign Gain

F = [0 -3/2 -1/2];
Q = [5 0 0;0 0 0;0 0 0];
W = [2 -1.4 -.45];
wsave = W;
W1 = W;
P = eye(2)*10;
Pold = -100*eye(2);
noise_on=1;

Psave=[];
Qsave=[];
ua=[];

Trjsave=[];
Tsave=[];

T=0.005;
xinit=-10;
rinit=5;
x=[xinit;zeros(9,1);rinit]';
for i=0:9
    %CXX=[];
    C1=[];
    C2=[];
    C3=[];
    CQ=[];
    %CZZ=[];
    %x=[1;zeros(9,1)]';
    for j=0:49
        [t,x] = ode45(@polysys,[j*T,j*T+T]+50*i*T,[x(end,1) zeros(1,9) x(end,11)]');
        %[t,x] = ode_yuri(j*T,j*T+T,[x(end,1) zeros(1,9)]',0.001);
        C1 = [C1; 1/2*(x(end,1)^2-x(1,1)^2) 1/3*(x(end,1)^3-x(1,1)^3) 1/4*(x(end,1)^4-x(1,1)^4)];
        C2 = [C2; x(end,2:6)];
        C3 = [C3; x(end,7:9)];
        CQ = [CQ; x(end,10)];
        Tsave=[Tsave;t];
        Trjsave=[Trjsave;x(:,[1 11])];
        for k=1:length(t)
            ua=[ua; W(:)'*x(k,1).^[1 2 3]'+(0.01*sin(10*t(k))+0.01*sin(3*t(k))+0.01*sin(100*t(k)))*noise_on];
        end
    end
    
    if norm(P(:)-Pold(:))>0.1
        cvx_begin sdp       
        variable Wn(3,1)
        variable dQs(3,3) symmetric
        Qv=-([C2 -C3]'*[C2 -C3])\[C2 -C3]'*(CQ+C1*Wn(:));
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
[t0,y0] = ode45(@polysys0,[0 Tsave(end)],[xinit, rinit]);
for i=1:length(t0)
 u0(i) = W1*y0(i,1).^[1 2 3]'; 
end
plot(Tsave,Trjsave(:,1), 'b-', t0,y0(:,1), 'r-.', 'linewidth', 2)
legend('With GRADP-based controller', 'With initial controller')
xlabel('time (sec)')
ylabel('\phi')
ylim([-12 5])
% Create textarrow
annotation(figure(1),'textarrow',[0.267857142857143 0.228200808625336],...
	[0.245238095238095 0.433832086450543],'String',{'1st iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(1),'textarrow',[0.351785714285714 0.291071428571429],...
	[0.511904761904762 0.666666666666667],'String',{'2nd iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(1),'textarrow',[0.38409703504043 0.386738544474392],...
	[0.763163519772001 0.700619878874244],'String',{'3rd iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(1),'textarrow',[0.603571428571428 0.535714285714285],...
	[0.533333333333334 0.673809523809524],'String',{'5th (final) iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(1),'textarrow',[0.480357142857143 0.427966101694915],...
	[0.528571428571429 0.669064748201439],'String',{'4th iteration'},...
	'FontSize',12);



figure(2)
plot(Tsave,Trjsave(:,2), 'b-', t0,y0(:,2), 'r-.', 'linewidth', 2)
legend('With GRADP-based controller', 'With initial controller')
xlabel('time (sec)')
ylabel('r')
ylim([-1 5])
% Create textarrow
annotation(figure(2),'textarrow',[0.243684992570579 0.1996691805209],...
	[0.676982591876209 0.278255039384132],'String',{'1st iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(2),'textarrow',[0.325408618127786 0.287523619950097],...
	[0.560928433268859 0.281012459180867],'String',{'2nd iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(2),'textarrow',[0.407132243684993 0.373938714289719],...
	[0.502901353965184 0.289622365749069],'String',{'3rd iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(2),'textarrow',[0.473997028231798 0.450802097059071],...
	[0.433268858800774 0.291719058253785],'String',{'4th iteration'},...
	'FontSize',12);

% Create textarrow
annotation(figure(2),'textarrow',[0.592867756315007 0.537890044576523],...
	[0.586073500967118 0.286266924564797],'String',{'5th (final) iteration'},...
	'FontSize',12);



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

vsxt=0;
for i=1:length(x)-1
y=x(i);
vn=[vn [y y^2]*P*[y ;y^2]];
v1=[v1 [y y^2]*P1*[y ;y^2]];
%vsx=(y^2 + 2)^(3/2)/15 - (2*2^(1/2))/15 - y^2/10;

vsx = 2.0*y*(0.25*y^4 + 1.5*y^3 + 2.25*y^2 + 5.0)^(1/2) - 3.0*y^2 - 1.0*y^3;
vsxt = vsxt+vsx*(x(i+1)-x(i));
% vsx=y^3/150 + (101*y^2 + 100)^(3/2)/15150 - 20/303;
vs=[vs vsxt];
end

figure(3)
plot(x(2:end),v1,'g:',x(2:end),vn,'r-.',x(2:end),vs+21.1238,'b','linewidth',2)
legend('V_1 : Initial cost', 'V^5: Improved cost', 'V^o: Optimal cost')
xlabel('\phi')

% figure(4)
% plot(x,u1,'g:',x,un,'r-.',x,us,'b','linewidth',2)
% legend('u_1: Initial control policy', 'u^5: Improved control policy', 'u^o: Optimal control policy')