function status = createAnimation(t, y, ~)
persistent pmgr
if isempty(pmgr)
	pmgr = paramMgr.getInstance();
end
status = 0;
if ~isempty(t)
	t = t(end);
	delta = y(1, end);
	delta1 = y(4, end);
	delta_ul = y(25, end);
	delta1_ul = y(28, end);
	hf = figure(1);
	ax1 = subplot(211);
	LocalCreateAnimation(t,pmgr.angle10,delta + pmgr.angle10, ...
		pmgr.angle20,delta1 + pmgr.angle20,true,1*(t>2)+1*(t>3), pmgr, hf, ax1, true, 'ADP Learning')
	
	ax2 = subplot(212);
	LocalCreateAnimation(t,pmgr.angle10,delta_ul + pmgr.angle10, ...
	  pmgr.angle20,delta1_ul + pmgr.angle20,true,0, pmgr, hf, ax2, false, 'Unlearned')
end
end


function LocalCreateAnimation(t,delta,delta0,delta1,delta10,ifct,ifadp, pmgr, hf, ax, ifshowtime, caller)
% utCreateAnimation creates the figure of two sync machines. 

B12 = pmgr.B12;
B11 = pmgr.B11;
Ef1 = pmgr.Ef1;
Ef2 = pmgr.Ef2;

% figure(1)
set(hf,'color',[255,255,255]/255)
%delta=1;
r1=0.15;
r2=0.1;
theta=0:2*pi/1000:2*pi;
xx1=r1*cos(theta);
yy1=r1*sin(theta);
xx2=r2*cos(theta);
yy2=r2*sin(theta);


x0=.4;
y0=.4;
x1=.85;
y1=.4;

plot(ax,xx1+x0,yy1+y0,'linewidth',1,'color',[1,0,0])
hold on

plot(ax,xx2+x1,yy2+y1,'linewidth',1,'color',[1,0,0])



line([x0,x0+r1*cos(delta)],[y0,y0+r1*sin(delta)],'linewidth',2,'color',[0,.41,0])
line([x0,x0+r1*cos(delta0)],[y0,y0+r1*sin(delta0)],'linewidth',3,'color',[0,0,255]/255)



%---------------------------------------------------
%----------- Draw the infinite Bus -----------------
line([.05,x0-r1],[y0,y0],'color',[0,0,0],'linewidth',2)
line([.05+0.05,x0-r1-0.05],[y0,y0],'color',[0,0,0],'linewidth',4)
line([.05,.05],[.2,.6],'color',[0,0,0],'linewidth',3)
line([0,.05],[.2,.25],'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05+0.05+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05+0.05+0.05+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05+0.05+0.05+0.05+0.05,'color',[0,0,0],'linewidth',.5)
line([0,.05],[.2,.25]+0.05+0.05+0.05+0.05+0.05+0.05+0.05,'color',[0,0,0],'linewidth',.5)
%===================================================

%---------------------------------------------------
% show time
if ifshowtime
	title(['t=',num2str(t,'%3.2f')],'fontsize',20)
end
%===================================================

ylabel(caller, 'fontsize',16)

%---------------------------------------------------
%   Draw the connection
lineColor = [0, 0, 0];
if ifct==1
	patch([.73 .73 .87 .87]+0.05,[0.2 .25 .25 .2],[0,0,0])
	% text(0.6,0.45,'Connected','fontsize',10,'color',[0,0.8,0])
	%draw the connection
	line([.13,0.13]+0.1,[y0,y0+0.16],'color',lineColor,'linewidth',2)
	line([.13,0.47]+0.1,[y0+0.16,y0+0.16],'color',lineColor,'linewidth',2)
	line([.47,0.47]+0.1,[y0+0.16,y0],'color',lineColor,'linewidth',2)
	line([.57,x1-r2],[y0,y0],'color',lineColor,'linewidth',2)
	line([.62,0.72],[y0,y0],'color',lineColor,'linewidth',4)
else
	patch([.73 .73 .87 .87]+0.05,[0.2 .25 .25 .2],[0,0,0])
	text(0.6,0.45,'Disonnected','fontsize',10,'color',[0,0,0])
	%draw the connection
	line([.13,0.13]+0.1,[y0,y0+0.16],'color',lineColor,'linewidth',2,'linestyle',':','color','c')
	line([.13,0.47]+0.1,[y0+0.16,y0+0.16],'color',lineColor,'linewidth',2,'linestyle',':','color','c')
	line([.47,0.47]+0.1,[y0+0.16,y0],'color',lineColor,'linewidth',2,'linestyle',':','color','c')
	line([.57,x1-r2],[y0,y0],'color',lineColor,'linewidth',2,'linestyle',':','color','c')
	%line([.62,0.72],[y0,y0],'color',lineColor,'linewidth',4,'linestyle',':')
end
%===================================================


if ifadp==0
	patch([.3 .3 .5 .5],[0.12 .2 .2 .12],[0,0,0])
end
if ifadp==1
	patch([.3 .3 .5 .5],[0.12 .2 .2 .12],[1,0,0])
end
if ifadp==2
	patch([.3 .3 .5 .5],[0.12 .2 .2 .12],[0,1,0])
end

line([0.3,0.3]+0.1,[.2,.25],'color',[0,0,0],'linewidth',4)
line([0.8,0.8]+0.05,[.25,.3],'color',[0,0,0],'linewidth',4)


Pe20 = 0.1;


if ifct==1
	Pe1=B11*Ef1*sin(delta0)+B12*Ef1*Ef2*sin(delta0-delta10);
	Pe2=-B12*Ef1*Ef2*sin(delta0-delta10)+Pe20;
	str_sta=['\Delta \delta_1=',num2str(delta0-delta,'%2.4f'),' rad'];
	str_p=['Pe_1=', num2str(Pe1,'%2.4f'), ' p.u.'];
	text(0.30,0.65,[str_sta])
	text(0.31,0.60,[str_p])
	
	line([x1,x1+r2*cos(delta1)],[y1,y1+r2*sin(delta1)],'linewidth',2,'color',[0,.41,0])
	line([x1,x1+r2*cos(delta10)],[y1,y1+r2*sin(delta10)],'linewidth',3,'color',[0,0,255]/255)
	
	str_sta=['\Delta \delta_2=',num2str(delta10-delta1,'%2.4f'),' rad'];
	text(0.70,0.60,str_sta)
	str_p=['Pe_2=', num2str(Pe2,'%2.4f'), ' p.u.'];
	text(0.71,0.55,[str_p])
else
	Pe1 = B11*Ef1*sin(delta0);
	Pe2 = Pe20;
	str_sta = ['\Delta \delta_1=',num2str(delta0-delta,'%2.4f'),' rad'];
	str_p = ['Pe_1=', num2str(Pe1,'%2.4f'), ' p.u.'];
	text(0.30,0.65,[str_sta])
	text(0.31,0.60,[str_p])
	line([x1,x1+r2*cos(delta10)],[y1,y1+r2*sin(delta10)],'linewidth',3,'color',[0,0,255]/255)
	
	str_p = ['Pe_2=', num2str(Pe2,'%2.4f'), ' p.u.'];
	text(0.71,0.55,[str_p])
end


if ifadp==0
	text(0.3,0.05,'Initial Control Law','fontsize',12,'color',[0,0,0])
end
if ifadp==1
	text(0.3,0.05,'Online Learning','fontsize',12,'color',[.5,0,0])
end
if ifadp==2
	text(0.3,0.05,'Post Learning','fontsize',12,'color',[0,.5,0])
end
% polylogo=imread('c:\logo.jpg');
% imagesc([0.72 1],[0.08 0],polylogo);
axis([0,1,0,1*.8])
axis equal 
axis off
hold off
drawnow
end