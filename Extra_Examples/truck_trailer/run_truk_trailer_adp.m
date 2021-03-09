%% Setup
p = get_default_truck_trailer_params();
p.noiseLevel = 0.01;
p.trailerWheelbase = 11;
forwardTarget = [20 2 0.3];
xx = [0 -5 0 0];
opt = odeset('Events', @obstacleEvents);
dt = 0.2;

%% setup viz
figure(101)
ax(1) = subplot(121);
tp1 = truck_trailer_plot(p, ax(1));
title('No learning')
ax(2) = subplot(122);
tp2 = truck_trailer_plot(p, ax(2));
title('Learning via ADP')
%% No learning
numPullups = 5;

tsave = [0];
ysave = [xx];

for ct = 1:numPullups

    
%     if ct <=3
%         p.forwardTarget = [20 3 0.1];
%     else
%         p.forwardTarget = [20 0 0];
%     end
%forward
if ysave(end,2) >= -0.3
   p.forwardTarget = [20 0 0];
else
   p.forwardTarget =forwardTarget;
end  

p.velocity = 1;
[t,y] = ode45(@(t, x)sys_truck_trailer_wrapper(t, x, p), ...
    tsave(end)+[0 40], [ysave(end,1:4) zeros(1,10)], opt);

tsave = [tsave; t];
ysave = [ysave; y(:,1:4)];

%backward
p.velocity = -1;
[t,y] = ode45(@(t, x)sys_truck_trailer_wrapper(t, x, p), ...
    tsave(end)+[0 40], [ysave(end,1:4) zeros(1,10)], opt);

tsave = [tsave; t];
ysave = [ysave; y(:,1:4)];

end


ts = (tsave(1):dt:tsave(end))';
[~, ix] = unique(tsave);
ys = interp1(tsave(ix), ysave(ix,:), ts, 'linear', 'extrap');
ts1 = ts;
ys1 = ys;
%%
% for jj = 1:numel(ts)
%     tp1.updateFig(ys(jj,:))    
%     drawnow
% end


%% With ADP


tsave = [0];
ysave = [xx];
psave = zeros(3,3,0);
ksave = zeros(3,1,0);

for ct = 1:numPullups
    
%     if ct <=3
%         p.forwardTarget = [20 3 0.1];
%     else
%         p.forwardTarget = [20 0 0];
%     end
    
if ysave(end,2) >= -0.3
   p.forwardTarget = [20 0 0];
else
   p.forwardTarget = forwardTarget;
end    
%forward
p.velocity = 1;
[t,y] = ode45(@(t, x)sys_truck_trailer_wrapper(t, x, p), ...
    tsave(end)+[0 40], [ysave(end,1:4) zeros(1,10)], opt);
tsave = [tsave; t];
ysave = [ysave; y(:,1:4)];

%backward
p.velocity = -1;
xbase = kron(xx(2:4), xx(2:4));
kbase = [0 0 0];
sbase = [0 0 0 0 0 0];
qbase = 0;

for ii = 1:30
    [t, y] = ode45(@(t, x)sys_truck_trailer_wrapper(t, x, p), ...
        tsave(end)+[0 0.5], [y(end,1:4) zeros(1,10)]);
    
    tsave = [tsave; t];
    ysave = [ysave; y(:,1:4)];
    
    xbase = [xbase;
        kron(y(end, 2:4), y(end, 2:4))];
    qbase = [qbase;
        y(end, 5)];
    kbase = [kbase;
        y(end, 6:8)];
    sbase = [sbase;
        y(end, 9:14)];
end

[t, y] = ode45(@(t, x)sys_truck_trailer_wrapper(t, x, p), ...
    tsave(end)+[0 20], [y(end,1:4) zeros(1,10)], opt);
tsave = [tsave; t];
ysave = [ysave; y(:,1:4)];

    
% Learning
xbase = xbase(2:end,:) - xbase(1:end-1,:);
kbase = kbase(2:end,:);
qbase = qbase(2:end,:);
sbase = sbase(2:end,:);
l = [xbase(:,[1 2 3 5 6 9]) kbase sbase] \ -qbase;
Ve = [l(1) l(2)/2 l(3)/2;
    l(2)/2 l(4) l(5)/2;
    l(3)/2 l(5)/2 l(6)]

Ke = (p.R \ l(7:9)) * 0.5

D1 = 3.0;     % trailer wheelbase
D2 = 11.0;   % tractor wheelbase
A = [  0    -1     0;
    0     0     0;
    0     0    1/p.trailerWheelbase];
B = [0;
    -1/p.truckWheelbase;
    1/p.truckWheelbase];

P = lyap((A-B*p.feedbackGain)',p.Q + p.feedbackGain'*p.R*p.feedbackGain)
Knext = p.R \ B'*P

psave(:,:,end+1) = P;
ksave(:,:,end+1) = p.feedbackGain;
p.feedbackGain = Knext;

end

%
% Resample t
ts = (tsave(1):dt:tsave(end))';
[~, ix] = unique(tsave);
ys = interp1(tsave(ix), ysave(ix,:), ts, 'linear', 'extrap');
ts2 =ts;
ys2 =ys;
% %%
% for jj = 1:numel(ts)
%     tp2.updateFig(ys(jj,:))    
%     drawnow
% end


[Ps, ~, Ks] = care(A,B,p.Q, p.R);
%% Parallel plot
%v = VideoWriter('newfile.mp4','MPEG-4');
%open(v)
for jj = 1:numel(ts1)
    tp1.updateFig(ys1(jj,:))    
    if jj <= numel(ts2)
        tp2.updateFig(ys2(jj,:))
    end    
    drawnow
 %   frame = getframe(gcf);
 %   writeVideo(v,frame);
end
%close(v)

%%
figure(301)
plot(1:5, norm(squeeze(ksave - Ks')))