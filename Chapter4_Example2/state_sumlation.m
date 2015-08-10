clc
close all
global all
K=Knew;
[t,y]=ode45(@jetsys0,[40 43], [-1 2]);



K=[K0 zeros(2,7)];
[t1,y1]=ode45(@jetsys0,[40 43], [-1 2]);




unew = y;
uold = y1;
for i=1:length(y)
 unew(i,:) = (Knew*sigma(y(i,1), y(i,2)))';
end
for i=1:length(y1)
 uold(i,:) = (K*sigma(y1(i,1), y1(i,2)))';
end
figure(1)
subplot(211)
plot(t,y(:,1),t1,y1(:,1))
legend('new', 'old')


subplot(212)
plot(t,y(:,2),t1,y1(:,2))
legend('new', 'old')

figure(2)
subplot(211)
plot(t,unew(:,1),t1,uold(:,1))
legend('new', 'old')


subplot(212)
plot(t,unew(:,2),t1,uold(:,2))
legend('new', 'old')