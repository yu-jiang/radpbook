b=10+(rand-0.5)*4;
tau=0.05+(rand-0.5)*0.04;
m=1.3+(rand-0.5)*0.6;
tmp=0;
for i=1:1000
A=[0 1 0;
   0 -b/m 1/m;
   0 0 1/tau;];

B=[0;0;1/tau];

k=[-100 -10 -10];

tmp=tmp+sum(real(eig(A+B*k))<0);

end
tmp
