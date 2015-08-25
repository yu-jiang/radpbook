%plot_movement_duration.m


para;
Tend=zeros(1,30);
L=Tend;

Kadp=K;
[Tend(1),L(1)]=movement_duration(Kadp);

Kadp=K0;
[Tend(2),L(2)]=movement_duration(Kadp);

for s=3:30
    disp(['Simulating the' num2str(s) '-th trial'])
    ADPlearning;
    %Kadp=inv(R)*B'*lyap((A-B*Kadp)',Q1+Kadp'*R*Kadp);
    [Tend(s),L(s)]=movement_duration(Kadp);
end

%%  
TendNF=movement_duration_NF() %calcualte the expected time duration in the NF
    
    
%%
x=1:30;
y=Tend;
%fit = polyfit(x, log(y), 1);
figure
plot(x,y,'-o',x,TendNF*ones(1,30),'r-')
axis([0 32 0.2 1.5])
xlabel('Number of trials')
ylabel('Movement duration (sec)')
legend('Movement duration in VF', 'Expected movement duration in NF')


%%
figure
plot(x,L,'-o', x, .24*ones(1,30),'r-')
axis([0 32 0.2 .4])
xlabel('Number of trials')
ylabel('Movement distance (m)')
legend('Movement distance in VF', 'Expected movement distance in NF')