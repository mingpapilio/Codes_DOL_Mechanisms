%This script creates maps of parameter space (group size and essentiality
%of cooperation) and evaluates which mechanism (coordinated specialisation
%or random specialisaiton) is favoured. We assume here the reproductive
%fecundity is decelerating in the proportion of helpers in the group.

clear all
close all

Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1); %discretisation for group size

alpha=1/(1.5); %shape of return from proportion of helpers (decelerating when alpha<1)

beta=0.01; %shape of coordination cost
lambda=0.1; %scale of coordination costs


essen=linspace(0.1,1,2*length(size)); %discretisation for essentiality of cooperation 


%Iterate through each possible parameter combination and evaluate 
%the fitness of fully coordinated groups and random groups. 
%Store in results. 
results=[]
for i=1:length(essen)
    disp(i)
    for j=1:length(size)
        e=essen(i);
        N=size(j);
        scale=1000;
        offset=0.25;
        Kthin=linspace(0,N,1000);
        K=linspace(0,N,N+1);
        trials=20;
        Wbroad=((N-(K))).*(1-e+e*((K./(N)).^alpha))/N;
        Kopt=K(max(Wbroad)==Wbroad);
        coordfit=(1-lambda.*(1-exp(-beta.*N.*(N-1)/2))).*max(Wbroad);
        stochfit=WK(mean(Kopt)/N,Wbroad);
        results(end+1,:)=[e N mean(Kopt)/N coordfit stochfit];
    end
end


close all


%determine relative fitness
matt=[];
for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        coordfit=results(results(:,2)==N&results(:,1)==e,4);
        stochfit=results(results(:,2)==N&results(:,1)==e,5);
        matt(i,j)=stochfit/coordfit;
    end
end

%generate figure

figure()
ftsz=50;
set(gcf, 'Position',  [0, 0, 800, 800])
y=[];
grid=3;
x=linspace(0.,1,grid)
y=[];
for i=1:grid
    y(end+1,:)=[x(i) x(i) x(i)]
end
y(2,:)=y(1,:)
colormap(y)


[M,c]=contourf(essen, size,1-matt,[0,1,1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;

set(gcf,'color','w');
ax = gca;
ax.XTick=[0 0.5 1];
ax.XTickLabel={'0', '1/2','1'};



ax.YTick=[2 20 40];
ax.FontSize = ftsz/2; 


axis([0 max(essen) min(size) max(size)])


