%This script generate figure S3A. We create a map of parameter space (group size and essentiality
%of cooperation) and evaluate which mechanism (fully coordinated specialisation
%or random specialisation) is favoured. We assume here the reproductive
%fecundity depends linearly on the number of helpers in the group.

clear all
close all

%Parameter space
Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1);  %discretisation for group size

essen=linspace(0.0,1,length(size)); %discretisation for essentiality of cooperation 

beta=0.01; %shape of coordination cost
alpha=0.1; %scale of coordination costs

%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers. Evaluate the fitness of a coordinated
%group, and random groups. Store in results. 
results=[]

for i=1:length(essen)
    disp(i)
    for j=1:length(size)
        e=essen(i);
        N=size(j);
        popt=max((e*(N+1)-1)/(2*e*N),0);
        cn=alpha.*(1-exp(-beta.*N.*(N-1)/2));
        cstar=(N*e-(1-e))/(N*((N-1)*e+1));
        results(end+1,:)=[e N popt cn cstar];
    
    end
end

%Setup for plots
y=[];
grid=3;
x=linspace(0,1,grid)
y=[];
for i=1:grid
    y(end+1,:)=[x(i) x(i) x(i)]
end
y(2,:)=y(1,:)
ftsz=50;
matt2=[];

for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        coordfit=results(results(:,2)==N&results(:,1)==e,4);
        stochcontfit=results(results(:,2)==N&results(:,1)==e,5);
        matt2(i,j)=stochcontfit/coordfit;
    end
end

for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        coordcost=results(results(:,2)==N&results(:,1)==e,4);
        randcost=results(results(:,2)==N&results(:,1)==e,5);
        matt3(i,j)=randcost/coordcost;
        if matt3(i,j)<0
            matt3(i,j)=0 %area where DOL not favoured for either mechanism
        end
    end
end



%Optimal mechanism when the probability of being a helper is equal to the 
%continuous optimal proportion of helpers
figure()
set(gcf, 'Position',  [0, 0, 800, 800])
colormap(flipud(y))
[M,c]=contourf(essen, size,matt3, [0, 1, 1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;
set(gcf,'color','w');
ax = gca;
ax.XTick=[0 0.5 1];
ax.YTick=[2 20 40];
ax.FontSize = ftsz/2; 
axis([min(essen) max(essen) min(size) max(size)])
ax.XTickLabel={'0', '1/2','1'};



endol=results(results(:,3)<=0,1);
nndol=results(results(:,3)<=0,2);
nmax=[];
for i=1:length(essen)
    e=essen(i);
    if sum(endol==e)>0
        nmax(i)=max(nndol(endol==e))
    else
        nmax(i)=0;
    end
end
bottom=zeros(1,length(essen));
hold on
fill([essen flip(essen)],[nmax bottom],'r','LineStyle','none')
