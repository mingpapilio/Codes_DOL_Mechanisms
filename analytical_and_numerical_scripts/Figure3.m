%This script generates figure 3 of the main text. We generate a
%discretisation of parameter space (group size and essentiality) and for
%each combination of parameter values, we evaluate whether random
%specialisation or fully coordinated specialisation has higher fitness. We
%do this for two different 'scales' of coordination costs.
clear all
close all

%parameter values and space discretisation

beta=0.01; %shape of costs
lambda=0.1; %scale of costs
Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,100); %discretisation for group size
essen=linspace(0.1,1,2*length(size)); %discretisation for essentiality of cooperation

matt=[]; %store relative fitness of mechanisms here
%iterate through parameter space
for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        stochcost=max((2.*e-1)./N,0); %analytical solution
        coordcost=(lambda.*(1-exp(-beta*N.*(N-1)/2))); %cost of coordination function
        matt(i,j)=stochcost/coordcost; %relative fitness of mechanisms
    end
end

%Generate figure
figure()
ftsz=50;
set(gcf, 'Position',  [0, 0, 800, 800])
y=[];
grid=3;
x=linspace(0.,.9,grid)
y=[];
for i=1:grid
    y(end+1,:)=1-[x(i) x(i) x(i)]
end
y(2,:)=y(1,:)
colormap(y)


v=[10]
[M,c]=contourf(essen, size,matt,[0,1,1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;

set(gcf,'color','w');
ax = gca;
ax.XTick=[0.5 0.75 1];
ax.XTickLabel={'1/2', '3/4','1'};



ax.YTick=[2 20 40];
ax.FontSize = ftsz/2; 


axis([0.5 max(essen) min(size) max(size)])

lambda=0.025; %plot alternate scale of costs
cn=lambda.*(1-exp(-beta.*size.*(size-1)./2));
thresh2=(cn.*size+1)./2
hold on
plot(thresh2, size, 'Color','k', 'Linestyle','--','Linewidth',3)

