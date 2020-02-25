%This script generates Figures S1B and S1C. We create a map of parameter space 
%(group size and essentiality of cooperation) and evaluate which mechanism 
%(fully coordinated specialisation or random specialisaiton) is favoured. We
%output figures for the cases where 1) the optimal proportion of
%helpers is constrained by the possible allocations of labour (not continuous), and 2) the
%optimal probability of becoming a helper is optimised for random groups. 

clear all
close all

%Parameters
Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1); %discretisation for group size
essen=linspace(0,1,length(size)); %discretisation for essentiality of cooperation 
beta=0.01; %shape of coordination cost
lambda=0.1; %scale of coordination costs


%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers as well as the optimal helper probability. Evaluate 
%the fitness of a fully coordinated group, and random groups for each of these 
%scenarios. Store in results. 
results=[];
for i=1:length(essen)
    disp(i)
    for j=1:length(size)
        e=essen(i);
        N=size(j);
        K=linspace(0,N,N+1);
        p=linspace(0.01,0.99,100);
        Wbroad=((1-(K./N))).*(1-e+e*((K./(N))));
        Kopt=K(max(Wbroad)==Wbroad);
        
        qopt=max(0,1+1/(2*(N-1))-N/(2*e*(N-1)));
        coordfit=(1-lambda.*(1-exp(-beta.*N.*(N-1)/2))).*max(Wbroad);
        stochdiscfit=WK(mean(Kopt)/N,Wbroad);
        stochoptfit=WK(qopt,Wbroad);
        results(end+1,:)=[e N mean(Kopt)/N qopt coordfit stochdiscfit stochoptfit];
    
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
matt1=[];
matt2=[];
matt3=[];
es=[];
ne=[];
for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        coordfit=results(results(:,2)==N&results(:,1)==e,5);
        stochdiscfit=results(results(:,2)==N&results(:,1)==e,6);
        stochoptfit=results(results(:,2)==N&results(:,1)==e,7);
        matt1(i,j)=stochdiscfit/coordfit;
        matt2(i,j)=stochoptfit/coordfit;
    end
end

%Optimal mechanism when the probability/poportion of helpers is set
%as the discrete optimal proportion of helpers
figure()
set(gcf, 'Position',  [0, 0, 800, 800])
colormap(y)
v=[0.00001 0.015 0.025 1]
[M,c]=contourf(essen, size,matt1,[0, 1, 1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;
set(gcf,'color','w');
ax = gca;
ax.XTick=[0 0.5 1];
ax.YTick=[2 20 40];
ax.XTickLabel={'0', '1/2','1'};
ax.FontSize = ftsz/2; 
axis([0 max(essen) min(size) max(size)])
hold on




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


%Optimal mechanism when the probability of being a helper is equal to the
%value that optimises the fitness of a random group
figure()
set(gcf, 'Position',  [0, 0, 800, 800])
colormap(y)
v=[0.00001 0.015 0.025 1]
[M,c]=contourf(essen, size,matt2,[0, 1, 1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;
set(gcf,'color','w');
ax = gca;
ax.XTick=[0 0.5 1];
ax.YTick=[2 20 40];
ax.XTickLabel={'0', '1/2','1'};
ax.FontSize = ftsz/2; 
axis([0 max(essen) min(size) max(size)])





endol=results(results(:,4)<=0,1);
nndol=results(results(:,4)<=0,2);
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


