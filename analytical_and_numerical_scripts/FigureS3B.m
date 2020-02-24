%This script generate figure S3B. We create a maps of parameter space (group size and essentiality
%of cooperation) and evaluate which mechanism (fully coordinated specialisation
%or random specialisation) is favoured. We assume here the reproductive
%fecundity depends linearly on the ratio of helpers to reproductives in the group.

clear all
close all

%Parameter space
Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1);  %discretisation for group size
essen=linspace(0.5,1,length(size)); %discretisation for essentiality of cooperation 

beta=0.01; %shape of coordination cost
lambda=0.1; %scale of coordination costs

%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers (discrete and continuous approximation) as well as
%the optimal helper probability. Evaluate the fitness of a coordinated
%group, and random groups for each of these scenarios. Store in results. 
results=[]

for i=1:length(essen)
    disp(i)
    for j=1:length(size)
        e=essen(i);
        N=size(j);
        K=linspace(0,N,N+1);
        p=linspace(0.01,0.99,100);
        Wbroad=(1-(K./N)).*(1-e)+e*(K/N);
        Wbroad(end)=0;
        Kopt=K(max(Wbroad)==Wbroad);
        
        popt=(N-1)./N;
        qopt=max(0,((2*e-1)/(N*e))^(1/(N-1)));
        coordfit=(1-lambda.*(1-exp(-beta.*N.*(N-1)/2))).*max(Wbroad);
        stochcontfit=WK(popt,Wbroad);
        stochdiscfit=WK(mean(Kopt)/N,Wbroad);
        stochoptfit=WK(qopt,Wbroad);
        results(end+1,:)=[e N mean(Kopt)/N qopt popt coordfit stochcontfit stochdiscfit stochoptfit];
    
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

matt3=[];
for i=1:length(size)
    for j=1:length(essen)
        e=essen(j);
        N=size(i);
        coordfit=results(results(:,2)==N&results(:,1)==e,6);
        stochcontfit=results(results(:,2)==N&results(:,1)==e,7);
        stochdiscfit=results(results(:,2)==N&results(:,1)==e,8);
        stochoptfit=results(results(:,2)==N&results(:,1)==e,9);
        matt3(i,j)=stochoptfit/coordfit;
       
    end
end


%Optimal mechanism when the probability of being a helper is equal to the
%value that optimises the fitness of a random group
figure()
set(gcf, 'Position',  [0, 0, 800, 800])
colormap(y)
v=[0.00000001 0.075 0.125 1]
[M,c]=contourf(essen, size,matt3,[0, 1, 1000],'Color','k','ShowText','off', 'Linestyle','-')
c.LineWidth = 2;
set(gcf,'color','w');
ax = gca;
ax.XTick=[0 0.5 1];
ax.XTickLabel={'0', '1/2','1'}
ax.YTick=[2 20 40];
ax.FontSize = ftsz/2; 
axis([0 max(essen) min(size) max(size)])


