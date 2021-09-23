%Generates Figure S7 of the Supp Info, showing how growing groups are
%expected to lead to larger clump sizes than non-growing groups.
%qmid.csv, smid.csv, dmid.csv, vmid.csv and filamentgrow2.csv are required
%in the root directory
clear all
close all
phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
reps=5;

%load optimal helper differentiation probabilities from previous simulations
qmid=csvread("qmid.csv")


%We first go run simulations of growing groups and store the average clump
%sizes across simulations
i=1;
j=1;
gamma1=gammas(j);
q=qmid(i,j);
phi=phis(i);
trials=1000;
L=50;
N=L;
theta=50;
cap1=100;
noise=0.01;
lamb=30;
lamb1=2.5;
lamb2=0.5;
lamb3=0.25;
growthmax=10;
gamma3=1;
gamma2=0.5;
start=[q 0 0 0];
params=[L phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
avcs1=zeros(1,trials);
fit1=zeros(1,trials);
if q>0.05
    for t=1:trials
        disp(t)
        [res1 stranded clumped cfit rolesc time replog grolog]=filamentgrow2(start,params);
        avcs1(t)=sum(clumped.*linspace(1,N,N)')/sum(clumped);
        fit1(t)=cfit;
    end
    c = polyfit(avcs1,fit1,1);
end

%Plot results for growing groups
figure()
scatter(avcs1,fit1)
grid=100;
c = polyfit(avcs1,fit1,1);
x=linspace(min(avcs1),max(avcs1),grid);
y_est = polyval(c,x);
hold on
plot(x,y_est,'r--','LineWidth',2)
figure()
edges=linspace(0.5,4.5,10)
histogram(avcs1,edges)
mavcs=mean(avcs1)
line([mavcs mavcs], [0, 120],'Color','r','LineWidth',4)

growthmax=10;
gamma3=1;
gamma2=0.5;

maxstrand=50;
start=[q 0 0 0];


%We now run simulations of non-growing groups and store average clump size
%per simulation

trials=500;
for t=1:trials
    roles=binornd(1,q,50,1);
    heps=find(roles==1);
    reps=find(roles==0);
    %map of all contiguous connections between helpers
    transfer=eye(L,L);
    transfer(1:end-1,2:end)=transfer(1:end-1,2:end)+eye(L-1,L-1);
    transfer(2:end,1:end-1)=transfer(2:end,1:end-1)+eye(L-1,L-1);
    for i=1:length(reps) %break connections wherever there is a reproductive
        transfer(reps(i),:)=zeros(L,1);
        transfer(:,reps(i))=zeros(L,1)';
    end
    %iterate through each helper and report size of the helper clump it is
    %in
    strands=zeros(length(heps),1);
    for i=1:length(heps)
        hep=heps(i);
        start=zeros(L,1);
        start(hep)=1;
        strands(i)=sum(sum((start.*(transfer)^maxstrand)>0));
    end
    
    %record the number of helper clumps of each possible size from 1 to
    %maxstrand
    hepstrand=zeros(maxstrand,1);
    for i=1:maxstrand
        hepstrand(i)=sum(strands==i)/i;
    end
    if isnan(mean(hepstrand))
        take=hepstrand;
    end
    avcs0(t)=sum(hepstrand.*linspace(1,L,L)')/sum(hepstrand);
end


%Plot results for non-growing groups
figure('Position', [10 10 800 800])
edges=linspace(1,5,20)
subplot(2,1,1)
histogram(avcs0,edges,'Normalization','probability','FaceColor','k')



%Plot results for growing groups
figure('Position', [10 10 800 800])
histogram(avcs1,edges,'Normalization','probability','FaceColor','k')
mavcs=mean(avcs1)
line([mavcs mavcs], [0, 120],'Color','k','LineWidth',4,'LineStyle','--')

box off
set(gcf,'color','w');
xticks([1 2 3 4 5])
yticks([0 0.125 0.25])
set(gca,'fontsize', 50)
set(gca,'linewidth',3)
set(gca,'TickDir','out')
axis([1 5 0 0.25])
hold on
figure('Position', [10 10 800 800])
histogram(avcs0,edges,'Normalization','probability','FaceColor','k')
mavcs=mean(avcs0)
line([mavcs mavcs], [0, 120],'Color','k','LineWidth',4,'LineStyle','--')

set(gca,'fontsize', 50)
axis([1 5 0 0.25])
box off
set(gcf,'color','w');
xticks([1 2 3 4 5])
yticks([0 0.125 0.25])
set(gca,'linewidth',3)
set(gca,'TickDir','out')