%This script generates Figure 5C. The script Figures5AandB.m must be run
%first so that the qmid.csv, smid.csv, dmid.csv and vmid.csv files all
%exist in the data directory. Filamentgrow3.m must be in the same directory


clear all
close all
phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
reps=5;


qmid=csvread("data/qmid.csv")
smid=csvread("data/smid.csv")
dmid=csvread("data/dmid.csv")
vmid=csvread("data/vmid.csv")


trials=100;
sam=10;
N=50;
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

samvar=zeros(length(phis),length(gammas));
exvar=zeros(length(phis),length(gammas));
tic
for j=1:length(gammas)
    disp(j)
    toc
    for i=1:length(phis)
        disp(i)
        gamma1=gammas(j);
        q=qmid(i,j);
        s=smid(i,j);
        d=dmid(i,j);
        v=vmid(i,j);
        
        phi=phis(i);
        
        current=[q s d v];
        params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
        res=filamentgrow3(current, params, trials);
        froles=res(:,11:end);
        p(i,j)=mean(sum(froles(:,2:end-1)==-1,2)/N);
        roles=res(:,12:11+sam-1);
        samvar(i,j)=var(sum(roles==-1,2));
        exvar(i,j)=p(i,j)*(1-p(i,j))*sam;
    end
end
coord=samvar./exvar;
figure('Position', [10 10 800 800])
nodol=find(qmid<0.02 &smid<0.001==1);
 coord2=coord;
 coord2(nodol)=-0.1;
 imagesc(phis,gammas,coord2')
  set(gcf,'color','w');    
set(gca,'YDir','normal')

col=[0 1 1; gray]
colormap(1-col)
cbh=colorbar
ax=gca
ax.CLim = [0.19 1];
box off
set(gcf,'color','w');
  xticks([0 5 10])
  yticks([0 0.5 1])
axis([-0.25 10 0.05 1.05])
set(gca,'fontsize', 25)
cbh.Ticks = [0.2 0.6 1]
cbh.TickLabels=[0.2 0.6 1]

csvwrite("data/precision.csv",coord)
