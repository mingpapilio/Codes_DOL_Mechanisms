%Generate the data that will produce Figure 6A and B
%filamentgrow2.m and qmid.csv are required in the root directory

clear all
close all
phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
reps=5;


qmid=csvread("qmid.csv")


trials=5;
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
cost=zeros(length(phis), length(gammas))
clump=zeros(length(phis), length(gammas))

tic
for j=1:length(gammas)
    disp(j)
    toc
    for i=1:length(phis)
        disp(i)
        gamma1=gammas(j);
        q=qmid(i,j);
        phi=phis(i);
        
        start=[q 0 0 0];
        params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
        if q>0.05
        for t=1:trials
            [res1 stranded clumped cfit rolesc time replog grolog]=filamentgrow2(start,params);
            avcs(t)=sum(clumped.*linspace(1,N,N)')/sum(clumped);
            maxclump(t)=max(find(clumped>0));
            fit(t)=cfit;
        end
        c = polyfit(avcs,fit,1);
        cost(i,j)=c(1);
        clump(i,j)=mean(avcs);
        end
    end
end
csvwrite("clumpcost.csv",cost)
csvwrite("clump.csv",clump)
