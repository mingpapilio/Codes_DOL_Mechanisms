%Generates Figure 6C and 6D of the main text.
clear all
close all
phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
reps=5;


qmid=csvread("qmid.csv")
smid=csvread("smid.csv")
dmid=csvread("dmid.csv")
vmid=csvread("vmid.csv")


%We first consider the case where cells specialise randomly
i=1;
j=1;
gamma1=gammas(j);
q=qmid(i,j);
s=smid(i,j);
v=vmid(i,j);
d=dmid(i,j);
phi=phis(i);
trials=50;
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
start=[q 0 0 0];
params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
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



%We now consider the case where cells coordinate labour
gamma1=gammas(j);
q=qmid(i,j);
phi=phis(i);
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
start=[q s d v];
params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
avcs2=zeros(1,trials);
fit2=zeros(1,trials);
if q>0.05
    for t=1:trials
        disp(t)
        [res1 stranded clumped cfit rolesc time replog grolog]=filamentgrow2(start,params);
        avcs2(t)=sum(clumped.*linspace(1,N,N)')/sum(clumped);
        fit2(t)=cfit;
    end
    c = polyfit(avcs2,fit2,1);
end


%Plot results
figure()
scatter(avcs2,fit2)
grid=100;
c = polyfit(avcs2,fit2,1);
x=linspace(min(avcs2),max(avcs2),grid);
y_est = polyval(c,x);
hold on
plot(x,y_est,'r--','LineWidth',2)
figure()
edges=linspace(0.5,4.5,10)
histogram(avcs2,edges)
mavcs=mean(avcs2)
line([mavcs mavcs], [0, 120],'Color','r','LineWidth',4)


close all
figure('Position', [10 10 800 800])
lw=5;
mksz=20;
scatter(avcs1,fit1,mksz,'filled','MarkerFaceColor',[1 0.25 0.25],'MarkerFaceAlpha',0.5)
hold on
scatter(avcs2,fit2,mksz,'filled','MarkerFaceColor',[0.25 0.25 1],'MarkerFaceAlpha',0.5)
grid=100;
c = polyfit(avcs1,fit1,1);
x=linspace(min(avcs1),max(avcs1),grid);
y_est = polyval(c,x);
hold on
plot(x,y_est,'k--','LineWidth',lw)
c = polyfit(avcs2,fit2,1);
x=linspace(min(avcs2),max(avcs2),grid);
y_est = polyval(c,x);
hold on
plot(x,y_est,'k--','LineWidth',lw)
box off
set(gcf,'color','w');
xticks([1 2 3 4 5])
yticks([0 1 2])
axis([1 5 0 2])
set(gca,'fontsize', 25)
edges=linspace(1,5,31)

figure('Position', [10 10 800 800])
histogram(avcs2,edges,'Normalization','probability','FaceColor',[0.25 0.25 1])

box off
set(gcf,'color','w');
xticks([1 2 3 4 5])
yticks([0 0.5 1])

axis([1 5 0 0.5])
line([mean(avcs2) mean(avcs2)],[0 0.5],"Color",'b',"LineWidth",4)
hold on
histogram(avcs1,edges,'Normalization','probability','FaceColor',[1 0.25 0.25])
line([mean(avcs1) mean(avcs1)],[0 0.5],"Color",'r',"LineWidth",4)

set(gca,'fontsize', 25)

