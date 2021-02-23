%This script generates Figure 5A and 5B. simrun.m must be run for all of
%the parameter combinations in the figure and the output stored in the
%folder data which must be in the same directory as Figures5AandB.m. The

clear all
close all

phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
reps=5;
q1=zeros(1,reps);
s1=zeros(1,reps);
d1=zeros(1,reps);
v1=zeros(1,reps);
f1=zeros(1,reps);

q=zeros(length(phis),length(gammas));
s=zeros(length(phis),length(gammas));
d=zeros(length(phis),length(gammas));
v=zeros(length(phis),length(gammas));
f=zeros(length(phis),length(gammas));
for j=1:length(gammas)
for i=1:length(phis)

for r=1:reps
dat=csvread("data/track_"+phis(i)+"_"+gammas(j)+"_"+r+".csv");
q1(r)=mean(dat(:,1));
s1(r)=mean(dat(:,2));
d1(r)=mean(dat(:,3));
v1(r)=mean(dat(:,4));
f1(r)=mean(dat(:,4));
end
q(i,j)=mean(q1);
s(i,j)=mean(s1);
d(i,j)=mean(d1);
v(i,j)=mean(v1);
f(i,j)=mean(f1);

end
end

ftsz=40

figure()
set(gcf, 'Position',  [0, 0, 800, 800])
imagesc(phis,gammas,s')
colormap(1-gray)
set(gca,'YDir','normal')
box off
cbh=colorbar
  set(gcf,'color','w');    
box off
set(gcf,'color','w');
  xticks([0 5 10])
  yticks([0 0.5 1])
axis([-0.25 10 0.05 1.05])
set(gca,'fontsize', ftsz)
caxis([0 0.07])
cbh.Ticks = [0 0.035 0.07]
cbh.TickLabels=[0 0.035 0.07]

figure()
set(gcf, 'Position',  [0, 0, 800, 800])
imagesc(phis,gammas,v')
colormap(1-gray)
set(gca,'YDir','normal')
box off
cbh=colorbar
  set(gcf,'color','w');    
box off
set(gcf,'color','w');
  xticks([0 5 10])
  yticks([0 0.5 1])
axis([-0.25 10 0.05 1.05])
set(gca,'fontsize', ftsz)
caxis([0 4])
cbh.Ticks = [0 2 4]
cbh.TickLabels=[0 2 4]


figure()
set(gcf, 'Position',  [0, 0, 800, 800])
imagesc(phis,gammas,q')
colormap(1-gray)
set(gca,'YDir','normal')
box off
cbh=colorbar
  set(gcf,'color','w');    
box off
set(gcf,'color','w');
  xticks([0 5 10])
  yticks([0 0.5 1])
axis([-0.25 10 0.05 1.05])
set(gca,'fontsize', ftsz)
caxis([0 0.5])
cbh.Ticks = [0 0.25 0.5]
cbh.TickLabels=[0 0.25 0.5]


figure()
set(gcf, 'Position',  [0, 0, 800, 800])
imagesc(phis,gammas,d')
colormap(1-gray)
set(gca,'YDir','normal')
box off
cbh=colorbar
  set(gcf,'color','w');    
box off
set(gcf,'color','w');
  xticks([0 5 10])
  yticks([0 0.5 1])
axis([-0.25 10 0.05 1.05])
set(gca,'fontsize', ftsz)
caxis([0 1.1])
cbh.Ticks = [0 0.55 1.1]
cbh.TickLabels=[0 0.55 1.1]

csvwrite("data/qmid.csv",q)
csvwrite("data/smid.csv",s)
csvwrite("data/vmid.csv",v)
csvwrite("data/dmid.csv",d)
