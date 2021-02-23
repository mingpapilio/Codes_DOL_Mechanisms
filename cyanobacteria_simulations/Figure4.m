%This file generates Figure 4 of the main text. Must be run in the same
%directory as filamentgrow2.m
clear all 
close all



%Generate Figure 4A
offf=20
figure('Position', [10 10 800 800])
ftsz=50;
axw=4
lw=4
sz=1000;
col1=[0 0 0]
col2=[0.76 0.76 0.76]
reps=[2 3 4 5];
heps=[1 6];
height=5;
offset=0.1;
hold on
N=length(reps)+length(heps);
cells=linspace(1,N,N)';
lamb=30;
gamma2=0.5;
signal=0.4;
gradflux=zeros(N, 1)
if length(heps)>0
    for v=1:length(heps)
        k=heps(v);
        dist=abs(cells-k);
        gradflux=gradflux+lamb*signal*(gamma2.^dist);
    end
end


hscale=2;
grad=gradflux/max(gradflux);
transp=0.4;
for i=1:length(reps)
    colb=1-[0 grad(reps(i)) grad(reps(i)) transp];
    rectangle('Position',[reps(i) height-hscale/4 1 hscale], 'FaceColor',colb, 'EdgeColor', colb)
end
for i=1:length(heps)
    colb=1-[0 grad(heps(i)) grad(heps(i)) transp];
    rectangle('Position',[heps(i) height-hscale/4 1 hscale], 'FaceColor',colb, 'EdgeColor',colb)
end

line([1+offset N+1-offset], [height+hscale/4 height+hscale/4], 'Color', 'k', 'LineWidth', lw)

for i=1:length(reps)
    colb=1-[0 grad(reps(i)) grad(reps(i)) transp];
    rectangle('Position',[reps(i)+offset/2  height+offset/2  1-offset 1-offset], 'FaceColor',col2, 'EdgeColor', [0 0 0],'Curvature',[1 1])
end
for i=1:length(heps)
    colb=1-[0 grad(heps(i)) grad(heps(i)) transp];
    rectangle('Position',[heps(i)+offset/2 height+offset/2  1-offset 1-offset], 'FaceColor',col1, 'EdgeColor', [0 0 0],'Curvature',[1 1])
end



reps=[2 3 5 6];
heps=[1 7];
new=[4];
%line([1 max([reps,heps])], [1 1],'Color',col1)
hold on

N=length(reps)+length(heps)+length(new);
cells=linspace(1,N,N)';
gradflux=zeros(N, 1)
if length(heps)>0
    for v=1:length(heps)
        k=heps(v);
        dist=abs(cells-k);
        gradflux=gradflux+lamb*signal*(gamma2.^dist);
    end
end


grad=gradflux/max(gradflux);
for i=1:length(reps)
    colb=1-[0 grad(reps(i)) grad(reps(i)) transp];
    rectangle('Position',[reps(i) 1-hscale/4 1 hscale], 'FaceColor',colb, 'EdgeColor', colb)
end
for i=1:length(heps)
    colb=1-[0 grad(heps(i)) grad(heps(i)) transp];
    rectangle('Position',[heps(i) 1-hscale/4 1 hscale], 'FaceColor',colb, 'EdgeColor',colb)
end
colb=1-[0 grad(new) grad(new) transp];
rectangle('Position',[new 1-hscale/4 1 hscale], 'FaceColor',colb, 'EdgeColor', colb)

line([1+offset N+1-offset], [1+hscale/4 1+hscale/4], 'Color', 'k', 'LineWidth', lw)

for i=1:length(reps)
    colb=1-[0 grad(reps(i)) grad(reps(i)) transp];
    rectangle('Position',[reps(i)+offset/2  1+offset/2  1-offset 1-offset], 'FaceColor',col2, 'EdgeColor', [0 0 0],'Curvature',[1 1])
end
for i=1:length(heps)
    colb=1-[0 grad(heps(i)) grad(heps(i)) transp];
    rectangle('Position',[heps(i)+offset/2 1+offset/2  1-offset 1-offset], 'FaceColor',col1, 'EdgeColor', [0 0 0],'Curvature',[1 1])
end
colb=1-[0 grad(new) grad(new) transp];
rectangle('Position',[new+offset/2 1+offset/2  1-offset 1-offset], 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0],'Curvature',[1 1])

parent=3;
line([parent+0.5 parent+0.5], [height+offset/2  1+hscale/2-offset/2], 'Color', 'k', 'LineWidth', lw/2, 'LineStyle','--')
line([parent+0.5 new+0.5], [height+offset/2  1+hscale/2-offset/2], 'Color', 'k', 'LineWidth', lw/2, 'LineStyle','--')

%box off
set(gca,'visible','off')
set(gcf,'color','w');
title("A)")

axis([0 max([reps,heps])+4 0 max([reps,heps])+1])



%Generate Figure 4B
axw=2
lw=5
figure('Position', [10 10 800 800])
grid=100000;
q=0.5;
p=linspace(0,10,grid);
senses0=[0 0.2 0.75];
optimp=5;
hold on
for i=1:length(senses0)
    sense=senses0(i);
    G=1-2./(1+exp(-(sense).*(p-optimp)));
    hold on
    if i==1
    plot(p,max(min(q+G,1),0), ':','LineWidth',lw,'Color','k')
    end
    if i==2
    plot(p,max(min(q+G,1),0), '--','LineWidth',lw,'Color','k')
    end
    if i==3
    plot(p,max(min(q+G,1),0), 'LineWidth',lw,'Color','k')
    end
    
end
set(gcf,'color','w');


set(gca,'ytick',[0 q 1])
set(gca,'yticklabels',{'0' 'q' '1'})

set(gca,'xtick',[0 optimp])
set(gca,'xticklabels',{'0', 'd'})
set(gca,'Fontsize',ftsz);
set(gca,'TickDir','out'); 
ax = gca;
ax.LineWidth = axw;

%Generate Figure 4C
figure('Position', [10 10 800 800])
trials=10000;

N=50;
sz=120

phi=1;
theta=10;
cap1=100;
noise=0.01;

lamb=30;


lamb1=1;
lamb2=1;
lamb3=0.05;
growthmax=10;
gamma3=1;

gamma1=0.5;
gamma2=0.5;

q=0.33;
s=0;
h=1;
v=0;
params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
roles=filamentgrow23([q s h v],params)
scal=150;

for n=1:N-3
    nowRoll=roles(:,n);
    heps=find(nowRoll==-1);
    reps=find(nowRoll==1);
    line([1 max([reps',heps'])], [scal*n scal*n],'Color',col1)
    hold on
    scatter(reps, scal*n*ones(1,length(reps)), sz,'filled','MarkerFaceColor',col2)
    hold on
    if length(heps)>0
        scatter(heps, scal*n*ones(1,length(heps)), sz,'filled','MarkerFaceColor',col1)
    end
end
axis([0 N 0 scal*N+offf])

set(gcf,'color','w');
ylabel("Time")
xlabel("Cell position in filament")
set(gca, 'YDir','reverse')
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'Visible','off')



%Generate Figure 4D
figure('Position', [10 10 800 800])

phi=1;
theta=10;
cap1=100;
noise=0.01;

lamb=30;


lamb1=1.5;
lamb2=1;
lamb3=0.05;
growthmax=10;
gamma3=1;

gamma1=0.5;
gamma2=0.5;

s=0.1;
h=1;
v=1.5;
params=[N phi theta gamma1 gamma2 cap1 growthmax gamma3 noise lamb lamb1 lamb2 lamb3];
roles=filamentgrow23([q s h v],params)
for n=1:N-3
    nowRoll=roles(:,n);
    heps=find(nowRoll==-1);
    reps=find(nowRoll==1);
    line([1 max([reps',heps'])], [scal*n scal*n],'Color',col1)
    hold on
    scatter(reps, scal*n*ones(1,length(reps)), sz,'filled','MarkerFaceColor',col2)
    hold on
    if length(heps)>0
        scatter(heps, scal*n*ones(1,length(heps)), sz,'filled','MarkerFaceColor',col1)
    end
end
axis([0 N 0 scal*N+offf])
set(gcf,'color','w');
ylabel("Time")
xlabel("Cell position in filament")
set(gca, 'YDir','reverse')
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'Visible','off')%

