% plot the data that will generate figure 6A and 6B of the main text
% clump.csv and clumpcost.csv as produced by Figure6AandB1.m are required
% in the root directory
clear all
close all

phis=linspace(0,9.5,20);
gammas=linspace(0.1,1,10);
clump=csvread("clump.csv")
clumpcost=csvread("clumpcost.csv")

figure('Position', [10 10 800 800])
imagesc(phis',gammas',clump')
cbh=colorbar
colormap(1-gray)
  box off
set(gcf,'color','w');
set(gca,'YDir','normal')
xticks([0 5 10])
yticks([0 0.5 1])
axis([0 10 0.1 1])
caxis([0 3])
set(cbh,'XTick',[0 1 2 3])
  set(gca,'fontsize', 50)  
  set(gca,'linewidth',3)
set(gca,'TickDir','out')
  
  
figure('Position', [10 10 800 800])
imagesc(phis',gammas',-clumpcost')
cbh=colorbar

colormap(1-gray)
  box off
set(gcf,'color','w');
set(gca,'YDir','normal')
xticks([0 5 10])
yticks([0 0.5 1])
axis([0 10 0.1 1])
caxis([0 0.3])
set(cbh,'XTick',[0 0.15 0.3])
  set(gca,'fontsize', 50)
set(gca,'linewidth',3)
set(gca,'TickDir','out')
