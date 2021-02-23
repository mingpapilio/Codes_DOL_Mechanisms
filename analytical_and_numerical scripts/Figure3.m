%This file generates Figure 3 of the main text. It discretises the relative
%importance of cooperation (h/b) and for two different values of the relative
%cost of coordination (gamma), it calculated the upper bound on the size of
%the group (n) that is required for coordinated specialisation to be
%favoured. 
close all
clear all

hbupper=10; %upper bound of hb
lw=3  %plot params
ftsz=50



hb=linspace(1,hbupper,100) %range of relative importanct of cooperation considered

C=0.002; %higher relative cost of coordination
line1=(hb-1)./((hb+1).*C)
bot=zeros(1,length(line1));
figure('Position', [10 10 800 800])
set(gca,'Color','w')
plot(hb,line1,'k-','LineWidth', lw)


x2=[hb, fliplr(hb)];
inBetween=[bot,fliplr(line1)];
fill(x2,inBetween, 'k')

C=0.001 %lower relative cost of coordination
hold on
line2=(hb-1)./((hb+1).*C)
plot(hb,line2,'k--','LineWidth', lw)
set(gcf,'Color','w')

axis([1 hbupper 0 1000])
xticks([1 hbupper/2 hbupper])
yticks([2 500 1000])

box off
set(gca,'fontsize', ftsz);