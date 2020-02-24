%This script generates Figure S4A. We plot an arbitrary asymmetric fitness
%function in order to illustrate that the optimal helper probability for
%random specialisers may differ from the optimal proportion of helpers. We
%do not actually calculate the optimal helper probability for this fitness
%function. We use this figure as an exaggerated 'cartoon' to demonstrate a
%concept.

clear all
close all
N=6; %group size
alpha=0.25 %sensitivity to number of helpers
beta=1 %seneitivity to number of reproductives


%Arbitrary fitness function
Kthin=linspace(0,N,100*N+1)
Kbroad=linspace(0,N,N+1)
Wthin=(Kthin.^alpha).*((1-Kthin/N)).^beta
Wbroad=(Kbroad.^alpha).*((1-Kbroad/N)).^beta


%Generate plot
figure()
sz=50;
ftsz=40
set(gcf, 'Position',  [1000, 0, 1000, 800])

plot(Kthin/N,(Wthin),'k', 'LineWidth', 2)
set(gcf,'color','w');
hold on
h=scatter(Kbroad/N,((Wbroad)),300,'k', 'filled')
axis([0 1 0 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
box off

Kopt=Kbroad(Wbroad==max(Wbroad))

line([Kopt/N Kopt/N], [0 max(Wbroad)],'Color','k', 'LineWidth', 3) %optimal proporiton of helpers
line([Kopt/N+1.9/N Kopt/N+1.9/N], [0 max(Wbroad)],'Color','k', 'LineWidth', 3,'LineStyle', '--') %optimal helper probability
%note-> we do not actually calculate the optimal helper probability for
%this fitness function

ax = gca;
ax.YTick=[0 1/2 1]
ax.YTickLabel={'0','1/2','1'}

ax.XTick=[0 1/2 1]
ax.XTickLabel={'0','1/2','1'}
ax.FontSize = ftsz; 