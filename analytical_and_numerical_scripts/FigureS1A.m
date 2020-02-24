%This script generates Figure S1A. We produce a plot that shows how the 
%discrete optimal proportion of helpers compares to the continuous 
%approximation optimal proportion of helpers and the optimal probability 
%of helpers as group size increases.
clear all
close all

%Parameter space
Nmax=40;
Nmin=2;
size=linspace(Nmin,Nmax,(Nmax-Nmin)+1); %discretisation for group size
e=0.75 %set essentiality to 3/4

qopt=[]; %optimal helper probability
popt=[]; %continuous approximation for optimal proportion of helpers
discopt=[]; %optimal proportion of helpers (sans continuous approximation)

%Iterate through group size and calculate the optimal
%proportion of helpers as well as the optimal helper probability (discrete and continuous). 
for j=1:length(size)
    N=size(j);
    K=linspace(0,N,N+1);
    p=linspace(0.01,0.99,100);
    Wbroad=((1-(K./N))).*(1-e+e*((K./(N))));
    Kopt=K(max(Wbroad)==Wbroad);
    popt(end+1)=max(0,(2*e-1)/(2*e))
    qopt(end+1)=max((2*N*e-N-e)./(2.*e.*(N-1)),0);
    discopt(end+1)=mean(Kopt)/N;
    
end
%Generate plot
figure()
set(gcf, 'Position',  [0, 0, 800, 800])
lw=3;
sz=150;
ftsz=50;

scatter(size,discopt,sz,'ko');
hold on
plot(size, popt,'k--','Linewidth',lw)
plot(size, qopt,'k-','Linewidth',lw)
axis([Nmin Nmax 0 0.5])

set(gcf,'color','w');
ax = gca;
ax.XTick=[2 20 40];
ax.YTick=[0 0.25 0.5]
ax.YTickLabel={'0', '1/4','1/2'};
ax.FontSize = ftsz/2; 
