%This script generates Figures S4B, C and D. We randomly simulate 100,000
%fitness sequences for a fixed group size. We plot the 99% confidence
%intervals for the sequences in Figure S4B, we quantify how the difference
%between the optimal proportion of helpers and the optimal helper
%probability corellated with the relative asymetry of the fitness sequences
%(Figure S4C). We also quantifed how the fitness gain from adopting the
%optimal helper probability (for random specialisers) correlated with the
%difference between the optimal helper probability and the optimal
%proportion of helpers (Figure S4D).

clear all
close all
N=6; %size of the group
tests=1000000; %number of simulations
trials=30; %number of times the optimal helper probability is numerically
%solved for (attempts before moving on).
K=linspace(0,N,N+1); %grid for number of helpers
data=[];
fit1=[];
fit2=[];
%Generate the sequences. For each sequence, quantify: i) the relative asymmetry
%of the  asymmetry of the sequence, ii) the difference between the optimal
%helper probability and the optimal proportion of helpers, iii) the difference in
%fitness attained by adopting the optimal helper probability.

for i=1:tests
    display(i)
    e=rand(1,1);
    increases=rand(1,N); 
    %generate random fitness sequence
    FK=[1-e, (1-e)+e.*(cumsum(increases)/sum(increases))]; %fecundity  sequence
    w=(N-K).*FK/N; %fitness sequence
    Kopt=K(w==max(w)); %optimal number of helpers
    if Kopt>0 && Kopt<N %if division favoured for coordinated specialisers, continue 
        
        K1=linspace(0, Kopt-1, Kopt);
        K2=linspace(Kopt+1, N, N-Kopt);
        left=w(1:Kopt);
        right=w(Kopt+2:end);
        leftmass=sum(left);
        rightmass=sum(right);
        leftward=leftmass-rightmass;
        
        fit1(end+1,:)=w;
        fit2(end+1,:)=FK;
        options = optimoptions('fsolve','Display','none', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-9);
        Stochopt=0;
        BinoOpt=binopdf(K,N,Stochopt/N);
        stochoptfit=sum(BinoOpt.*w);
        for t=1:trials  %determine the optimal helper probability numerically
            guess=rand();
            point=fsolve(@dpdWK, guess,options,w);
            if abs(dpdWK(point,w))<10e-4
                temp=max(min(point,1),0);
                if temp<1 && temp>0
                    BinoOpt=binopdf(K,N,temp);
                    stochoptfit2=sum(BinoOpt.*w);
                    if stochoptfit2>stochoptfit
                        Stochopt=N*temp;
                        stochoptfit=stochoptfit2;
                    end
                end
            end
        end
        if Stochopt>0 %if division favoured, store relevant quantities
            Bino=binopdf(K,N,Kopt/N);
            BinoOpt=binopdf(K,N,Stochopt/N);
            
            stochfit=sum(Bino.*w);
            stochoptfit=sum(BinoOpt.*w);
            %quantify asymmetry of sequence
            leftward2=sum(left.*binopdf(K1,N,Kopt/N))/(sum(left.*binopdf(K1,N,Kopt/N))+sum(right.*binopdf(K2,N,Kopt/N)));
            data(end+1,:)=[(Kopt-Stochopt)/N leftward leftward2 stochoptfit-stochfit Stochopt/N w];
        end
    end
end
close all

sz=50;

ftsz=35


%Plot the 99% confidence intervals for the generated fitness sequences
%(Figure S4B)
k=K;
n=N;
f1=mean(fit1)
d1=sqrt(var(fit1));
sz=400;
SEM=d1./sqrt(length(d1));
ts = tinv([0.01  0.99],length(d1)-1);
l1=f1 + ts(1)*SEM;
u1=f1 + ts(2)*SEM;

lightgrey=0.75*[1 1 1];
darkgrey=0.35*[1 1 1];
multi=2;
f1=mean(fit1)
figure()
set(gcf, 'Position',  [1000, 0, 1000, 800])
set(gcf,'color','w');

hold on

k2 = [k/n, fliplr(k/n)];
inBetween = [l1, fliplr(u1)];
h=fill(k2, inBetween, lightgrey);
scatter(k/n,f1,sz,'k.')
set(h,'facealpha',.75)
plot(k/n,l1,'Color',darkgrey)
plot(k/n,u1,'Color',darkgrey)

ax = gca;
ax.YTick=[0 1/4 1/2]
ax.YTickLabel={'0','1/4','1/2'}
ax.XTick=[0 1/2 1]
ax.XTickLabel={'0','1/2','1'}
ax.FontSize = ftsz;
set(gcf,'color','w');

sz=10;

%Plot the how the difference between the optimal proportion of helpers and
%the optimal helper probability corellated with the relative asymetry of the
%fitness sequences(Figure S4C)
figure()
set(gcf, 'Position',  [1000, 0, 1000, 800])
set(gcf,'color','w');

line([-N N],[0 0], 'color', 'k')
hold on
line([0.5 0.5],[-N N], 'color', 'k')
scatter((data(:,3)), data(:,1), sz, darkgrey,'filled')
axis([min(data(:,3)), max(data(:,3)), min(data(:,1)), max(data(:,1))])
ax = gca;
ax.YTick=linspace(-N, N, 2*N+1)/N
ax.YTick=[-1/2 -1/4 0 1/4 1/2 3/4 1]
ax.YTickLabel={'-1/2','-1/4','0','1/4','1/2', '3/4', '1'}
ax.XTick=[-1/2 -1/4 0 1/4 1/2 3/4 1]
ax.XTickLabel={'-1/2','-1/4','0','1/4','1/2', '3/4', '1'}
ax.FontSize = ftsz;

coefficients = polyfit((data(:,3)), data(:,1), 1);
xFit = linspace(min((data(:,3))), max((data(:,3))), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'k-', 'LineWidth', 2);
grid off;


%Plot the how how the fitness gain from adopting the
%optimal helper probability (for random specialisers) correlated with the
%difference between the optimal helper probability and the optimal
%proportion of helpers (Figure S4D).
figure()
set(gcf, 'Position',  [1000, 0, 1000, 800])
set(gcf,'color','w');

line([0 1],[0 0], 'color', 'k')
hold on
line([0 0],[0 max(data(:,4))], 'color', 'k')
scatter(abs(data(:,1)), data(:,4), sz, darkgrey,'filled')
axis([min(abs(data(:,1))), max(abs(data(:,1))), min(data(:,4)), max(data(:,4))])
ax = gca;
ax.YTick=[1/1000 1/30 2/30 3/30]
ax.YTickLabel={'0','1/20','2/20', '3/20'}
ax.XTick=[-1/2 -1/4 0 1/4 1/2 3/4 1]
ax.XTickLabel={'-1/2','-1/4','0','1/4','1/2', '3/4', '1'}
ax.FontSize = ftsz;
grid off
coefficients = polyfit(abs(data(:,1)), data(:,4), 1);
xFit = linspace(min(abs(data(:,1))), max(abs(data(:,1))), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'k-', 'LineWidth', 2);
grid off;