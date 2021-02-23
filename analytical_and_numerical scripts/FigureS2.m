%This script generates Figure S2. It creates maps of parameter space (group size and essentiality
%of cooperation) and evaluates which mechanism (coordinated specialisation
%or random specialisaiton) is favoured. We assume here the reproductive
%fecundity is decelerating (accelerating) in the proportion of helpers in the group.

clear all
close all

Nmax=50;
Nmin=2;
HBmax=10;
HBmin=1;
size=linspace(Nmin,Nmax,(Nmax-Nmin)+1); %discretisation for group size
HBr=linspace(1,HBmax,length(size))
HBl=fliplr(1./HBr(2:end))
HB=[HBl HBr] %discretisation for relative importance of cooperation 
alpha=1/(3/4)  %shape of return from proportion of helpers (decelerating when alpha<1)
ftsz=50

%Iterate through each possible parameter combination and evaluate 
%the fitness of fully coordinated groups and random groups. 
%Store in results. 
results=[]
matt0=zeros(length(size), length(HB));
matt1=zeros(length(size), length(HB));
matt2=zeros(length(size), length(HB));
matt3=zeros(length(size), length(HB));
matt4=zeros(length(size), length(HB));
matt5=zeros(length(size), length(HB));
for i=1:length(size)
    disp(i)
    for j=1:length(HB)
        hb=HB(j);
        N=size(i);
        scale=1000;
        offset=0.25;
        Kthin=linspace(0,N,1000);
        K=linspace(0,N,N+1);
        trials=20;
        Wbroad=((1-(K/N))).*(1+hb*((K./(N)).^alpha));
        Kopt=K(max(Wbroad)==Wbroad);
        popt=mean(Kopt)/N;
        coordfit=max(Wbroad);
        stochfit=WK(mean(Kopt)/N,Wbroad);
        results(end+1,:)=[hb N mean(Kopt)/N coordfit stochfit];
        matt1(i,j)=coordfit;
        matt3(i,j)=stochfit;
        matt4(i,j)=popt;
    end
end


for i=1:length(HB)
    m11=(1-2E-2)*matt1(:,i);
    m12=(1-4E-2)*matt1(:,i);
    m2=matt3(:,i);
    popt=matt4(:,i);
    if length(find(m11>m2))>0
        track1(i)=max(size(find(m11>m2)));
    else
        track1(i)=0;
    end
    if length(find(m12>m2))>0
        track2(i)=max(size(find(m12>m2)));
    else
        track2(i)=0;
    end
    if m2>m12 & length(find(popt<=0))>0
        nodol(i)=max(size(find(popt<=0)));
    else
        nodol(i)=0;
    end
end
curve2=track2;
curve1=zeros(1,length(track2));
x=log(HB);
figure()
set(gcf, 'Position',  [0, 0, 800, 800])

plot(x, curve1, 'k', 'LineWidth', 1);
hold on;
plot(x, curve2, 'k', 'LineWidth', 1);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k ');
hold on
plot(x,track1,'k--')

box off
  set(gcf,'color','w');    
ax = gca;
ax.XTick=[log(1/HBmax) 0 log(HBmax)];
ax.XTickLabels={"1/"+HBmax "1" HBmax};
ax.YTick=[Nmin Nmax/2 Nmax];
ax.FontSize = ftsz/2; 
axis([log(1/HBmax) log(HBmax) min(size) max(size)])



grid=length(HB);
x=log(HB)
curve1=zeros(1,grid);
curve2=nodol;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);
fill(x2, inBetween, 'r ');
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);

plot(x,track1,'k--')

