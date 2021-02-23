%This script generates figure S3A. We create a map of parameter space (group size and essentiality
%of cooperation) and evaluate which mechanism (fully coordinated specialisation
%or random specialisation) is favoured. We assume here the reproductive
%fecundity depends linearly on the number of helpers in the group.

clear all

%Parameter space

Nmax=1000;
Nmin=2;
HBmax=50;
HBmin=0;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1); %discretisation for group size
HBr=linspace(1,HBmax,length(size))
HBl=fliplr(1./HBr(2:end))
HB=[HBl HBr]%discretisation for relative importance of cooperation 
ftsz=50;

%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers. Evaluate the fitness of a coordinated
%group, and random groups and determine relative fitness
matt=zeros(length(size),length(HB));
matt2=zeros(length(size),length(HB));

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
        Wbroad=((1-(K/N))).*(1+hb*((K)));
        Kopt=(hb*N-1)/(2*hb*N);
        popt=(Kopt)/N;
        coordfit=((1-(Kopt/N))).*(1+hb*((Kopt))); 
        stochfit=((1-(Kopt/N))).*(1+hb*((Kopt)))-hb*(Kopt/N)*(1-Kopt/N);
        matt1(i,j)=coordfit;
        matt3(i,j)=stochfit;
        matt4(i,j)=popt;
    end
end


for i=1:length(HB)
    m11=(1-1E-3)*matt1(:,i);
    m12=(1-2E-3)*matt1(:,i);
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


