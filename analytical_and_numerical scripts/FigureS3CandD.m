%This script generate figure S3C or D. We create a map of parameter space (group size and essentiality
%of cooperation) and evaluate which mechanism (fully coordinated specialisation
%or random specialisation) is favoured. We assume here that the degree of helper
%specialisation, x, may vary.

clear all

%Parameter space

Nmax=1000;
Nmin=2;
HBmax=50;
HBmin=0;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1); %discretisation for group size
HBr=linspace(1,HBmax,length(size))
HBl=fliplr(1./HBr(2:end))
HB=[HBl HBr] %discretisation for relative importance of cooperation 
ftsz=50;
%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers. Evaluate the fitness of a coordinated
%group, and random groups. Store in results. 

x=0.75;  %degree of helper specialisation
gamma=2;  %efficiency benefits from specialisation
%determine relative fitness
matt=[];

results=[]
coordbenss=zeros(length(size), length(HB));
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
        popt=max(0,(hb*(x^(gamma-1))-1)/(2*hb*(x^gamma)));
        matt4(i,j)=popt;
        if hb*(x^(gamma-1))>1
            coordbens(i,j)=(-(hb*(x^(gamma-1))-1)^2+(hb*(x^(gamma-1))-1)*(2*hb*(x^gamma)))/(N*((hb*(x^(gamma-1))+1)^2));
        else
            coordbens(i,j)=0;
        end
    end
end



for i=1:length(HB)
    c1=6E-4;
    c2=1E-3;
    m=coordbens(:,i);
    popt=matt4(:,i);
    if length(find(m>c1))>0
        track1(i)=max(size(find(m>c1)));
    else
        track1(i)=0;
    end
    if length(find(m>c2))>0
        track2(i)=max(size(find(m>c2)));
    else
        track2(i)=0;
    end
    if length(find(popt<=0))>0
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


