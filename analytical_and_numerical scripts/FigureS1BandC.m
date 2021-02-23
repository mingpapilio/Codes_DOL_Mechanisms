%This script generates Figures S1B and S1C. We create a map of parameter space 
%(group size and relative importance of cooperation) and evaluate which mechanism 
%(fully coordinated specialisation or random specialisaiton) is favoured. We
%output figures for the cases where 1) the optimal proportion of
%helpers is constrained by the possible allocations of labour (not continuous), and 2) the
%optimal probability of becoming a helper is optimised for random groups. 

clear all
close all

%Parameters
Nmax=1000;
Nmin=2;
HBmax=10;
HBmin=0
ftsz=50;
size=linspace(Nmin,Nmax,(Nmax-Nmin)/2+1); %discretisation for group size
HBr=linspace(1,HBmax,length(size));
HBl=fliplr(1./HBr(2:end));
HB=[HBl HBr]; %discretisation for relative importance of cooperation 

%Iterate through each possible parameter combination, calculate the optimal
%proportion of helpers as well as the optimal helper probability. Evaluate 
%the fitness of a fully coordinated group, and random groups for each of these 
%scenarios. Store in matt arrays. 
matt0=zeros(length(size), length(HB));
matt1=zeros(length(size), length(HB));
matt2=zeros(length(size), length(HB));
matt3=zeros(length(size), length(HB));
matt4=zeros(length(size), length(HB));
matt5=zeros(length(size), length(HB));

for i=1:length(size)
    disp(i);
    for j=1:length(HB)
        hb=HB(j);
        N=size(i);
        K=linspace(0,N,N+1);
        p=linspace(0.01,0.99,100);
        Wbroad=((1-(K./N))).*(1+hb*((K./(N))));
        Kopt=K(max(Wbroad)==Wbroad);
        popt=mean(Kopt)/N;
        qopt=max(0,(hb-(N/(N-1)))/(2*hb));
        pcopt=max(0,(hb-1)/(2*hb));
        coorddiscfit=max(Wbroad); %fitness of coordination at discrete coordination optimum
        cordcontfit=(1-pcopt).*(1+hb*pcopt); %fitness of coordination at continuous coordination optimum
        stochdiscfit=((1-popt)).*(1+hb*(popt))-hb*popt*(1-popt)/N; %fitness of random at discrete optimum
        stochoptfit=((1-qopt)).*(1+hb*(qopt))-hb*qopt*(1-qopt)/N; %fitness of random at random optimum
        matt0(i,j)=cordcontfit;
        matt1(i,j)=coorddiscfit;
        matt3(i,j)=stochoptfit;
        matt2(i,j)=stochdiscfit;
        matt4(i,j)=popt;
        matt5(i,j)=qopt;
    end
end

%Generate Figure S1B
for i=1:length(HB)
    m11=(1-1E-3)*matt1(:,i);
    m12=(1-2E-3)*matt1(:,i);
    m2=matt2(:,i);
    popt=matt5(:,i);
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


%Generate Figure S1C

for i=1:length(HB)
    m11=(1-1E-3)*matt0(:,i);
    m12=(1-2E-3)*matt0(:,i);
    qopt=matt4(:,i);
    m2=matt3(:,i);
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
    if m2>m12 & length(find(qopt<=0))>0
        nodol(i)=max(size(find(qopt<=0)));
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



