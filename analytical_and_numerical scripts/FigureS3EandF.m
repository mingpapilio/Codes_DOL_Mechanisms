%This script generate figure S3E and S3F. We create a map of parameter space (group size and essentiality
%of cooperation) and evaluate which mechanism (fully coordinated specialisation
%or random specialisation) is favoured. We assume here that groups grow, that division of labour
%is required each generation of the group growth and
%that a cost of coordination is paid each generation of the group growth.

clear all

Nmin=2;
Nmax=250; %Maximmum size that groups grow to



HBmax=50; %Maximum importance of cooperation

n=linspace(Nmin,Nmax,(Nmax-Nmin)/4+1)
HBr=linspace(1,HBmax,length(n))
HBl=fliplr(1./HBr(2:end))
hb=[HBl HBr] 
wcwr=[];
trials=100; %Number of independent simulations of group growth to estimate fitness of each strategy


track1=[];
track2=[];


%iterate through parameter space
for x=1:length(n)
    disp(x);
    for y=1:length(hb)
        N=n(x);
        HB=hb(y);
        cells=linspace(1,N,N);

        optp=max(0,(HB-1)/(2*HB));
        
        for k=1:trials
            roles=[];
            fit=[];
            roles(end+1)=1;
            
            
            for i=1:(N-1)
                heps=sum(roles==-1);
                port=heps/i;
                fit(end+1)=(1-port)*(1+HB*port);
                if rand<optp
                    roles(end+1)=-1;
                else
                    roles(end+1)=1;
                end
            end
            heps=sum(roles==-1);
            port=heps/(N);
            fit(end+1)=(1-port)*(1+HB*port);
            
            track2(k+1,1)=prod(fit);
            track2(k+1,2)=port;
        end
        
        for k=1:trials
            roles=[];
            fit=[];
            roles(end+1)=1;
            optp=max(0,(HB-1)/(2*HB));
            
            
            for i=1:(N-1)
                heps=sum(roles==-1);
                port=heps/i;
                fit(end+1)=(1-port)*(1+HB*port);
                next1=(heps+1)/(i+1);
                next2=(heps)/(i+1);
                dist1=abs(next1-optp);
                dist2=abs(next2-optp);
                if dist1<dist2
                    roles(end+1)=-1;
                elseif dist2<dist1
                    roles(end+1)=1;
                else
                    if rand<0.5
                        roles(end+1)=-1;
                    else
                        roles(end+1)=1;
                    end
                end
                
            end
            heps=sum(roles==-1);
            port=heps/(N);
            fit(end+1)=(1-port)*(1+HB*port);
            
            track1(k+1,1)=prod(fit);
            track1(k+1,2)=port;
            
            
        end
        
        wc(x,y)=mean(track1(:,1));
        wr(x,y)=mean(track2(:,1));
    end
end



ftsz=50;
cc2=0.015;
cc1=0.03;
cost1=ones(size(wc,1), size(wc,2));
cost2=ones(size(wc,1), size(wc,2));
for i=1:length(n)
    for j=1:length(hb)
        cost1(i,j)=(1-cc1)^(n(i)-1);
        cost2(i,j)=(1-cc2)^(n(i)-1);
    end
end

%Plot results for high cost of coordination
figure('Position', [10 10 800 800])
cmax=10
imagesc(log(hb), n, (cost1.*wc>wr))
set(gca,'YDir','normal')
col=[gray]
colormap(1-col)
set(gca,'fontsize', ftsz)
ax = gca;
ax.FontSize = ftsz/2; 
ax.XTick=[log(1/HBmax) 0 log(HBmax)];
ax.XTickLabels={"1/"+HBmax "1" HBmax};
ax.YTick=[Nmin Nmax/2 Nmax];
box off
  set(gcf,'color','w');  
  hold on
x=log(linspace(min(hb)-0.01,1));
curve2=Nmax*ones(1,length(x));
curve1=-2*ones(1,length(x));
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);
fill(x2, inBetween, 'r ');
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);
line([0 0], [-2, Nmax], "Color",'r');

%Plot results for low cost of coordination
figure('Position', [10 10 800 800])
cmax=10
imagesc(log(hb), n, (cost2.*wc>wr))
set(gca,'YDir','normal')
col=[gray]
colormap(1-col)
set(gca,'fontsize', ftsz)
ax = gca;
ax.XTick=[log(1/HBmax) 0 log(HBmax)];
ax.XTickLabels={"1/"+HBmax "1" HBmax};
ax.YTick=[Nmin Nmax/2 Nmax];
box off
  set(gcf,'color','w'); 
  ax.FontSize = ftsz/2; 
  hold on
x=log(linspace(min(hb)-0.01,1));
curve2=Nmax*ones(1,length(x));
curve1=-2*ones(1,length(x));
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);
fill(x2, inBetween, 'r ');
plot(x, curve1, 'r', 'LineWidth', 1);
plot(x, curve2, 'r', 'LineWidth', 1);
line([0 0], [-2, Nmax], "Color",'r');


