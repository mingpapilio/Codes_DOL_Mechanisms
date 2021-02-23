%Script to run the simulation of a growing cyanobacteria filament, see
%filamentgrow1.m for more details. The output of this script is the
%phenotyps of individual cells at the end of the filament growth. 

function res = filamentgrow2(start, params)

q=start(1);
signal=start(2);
thresh=start(3);
sense=start(4);

N=params(1);
phi=params(2);
theta=params(3);
gamma1=params(4);
gamma2=params(5);
cap1=params(6);
growthmax=params(7);
gamma3=params(8);
noise=params(9);
lamb=params(10);
lamb1=params(11);
lamb2=params(12);
lamb3=params(13);

slow=0;


cap=cap1+exp(lamb2*sense)-1;




sz=8000/N;
offset=0.15;
repcol=[0.5 0.5 1];
hepcol=[0.5 1 0.5];
bord=[0.33 0.33 0.33];

res=zeros(N,N-3);
cells=linspace(1,N,N)';
roles=zeros(N,N-3);
    nowRoll=zeros(N,1);
    nextRoll=zeros(N,1);
    nowRoll(2)=1;
    nowRoll(3)=1;
    nowRoll(1)=-1;
    nowRoll(4)=-1;
roles(:,1)=nowRoll;
    reps=[2 3];
    heps=[1 4];
    goodflux=phi*ones(N,1);
    gradflux=zeros(N,1);
    growth=zeros(N,N);
    nowGrow=zeros(N,1);
    nextGrow=zeros(N,1);
    Gs=zeros(N-1,1);
    timelog=-1*ones(N-3,1);
    timelog(1)=0;
    fil=length(heps)+length(reps);
    if length(heps)>0
    for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
                    goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
    end
    end
  
        for n=1:(N-4)
        growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
        times=(cap-nowGrow(cells))./growthrate;
        rep=datasample(reps(find(min(times(reps))==times(reps))),1);
        temp=min(times);%max(min(stimes),minrep);
        timelog(n+1)=timelog(n)+temp;
        nextGrow(1:rep-1)=min(nowGrow(1:rep-1)+growthrate(1:rep-1)*temp,cap);
        nextGrow(rep+2:N)=min(nowGrow(rep+1:(N-1))+growthrate(rep+1:(N-1))*temp,cap);
        gradflux=zeros(N,1);
        if rand<0.5
            
            nextRoll(rep+2:N)=nowRoll(rep+1:(N-1));
            nextRoll(1:rep)=nowRoll(1:rep);
            heps=find(nextRoll==-1);
            nextRoll(rep+1)=1;
            reps=find(nextRoll==1);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(gamma2.^dist);
                end
            end
            receive=max(0,normrnd(gradflux(rep+1),noise,1));
            G=1-2./(1+exp(-(sense).*(receive-thresh)));
            flip=min(max(q+G,0),1);
            Gs(n)=flip;
            if rand<flip
                nextRoll(rep+1)=-1;
            else
                nextRoll(rep+1)=1;
            end
            
         
            reps=find(nextRoll==1);
            heps=find(nextRoll==-1);
            none=find(nextRoll==0);
            goodflux=phi*ones(N,1);
            gradflux=zeros(N,1);
            fil=length(heps)+length(reps);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
                    goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(gamma2.^dist);
                end
            end
                      
            nextGrow(none)=0;
            nextGrow(heps)=0;
            nextGrow(rep)=0;
            
        else
            
            nextRoll(rep+1:N)=nowRoll(rep:(N-1));
            nextRoll(1:rep)=nowRoll(1:rep);
            nextRoll(rep)=1;
            heps=find(nextRoll==-1);
            reps=find(nextRoll==1);

            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(gamma2.^dist);
                end
            end
            
           
            receive=max(0,normrnd(gradflux(rep),noise,1));
            G=1-2./(1+exp(-(sense).*(receive-thresh)));
            flip=min(max(q+G,0),1);
            Gs(n)=flip;
            if rand<flip
                nextRoll(rep)=-1;
            else
                nextRoll(rep)=1;
            end
            reps=find(nextRoll==1);
            heps=find(nextRoll==-1);
            none=find(nextRoll==0);
            nextGrow(none)=0;
            nextGrow(heps)=0;
            
            goodflux=phi*ones(N,1);
            gradflux=zeros(N,1);
            fil=length(heps)+length(reps);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
                    goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(gamma2.^dist);
                end
            end
            

        end
        nowRoll=nextRoll;
        nowGrow=nextGrow;
        nextGrow=zeros(N,1);
        nextRoll=zeros(N,1);
        growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));


    
        roles(:,n+1)=nowRoll;

        end
            res=roles;
    
  