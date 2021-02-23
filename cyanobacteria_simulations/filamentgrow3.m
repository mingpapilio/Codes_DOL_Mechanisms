%This function simulates the growth of a (trials) number cyanobacteria filaments from a spore
%to a maximum size of (L) cells. Input arguments are the mechanism stratetgy
%(q, s, d, v) and the parameters of the model. The output is
%characteristics of each filament growth such as fitness, number of
%reproductive cells, growth time, ect.

function res = filamentgrow3(start, params, deets)


%trait values
q=start(1);
signal=start(2);
thresh=start(3);
sense=start(4);

%parameter values
L=params(1);
phi=params(2);
theta=params(3);
eta1=params(4);
eta2=params(5);
cap1=params(6);
growthmax=params(7);
eta3=params(8);
noise=params(9);
lamb=params(10);
lamb1=params(11);
lamb2=params(12);
lamb3=params(13);

%number of independent simulations
trials=deets(1);
slow=0;


%growth cap of filaments increases if cells are more sensitive
cap=cap1+exp(lamb2*sense)-1;



res=zeros(trials, 10+L);
cells=linspace(1,L,L)';

%iterate through each simulation
for t=1:trials
    if mod(t,1000)==0
        disp(t)
    end
    %initialise a filament with 2 helpers and 2 reproductives
    nowRoll=zeros(L,1); %keep track of phenotypes along filament
    nextRoll=zeros(L,1);
    nowRoll(2)=1;
    nowRoll(3)=1;
    nowRoll(1)=-1;
    nowRoll(4)=-1;
    reps=[2 3]; %position of reps in filament
    heps=[1 4]; %position of helpers in filament
    
    goodflux=phi*ones(L,1);
    gradflux=zeros(L,1);
    growth=zeros(L,L);
    nowGrow=zeros(L,1);
    nextGrow=zeros(L,1);
    Gs=zeros(L-1,1);
    timelog=-1*ones(L-3,1);
    timelog(1)=0;
    fil=length(heps)+length(reps);
    if length(heps)>0
        
        for v=1:length(heps)
            k=heps(v);
            dist=abs(cells-k);
            x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(eta1.^(dist(cells(1:fil))));
            goodflux(reps)=goodflux(reps)+x*(eta1.^dist(reps)).*(eta3.^dist(reps));
        end
    end
    for n=1:(L-4)
        growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
        times=(cap-nowGrow(cells))./growthrate;
        rep=datasample(reps(find(min(times(reps))==times(reps))),1);
        temp=min(times);%max(min(stimes),minrep);
        timelog(n+1)=timelog(n)+temp;
        nextGrow(1:rep-1)=min(nowGrow(1:rep-1)+growthrate(1:rep-1)*temp,cap);
        nextGrow(rep+2:L)=min(nowGrow(rep+1:(L-1))+growthrate(rep+1:(L-1))*temp,cap);
        gradflux=zeros(L,1);
        if rand<0.5
            
            nextRoll(rep+2:L)=nowRoll(rep+1:(L-1));
            nextRoll(1:rep)=nowRoll(1:rep);
            heps=find(nextRoll==-1);
            nextRoll(rep+1)=1;
            reps=find(nextRoll==1);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(eta2.^dist);
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
            goodflux=phi*ones(L,1);
            gradflux=zeros(L,1);
            fil=length(heps)+length(reps);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(eta1.^(dist(cells(1:fil))));
                    goodflux(reps)=goodflux(reps)+x*(eta1.^dist(reps)).*(eta3.^dist(reps));
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(eta2.^dist);
                end
            end
            
            nextGrow(none)=0;
            nextGrow(heps)=0;
            nextGrow(rep)=0;
            
        else
            
            nextRoll(rep+1:L)=nowRoll(rep:(L-1));
            nextRoll(1:rep)=nowRoll(1:rep);
            nextRoll(rep)=1;
            heps=find(nextRoll==-1);
            reps=find(nextRoll==1);
            
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(eta2.^dist);
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
            
            goodflux=phi*ones(L,1);
            gradflux=zeros(L,1);
            fil=length(heps)+length(reps);
            if length(heps)>0
                for v=1:length(heps)
                    k=heps(v);
                    dist=abs(cells-k);
                    x=theta*((1-(signal-slow)/(1-slow))^lamb1)/sum(eta1.^(dist(cells(1:fil))));
                    goodflux(reps)=goodflux(reps)+x*(eta1.^dist(reps)).*(eta3.^dist(reps));
                    dist=abs(cells-k);
                    gradflux=gradflux+lamb*signal*(eta2.^dist);
                end
            end
            
            
        end
        nowRoll=nextRoll;
        nowGrow=nextGrow;
        nextGrow=zeros(L,1);
        nextRoll=zeros(L,1);
        growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
        reps=find(nowRoll==1);
        heps=find(nowRoll==-1);
        inter=heps(2:end)-heps(1:end-1);
        holder=zeros(1,L);
        holder(1:length(inter))=inter;
        res(t,:)=[timelog(end) length(reps) mean(goodflux(reps)) (1/(timelog(end)))*sum(growthrate(reps)) var(goodflux(reps)) sum(Gs.*(1-Gs)) mean(growthrate(reps)) var(growthrate(reps)) mean(inter) var(inter) nowRoll'];
        
    end
    
end
