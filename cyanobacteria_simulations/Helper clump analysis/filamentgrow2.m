%This function simulates the growth of a cyanobacteria filaments from a spore
%to a maximum size of (L) cells. Input arguments are the mechanism stratetgy
%(q, s, d, v) and the parameters of the model. The output is
%characteristics of the filament growth such as fitness, number of
%reproductive cells, growth time, and number of helper clunmps of each possible size, ect.

function [grow repstrand hepstrand fit roles time replog grolog] = filamentgrow2(start, params)

%a lot of the below is the same as for the filament grow function
q=start(1);
signal=start(2);
thresh=start(3);
sense=start(4);

L=params(1);
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



cap=cap1+exp(lamb2*sense)-1;



cells=linspace(1,L,L)';

nowRoll=zeros(L,1);
nextRoll=zeros(L,1);
nowRoll(2)=1;
nowRoll(3)=1;
nowRoll(1)=-1;
nowRoll(4)=-1;
reps=[2 3];
heps=[1 4];
goodflux=phi*ones(L,1);
gradflux=zeros(L,1);
growth=zeros(L,L);
nowGrow=zeros(L,1);
nextGrow=zeros(L,1);
Gs=zeros(L-1,1);
timelog=-1*ones(L-3,1);
timelog(1)=0;
replog=-1*ones(L-3,1);
grolog=-1*ones(L-3,1);
fil=length(heps)+length(reps);
for v=1:length(heps)
    k=heps(v);
    dist=abs(cells-k);
    x=theta*((1-signal)^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
    goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
end

for n=1:(L-4)
    growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
    times=(cap-nowGrow(cells))./growthrate;
    rep=datasample(reps(find(min(times(reps))==times(reps))),1);
    temp=min(times);%max(min(stimes),minrep);
    timelog(n+1)=timelog(n)+temp;
    replog(n)=length(reps);
    grolog(n)=mean(growthrate(reps));
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
        goodflux=phi*ones(L,1);
        gradflux=zeros(L,1);
        fil=length(heps)+length(reps);
        if length(heps)>0
            for v=1:length(heps)
                k=heps(v);
                dist=abs(cells-k);
                x=theta*((1-signal)^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
                goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
                dist=abs(cells-k);
                gradflux=gradflux+lamb*signal*(gamma2.^dist);
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
        
        goodflux=phi*ones(L,1);
        gradflux=zeros(L,1);
        fil=length(heps)+length(reps);
        if length(heps)>0
            for v=1:length(heps)
                k=heps(v);
                dist=abs(cells-k);
                x=theta*((1-signal)^lamb1)/sum(gamma1.^(dist(cells(1:fil))));
                goodflux(reps)=goodflux(reps)+x*(gamma1.^dist(reps)).*(gamma3.^dist(reps));
                dist=abs(cells-k);
                gradflux=gradflux+lamb*signal*(gamma2.^dist);
            end
        end
        
        
    end
    nowRoll=nextRoll;
    nowGrow=nextGrow;
    nextGrow=zeros(L,1);
    nextRoll=zeros(L,1);
    growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
    
    stranded=zeros(length(reps),1);  %for each reproductive, store the minimum distance to a helper
    for r=1:length(reps)
        k=reps(r);
        stranded(r)=min(abs(k-heps));
    end
    
    
    clumped=zeros(length(heps),1); %create a map of the contiguous helper (to the left) that each helper is linked to
    for r=1:length(heps)
        k=heps(r);
        temp=sort(abs(k-heps));
        if length(temp)>1
            clumped(r)=temp(2);
        end
    end
    
    
    grow=growthrate(reps);
    roles=nowRoll;
    fit=(1/(timelog(end)))*sum(grow);
    
    %map of all contiguous connections between reproductives
    transfer=eye(L,L);
    transfer(1:end-1,2:end)=transfer(1:end-1,2:end)+eye(L-1,L-1);
    transfer(2:end,1:end-1)=transfer(2:end,1:end-1)+eye(L-1,L-1); %full cell connections
    for i=1:length(heps) %break connections wherever there is a helper
        transfer(heps(i),:)=zeros(L,1);
        transfer(:,heps(i))=zeros(L,1)';
    end
    
    strands=zeros(length(reps),1);
    maxstrand=L;
    %iterate through reproductives and report the size of the reproductive
    %clump it is in
    for i=1:length(reps)
        rep=reps(i);
        start=zeros(L,1);
        start(rep)=1;
        strands(i)=sum(sum((start.*(transfer)^maxstrand)>0));
    end
    %record the number of reproductive clumps of each possible size from 1 to
    %maxstrand
    repstrand=zeros(maxstrand,1);
    for i=1:maxstrand
        repstrand(i)=sum(strands==i)/i;
    end
    
    
    %map of all contiguous connections between helpers
    transfer=eye(L,L);
    transfer(1:end-1,2:end)=transfer(1:end-1,2:end)+eye(L-1,L-1);
    transfer(2:end,1:end-1)=transfer(2:end,1:end-1)+eye(L-1,L-1);
    for i=1:length(reps) %break connections wherever there is a reproductive
        transfer(reps(i),:)=zeros(L,1);
        transfer(:,reps(i))=zeros(L,1)';
    end
    %iterate through each helper and report size of the helper clump it is
    %in
    strands=zeros(length(heps),1);
    for i=1:length(heps)
        hep=heps(i);
        start=zeros(L,1);
        start(hep)=1;
        strands(i)=sum(sum((start.*(transfer)^maxstrand)>0));
    end
    
    %record the number of helper clumps of each possible size from 1 to
    %maxstrand
    hepstrand=zeros(maxstrand,1);
    for i=1:maxstrand
        hepstrand(i)=sum(strands==i)/i;
    end
    time=timelog;
    
end
    growthrate=growthmax*(1-exp(-lamb3.*goodflux(cells)));
    grolog(n+1)=mean(growthrate(reps));
    replog(n+1)=length(reps);



