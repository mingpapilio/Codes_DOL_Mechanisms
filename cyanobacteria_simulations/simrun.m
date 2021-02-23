%This function takes as input which parameter combination of phi and eta
%to examine. And runs a replicates number of independent evolutionary
%simulations to estimate the long term evolutionary strategy (q,s,v,d). The
%output are seperate files for each replicate, recording the trait values 
%and estimated fitness in the last take=2000 generations of the simulation.

%This function must be called in the same directory as filamentgrow.m and
%input.csv. There must be a folder called data that the output is directed
%to.

clear all
close all
input=csvread("input.csv"); %read in the parameter combination to consider

%parameter values of the model
phik=input(1); %background density of fixed N2
eta1k=input(2); %cooperation range
phis=[0 1 2 3 4 5 6 7 8 9];    %Parameter discretization of phi and eta that is considered
etas=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
phi=phis(phik);
eta1=etas(eta1k);
N=50;
theta=50;
cap1=100;
noise=0.01;
lamb=30;
lamb1=2.5;
lamb2=0.5;
lamb3=0.25;
growthmax=10;
eta3=1;
eta2=0.5;


%parameters of the evolutionary simulation
trials=4000;
take=1;
reps=200;
replicates=5;


%iterate through each independent simulation
for r=1:replicates

%set the starting trait values
q=0;
s=0;
h=0;
v=0;
track=zeros(trials-take+1, 5);

mu=1;

%only evolve q at first.
sigq=0.001;
sigs=0;
sigv=0;
sigh=0;
sig=[sigq sigs sigh sigv];
current=[q s h v];
params=[N phi theta eta1 eta2 cap1 growthmax eta3 noise lamb lamb1 lamb2 lamb3];
qtrack=zeros(1,500-250);
res=filamentgrow(current,params,reps);
fitold=mean(res(:,4));

%iterate through the generations
for t=1:trials

    %record q values of first 500 generations
    if t<500 & t>=250
        qtrack(1,t-250+1)=current(1);
    end
    
    %at generation 500 change so that only s,v and d evolve
    if t==500
        sigs=0.00001;
        sigh=0.01;
        sigv=0.05;
        sigq=0;
        sig=[sigq sigs sigh sigv];
        randres=mean(qtrack);
        q=mean(randres(:,1));
        current=[q s h v];
    end
    %if mutation occurs
    if rand<mu
        next=max(mvnrnd(current,sig),0); %invading mutant strategy
        if next(1)>1
            next(1)=1;
        end
        if next(2)>1
            next(2)=1;
        end
        if next(4)==0 || next(2)==0
            next(3)=0;
        end
        %evaluate relative fitness of mutant and resident strategies
        run=filamentgrow(next,params,reps);
        old=filamentgrow(current,params,reps);
        fitnew=mean(run(:,4));
        fitold=mean(old(:,4));
        
        %update resident strategy if mutant has higher average fitness
        if fitnew>fitold
            current=next;
            %fitold=fitnew;
            germ=mean(run(:,2));
            germvar=var(run(:,2));
            tim=mean(run(:,1));
            cover=mean(run(:,3));
            covervar=mean(run(:,5));
            deter=mean(run(:,6));     
            inter1=mean(run(:,10)./run(:,9));     
            inter2=mean(run(:,10));  
            inter3=mean(sqrt(run(:,10))./run(:,9));
        else
            %fitold=fitnew;
            germ=mean(old(:,2));
            germvar=var(old(:,2));
            tim=mean(old(:,1));
            cover=mean(old(:,3));
            covervar=mean(old(:,5));
            deter=mean(old(:,6)); 
            inter1=mean(old(:,10)./old(:,9));     
            inter2=mean(old(:,10));
            inter3=mean(sqrt(old(:,10))./old(:,9));
        end
   
    end
    %record resident strategy in last take generations
    if t>take
        track(t-take+1,:)=[current fitold];
    end
    disp([t, current, fitold])
end
%write outupt in seperate files for each independent replicate
csvwrite("data/track_"+phi+"_"+eta1+"_"+r+".csv",track)

end


%toc





