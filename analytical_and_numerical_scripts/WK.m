%function returning the expected fitness of a random group with fecundity
%per group composition given by params and with probability of being a
%helper given by p

function valuef = WK(p, params)
FK=params;
N=length(FK)-1;
value=[];
for k=0:N
value(end+1,:)=FK(k+1).*(nchoosek(N,k)).*((p.^(k)).*((1-p).^(N-k)));
end
valuef=sum(value);