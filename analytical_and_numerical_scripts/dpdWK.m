%Returns the derivative of the WK function with respect to the probability
%of becoming a helper. This is needed to calculate the optimal probability
%of adopting a helper role for random specialisers. 

function valuef = dpdWK(p, params)
FK=params;
N=length(FK)-1;
value=[];
for k=0:N
value(end+1,:)=FK(k+1).*(nchoosek(N,k)).*(k.*(p.^(k-1)).*((1-p).^(N-k))-(N-k).*(p.^(k)).*((1-p).^(N-k-1)));
end
valuef=sum(value);