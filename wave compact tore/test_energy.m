% test IPP
clc; clear all; close all;

n=1023;
Q=diag(ones(n,1),1)-diag(ones(n,1),-1);
Q(1,end)=-1;
Q(end,1)=1;
Q=sparse(Q);

P=4/6*diag(ones(n+1,1))+1/6*(diag(ones(n,1),1)+diag(ones(n,1),-1));
P(1,end)=1/6;
P(end,1)=1/6;
P=sparse(P);

B=P\Q;

for ite=1:1000
    clc; ite
    
    u=rand(n+1,1);
    h=rand(n+1,1);

    s=B*(h.*u)-h.*(B*u)-u.*(B*h);
    maxi(ite)=max(s);
    somme(ite)=sum(s);
    
end

figure(1)
plot(somme)
title('err. dans IPP form.')
