clc;

n=10000;
h=1./(n+1);
x=[h:h:1]';

P=4/6*diag(ones(n+1,1))+1/6*(diag(ones(n,1),1)+diag(ones(n,1),-1));
P(1,end)=1/6;
P(end,1)=1/6;
P=sparse(P);

Q=diag(ones(n,1),1)-diag(ones(n,1),-1);
Q(1,end)=-1;
Q(end,1)=1;
Q=sparse(Q)/(2*h);

u=cos(2*pi*x);
du=-2*pi*sin(2*pi*x);

w=Q*u;
dp=P*w;
e=norm(du-dp,2)