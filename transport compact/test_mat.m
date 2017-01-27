clc;

n=10;
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

id=speye(size(Q));

D=P\Q;
D1=kron(id,D);

P2=kron(id,P);
Q2=kron(id,Q);
D2=P2\Q2;

surf(D2-D1)