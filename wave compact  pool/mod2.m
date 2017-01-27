% MOD WAVE EQUATION 2D
global n dx
global X Y
global Qx Qy Px Py lap
global h0 radius gp length mu
global cgrav

h0=.9;
radius=6.37122d+06;
gp=9.80616;
cgrav=sqrt(gp*h0);
length=1;
mu=10;

dx=length./(n+1);
x=dx:dx:length-dx;
[X,Y]=meshgrid(x,x);

Q=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
Q=3/2*Q;
Q(1,1)=-2;
Q(1,2)=3;
Q(1,3)=-2/3;
Q(1,4)=1/8;

Q(end,end)=2;
Q(end,end-1)=-3;
Q(end,end-2)=2/3;
Q(end,end-3)=-1/8;
Q=Q./(2*dx);

P=1/4*diag(ones(n-1,1),1)+1/4*diag(ones(n-1,1),-1)+speye(n,n);

id=speye(n,n);

Px=kron(P,id);
Py=kron(id,P);
Qx=kron(Q,id);
Qy=kron(id,Q);


 N=2*sparse(eye(n,n))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
N=(6/5)*N;
N(1,1)=2681/480;
N(1,2)=-23/3;
N(1,3)=113/40;
N(1,4)=-13/15;
N(1,5)=59/480;
%N(1,6)=33/40;
N(end,end)=N(1,1);
N(end,end-1)=N(1,2);
N(end,end-2)=N(1,3);
N(end,end-3)=N(1,4);
N(end,end-4)=N(1,5);
%N(end,end-5)=N(1,6);
M=sparse(eye(n,n))+1/10*diag(ones(n-1,1),1)+1/10*diag(ones(n-1,1),-1);

lap=M\N;
lap=kron(lap,id)+kron(id,lap);