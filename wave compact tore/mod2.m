% MOD TRANSPORT 2D
global n dx
global X Y
global Qx Qy Px Py Fx Fy
global h0 radius gp omega length
global cgrav ccor

h0=100;
radius=6.37122d+06;
gp=9.80616;
cgrav=sqrt(gp*h0);
ccor=h0*omega;
omega=0*7.292d-05;
length=10;

dx=length./(n+1);
x=dx:dx:length;
[X,Y]=meshgrid(x,x);

cor=2.*omega./(2*pi/length).*Y;

P=4/6*diag(ones(n+1,1))+1/6*(diag(ones(n,1),1)+diag(ones(n,1),-1));
P(1,end)=1/6;
P(end,1)=1/6;
P=sparse(P);

Q=diag(ones(n,1),1)-diag(ones(n,1),-1);
Q(1,end)=-1;
Q(end,1)=1;
Q=sparse(Q)/(2*dx);

id=speye(n+1,n+1);

Px=kron(P,id);
Py=kron(id,P);
Qx=kron(Q,id);
Qy=kron(id,Q);


% f0=1; f1=0; f2=0; f3=0; f4=0; f5=0;
% f0=1/2;     f1=1/4;    f2=0;    f3=0;    f4=0;    f5=0;
% f0=10/16;     f1=4/16;    f2=-1/16;    f3=0;    f4=0;    f5=0;
% f0=44/64;     f1=15/64;    f2=-6/64;    f3=1/64;    f4=0;    f5=0;
% f0=186/256;     f1=56/256;    f2=-28/256;    f3=8/256;    f4=-1/256;    f5=0;
f0=772/1024; f1=210/1024; f2=-120/1024; f3=45/1024; f4=-10/1024; f5=1/1024;

lig1=[0,1, zeros(1,n-1)];
col1=[zeros(n,1);1];
sh1=toeplitz(col1,lig1);
sh1i=inv(sh1);
sh12=sh1^2;
sh1i2=sh1i^2;
sh13=sh12*sh1;
sh1i3=sh1i2*sh1i;
sh14=sh13*sh1;
sh1i4=sh1i3*sh1i;
sh15=sh14*sh1;
sh1i5=sh1i4*sh1i;
ftr=f0*eye(n+1)+f1*(sh1+sh1i)+f2*(sh12+sh1i2)+f3*(sh13+sh1i3)+f4*(sh14+sh1i4)+f5*(sh15+sh1i5);
Fx=sparse(kron(ftr,id));
Fy=sparse(kron(id,ftr));
