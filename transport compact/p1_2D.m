clc; clear all; close all;

n=100;
h=1./(n+1);
x=h:h:1;
[X,Y]=meshgrid(x,x);

P=4/6*diag(ones(n+1,1))+1/6*(diag(ones(n,1),1)+diag(ones(n,1),-1));
P(1,end)=1/6;
P(end,1)=1/6;
P=sparse(P);
id=speye(size(X));

Q=diag(ones(n,1),1)-diag(ones(n,1),-1);
Q(1,end)=-1;
Q(end,1)=1;
Q=sparse(Q)/(2*h);

D=P\Q;
Dx=kron(D,id);


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



%% initial function and velocity
R=0.1;
r=sqrt((X-.5).^2+(Y-.5).^2);
u=.5*(1+cos(pi*r/R)).*(r<R);
u=reshape(u,[],1);

cx=.03;
cy=0;
c=max(abs(cx),abs(cy));


%% time data
cfl=.001;
ddt=cfl*h/c;
Tmax=100;

%% iterations
t=0;
while t<Tmax
    
    t=t+ddt;
    clc; disp([t Tmax]);
    
    % K1
    w1=u;
    k1=-cx*Dx*w1;
    
    % K2
    w2=u+ddt/2*k1;
    k2=-cx*Dx*w2;
    
    % K3
    w3=u+ddt/2*k2;
    k3=-cx*Dx*w3;
    
    % K4
    w4=u+ddt*k3;
    k4=-cx*Dx*w4;
    
    % assemblage
    ut=u+ddt/6*(k1+2*k2+2*k3+k4);
    
    % filtrage
    u=Fx*u;
    
    
    % figure
    pause(0.01)
    im_u=reshape(u,size(X));
    figure(1)
    contourf(X,Y,im_u)
    colorbar
    axis([0 1 0 1])
    
end
    
    
    
    
