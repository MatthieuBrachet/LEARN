clc; clear all; close all;

n=100;
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
F=sparse(ftr);



%% initial function and velocity
R=0.1;
r=sqrt((x-.5).^2);
u=.5*(1+cos(pi*r/R)).*(r<R);

c=.03;

plot(x,u)

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
    w=Q*u;
    k1=-c*P\w;
    
    % K2
    w=Q*(u+ddt/2*k1);
    k2=-c*P\w;
    
    % K3
    w=Q*(u+ddt/2*k2);
    k3=-c*P\w;
    
    % K4
    w=Q*(u+ddt*k3);
    k4=-c*P\w;
    
    % ASSEMBLAGE
    ut=u+ddt/6*(k1+2*k2+2*k3+k4);
    
    % FILTRAGE
    u=F*ut;
    
    pause(0.01)
    figure(1)
    plot(x,u)
    axis([0 1 0 1])
    
end
    
    
    
    
