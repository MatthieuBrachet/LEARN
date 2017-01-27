clc; clear all; close all;

gp=9.81;
hp=10;

n=2047;
dx=1./(n+1);
x=[dx:dx:1]';

cgrav=sqrt(gp*hp);
cfl=.9;
ddt=cfl*dx/cgrav;

t=0;

Q=diag(ones(n,1),1)-diag(ones(n,1),-1);
Q(1,end)=-1;
Q(end,1)=1;
Q=sparse(Q)/(2*dx);

P=4/6*diag(ones(n+1,1))+1/6*(diag(ones(n,1),1)+diag(ones(n,1),-1));
P(1,end)=1/6;
P(end,1)=1/6;
P=sparse(P);

B=P\Q;

r=sqrt((x-.5).^2);
h=exp(-100*(r.^2));
u=zeros(size(x));

E=[];M=[];
while t<1
    clc; t
    
    t=t+ddt;
    
    %% K1
    K1h=-hp*B*u;
    K1v=-gp*B*h;
    
    %% K2
    K2h=-hp*B*(u+ddt/2*K1v);
    K2v=-gp*B*(h+ddt/2*K1h);
    
    %% K3
    K3h=-hp*B*(u+ddt/2*K2v);
    K3v=-gp*B*(h+ddt/2*K2h);
    
    %% K4
    K4h=-hp*B*(u+ddt*K3v);
    K4v=-gp*B*(h+ddt*K3h);
    
    %% assemblage
    h=h+ddt/6*(K1h+2*K2h+2*K3h+K4h);
    u=u+ddt/6*(K1v+2*K2v+2*K3v+K4v);
    
%     pause(0.0001)
%     clf
%     figure(1)
%     plot(x,h)
%     axis([0 1 0 1])
    
    nrj=sum(gp.*h.^2+hp.*u.^2);
    mass=sum(h);
    E=[E nrj];
    M=[M mass];
    

end

figure(2)
plot(E./E(end)-1); hold on
plot(M./M(end)-1); hold off
legend('energy','mass')

