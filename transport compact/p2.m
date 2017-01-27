clc; clear all; close all;
global n
global X Y

n=100;
mod2;

%% velocity data
cx=0.1*ones(n+1,n+1);
cy=0.16*ones(n+1,n+1);

%% time data
Tmax=100;
cfl=0.9;
c=max(max([cx cy]));
ddt=cfl*h/c;

%% initial function
r=sqrt((X-.5).^2+(Y-.5).^2);
%u=.5*(1+cos(pi*r/.25)).*(r<.25);
%u=exp(-100.*r.^2);
u=(r<0.25);

%% iterations
t=0;
while t<Tmax
    t=t+ddt;
    clc; disp([t Tmax])
    
    % K1
    w=u;
    [gradx,grady] = grad2(w);
    k1=-cx.*gradx-cy.*grady;
    
    % K2
    w=u+ddt/2*k1;
    [gradx,grady] = grad2(w);
    k2=-cx.*gradx-cy.*grady;
    
    % K3
    w=u+ddt/2*k2;
    [gradx,grady] = grad2(w);
    k3=-cx.*gradx-cy.*grady;
    
    % K4
    w=u+ddt*k3;
    [gradx,grady] = grad2(w);
    k4=-cx.*gradx-cy.*grady;
    
    % Assemblage
    ut=u+ddt/6*(k1+2*k2+2*k3+k4);
    
    % filtrage
    ut=reshape(ut,[],1);
    u=.5*(Fx*Fy*ut+Fy*Fx*ut);
    u=reshape(u,n+1,n+1);
    
    
    % figure
    pause(0.001)
    figure(1)
    contourf(X,Y,u)
    title(['solution at time : ' num2str(t)])
    colorbar
    
end