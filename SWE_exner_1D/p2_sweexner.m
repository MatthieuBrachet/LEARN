clc; clear all; close all;

global n x dx cgrav
global h0 length gp


n=1028;
mod2;

%% *** options
opt_film='no';

%% *** initial time
t=0;
cfl=.9;
ddt=cfl*dx/cgrav;
tmax=500;

%% *** initial function
dataR=0.6;
r=sqrt((x-.2*length).^2);
h=h0+.3*h0*(1+cos(pi*r/dataR)).*(r<dataR);
r=sqrt((x-.7*length).^2);
h=h+.1*h0*(1+cos(pi*r/dataR)).*(r<dataR);

dataR=0.1*length;
r=sqrt((x-0.5*length).^2);
hs=0.05*h0*(1+cos(pi*r/dataR)).*(r<dataR)+0.01*rand(size(x))+0.25;

hstar=h-hs;
consref=sum(hstar);

u=zeros(size(x));

%% *** iterations
cons=[];
maxi=[];
while t<tmax
    clc; t=t+ddt;
    [sum(hstar)./consref t]
    
    % K1
    hhstar=hstar;
    hhs=hs;
    uu=u;
    
    K1hstar=-P\(Q*(hhstar.*uu));
    K1u=-P\(Q*(.5*uu.^2+gp.*(hhstar+hhs)));
    [ Qu ] = exnerlaw( uu, hhstar, t );
    K1hs=-P\(Q*Qu);
    
    % K2
    hhstar=hstar+ddt/2*K1hstar;
    hhs=hs+ddt/2*K1hs;
    uu=u+ddt/2*K1u;
    
    K2hstar=-P\(Q*(hhstar.*uu));
    K2u=-P\(Q*(.5*uu.^2+gp.*(hhstar+hhs)));
    [ Qu ] = exnerlaw( uu, hhstar, t+ddt/2 );
    K2hs=-P\(Q*Qu);
    
    % K3
    hhstar=hstar+ddt/2*K2hstar;
    hhs=hs+ddt/2*K2hs;
    uu=u+ddt/2*K2u;
    
    K3hstar=-P\(Q*(hhstar.*uu));
    K3u=-P\(Q*(.5*uu.^2+gp.*(hhstar+hhs)));
    [ Qu ] = exnerlaw( uu, hhstar, t+ddt/2 );
    K3hs=-P\(Q*Qu);
    
    % K4
    hhstar=hstar+ddt*K3hstar;
    hhs=hs+ddt*K3hs;
    uu=u+ddt*K3u;
    
    K4hstar=-P\(Q*(hhstar.*uu));
    K4u=-P\(Q*(.5*uu.^2+gp.*(hhstar+hhs)));
    [ Qu ] = exnerlaw( uu, hhstar, t+ddt );
    K4hs=-P\(Q*Qu);
    
    % Assemblage
    hstar=hstar+ddt/6*(K1hstar+2*K2hstar+2*K3hstar+K4hstar);
    hs=hs+ddt/6*(K1hs+2*K2hs+2*K3hs+K4hs);
    u=u+ddt/6*(K1u+2*K2u+2*K3u+K4u);
    
    % filtrage
    hstar=F*hstar;
    hs=F*hs;
    u=F*u;
    
    cons=[cons sum(hstar)];
    maxi=[maxi max(hstar+hs)];
    
    % figure
    if strcmp(opt_film,'yes')==1
        pause(10^-16);
        figure(100)
        plot(x,hs,x,hstar+hs);
        legend('botom topography','surface water')
        axis([0 length -0.3 2.8*h0])
    end
end
    
figure(1)
plot(x,hs,x,hstar+hs);
legend('bottom topography','surface water')
axis([0 length -0.3 2*h0])
    
figure(2)
plot(x,hs)

figure(3)
plot(cons./consref-1)

figure(4)
plot(x,Qu)

figure(5)
plot(maxi)

fig_placier;


