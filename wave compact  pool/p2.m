clc; clear all; close all;
% solve the wave equation on a swimming pool with hardly viscous term on each equation.
global n dx
global X Y
global h0 gp mu cgrav length 

n=255;
mod2;
video='yes';

%% time data
Tmax=.6;
cfl=.9;
ddt=cfl*dx/cgrav;

%% initial function
test = 2;
if test == 1
    r1=sqrt((X-.15*length).^2+(Y-.1*length).^2);
    h1=h0.*exp(-(r1.^2./(.01*length)));
    r2=sqrt((X-.8*length).^2+(Y-.5*length).^2);
    h2=.75*h0.*exp(-(r2.^2./(.01*length)));
    h=h1+h2;
elseif test == 2
    r=sqrt((X-.3*length).^2+(Y-.5*length).^2);
    h=h0.*exp(-(r.^2./(.01*length)));
end
ux=zeros(size(X));
uy=zeros(size(X));

mat(1)=1;
energy(1)=1;
time(1)=0;
ref_mat=sum(sum(h));
ref_energy=gp*sum(sum(h.^2))+h0.*sum(sum(abs(ux.^2+uy.^2)));

%% iterations
if strcmp(video,'yes')==1
    vidObj=VideoWriter(['video_wave_test' test '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

t=0;iter=1;
while t<Tmax
    iter=iter+1;
    t=t+ddt;
    clc; disp([t Tmax mat(end) energy(end)])
    
    %% K1
    uux=ux;
    uuy=uy;
    hh=h;
    [div] = div2(uux,uuy);
    [visc] = visc2(hh);
    kh1=-h0*div-mu*visc;
    [gradx,grady] = grad2(hh);
    [viscvx] = visc2(uux);
    [viscvy] = visc2(uuy);
    ku1x=-gp.*gradx-mu*viscvx;
    ku1y=-gp.*grady-mu.*viscvy;
    
    %% K2
    uux=ux+ddt/2*ku1x;
    uuy=uy+ddt/2*ku1y;
    hh=h+ddt/2*kh1;
    [div] = div2(uux,uuy);
    [visc] = visc2(hh);
    kh2=-h0*div-mu*visc;
    [gradx,grady] = grad2(hh);
    [viscvx] = visc2(uux);
    [viscvy] = visc2(uuy);
    ku2x=-gp.*gradx-mu*viscvx;
    ku2y=-gp.*grady-mu.*viscvy;
    
    %% K3
    uux=ux+ddt/2*ku2x;
    uuy=uy+ddt/2*ku2y;
    hh=h+ddt/2*kh2;
    [div] = div2(uux,uuy);
    [visc] = visc2(hh);
    kh3=-h0*div-mu*visc;
    [gradx,grady] = grad2(hh);
    [viscvx] = visc2(uux);
    [viscvy] = visc2(uuy);
    ku3x=-gp.*gradx-mu*viscvx;
    ku3y=-gp.*grady-mu.*viscvy;
    
    %% K4
    uux=ux+ddt*ku3x;
    uuy=uy+ddt*ku3y;
    hh=h+ddt*kh3;
    [div] = div2(uux,uuy);
    [visc] = visc2(hh);
    kh4=-h0*div-mu*visc;
    [gradx,grady] = grad2(hh);
    [viscvx] = visc2(uux);
    [viscvy] = visc2(uuy);
    ku4x=-gp.*gradx-mu*viscvx;
    ku4y=-gp.*grady-mu.*viscvy;
    
    %% Assemblage
    ux=ux+ddt/6*(ku1x+2*ku2x+2*ku3x+ku4x);
    uy=uy+ddt/6*(ku1y+2*ku2y+2*ku3y+ku4y);
    h=h+ddt/6*(kh1+2*kh2+2*kh3+kh4);
    

    
    %% figure
    if strcmp(video,'yes')==1
        pause(10e-10)
        figure(1)
        surf(X,Y,h)
        shading interp;
        axis([0 length 0 length -h0/3 h0])
        title(['solution at time : ' num2str(t)])
        colorbar
        caxis([-h0/3 h0*2/3])
        
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
    mat(iter)=sum(sum(h))./ref_mat;
    energy(iter)=(gp*sum(sum(h.^2))+h0.*sum(sum(abs(ux.^2+uy.^2))))./ref_energy;
    time(iter)=t;
end
if strcmp(video,'yes') == 1
    close(vidObj);
end


figure(1)
subplot(121)
plot(time,mat-1)
title('relative mass')
xlabel('time (sec.)')

subplot(122)
plot(time,energy-1)
title('relative energy')
xlabel('time (sec.)')

figure(2)
surf(X,Y,h)
shading interp;
axis([0 length 0 length -h0/3 h0])
title(['solution at time : ' num2str(t)])