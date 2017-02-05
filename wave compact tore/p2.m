clc; clear all; close all;
% solve the wave equation on a periodic square
global n dx
global X Y
global h0 gp cgrav length

n=511;
mod2;
video='no';

%% time data
Tmax=1;
cfl=1;
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
    r=sqrt((X-.5*length).^2+(Y-.5*length).^2);
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
    vidObj=VideoWriter(['video_wave.avi']);
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
    kh1=-h0*div;
    [gradx,grady] = grad2(hh);
    ku1x=-gp.*gradx-cor.*uuy;
    ku1y=-gp.*grady+cor.*uux;
    
    %% K2
    uux=ux+ddt/2*ku1x;
    uuy=uy+ddt/2*ku1y;
    hh=h+ddt/2*kh1;
    [div] = div2(uux,uuy);
    kh2=-h0*div;
    [gradx,grady] = grad2(hh);
    ku2x=-gp.*gradx-cor.*uuy;
    ku2y=-gp.*grady+cor.*uux;
    
    %% K3
    uux=ux+ddt/2*ku2x;
    uuy=uy+ddt/2*ku2y;
    hh=h+ddt/2*kh2;
    [div] = div2(uux,uuy);
    kh3=-h0*div;
    [gradx,grady] = grad2(hh);
    ku3x=-gp.*gradx-cor.*uuy;
    ku3y=-gp.*grady+cor.*uux;
    
    %% K4
    uux=ux+ddt*ku3x;
    uuy=uy+ddt*ku3y;
    hh=h+ddt*kh3;
    [div] = div2(uux,uuy);
    kh4=-h0*div;
    [gradx,grady] = grad2(hh);
    ku4x=-gp.*gradx-cor.*uuy;
    ku4y=-gp.*grady+cor.*uux;
    
    %% Assemblage
    utx=ux+ddt/6*(ku1x+2*ku2x+2*ku3x+ku4x);
    uty=uy+ddt/6*(ku1y+2*ku2y+2*ku3y+ku4y);
    ht=h+ddt/6*(kh1+2*kh2+2*kh3+kh4);
    
    %% filtrage
    utx=reshape(utx,[],1);
    uty=reshape(uty,[],1);
    ht=reshape(ht,[],1);
    ux=.5*(Fx*Fy*utx+Fy*Fx*utx);
    uy=.5*(Fx*Fy*uty+Fy*Fx*uty);
    h=.5*(Fx*Fy*ht+Fy*Fx*ht);
    ux=reshape(ux,n+1,n+1);
    uy=reshape(uy,n+1,n+1);
    h=reshape(h,n+1,n+1);
    
    %% figure
    if strcmp(video,'yes')==1
        figure(1)
        surfl(X,Y,h)
        shading interp;
        alpha 0.6
        axis([0 length 0 length -h0/3 h0])
        title(['solution at time : ' num2str(t)])
        view([1 1 1.5])
        % caxis([-h0/3 1*h0])
        % colormap winter
        
        
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