clc; clear all; close all;
% solve the Shallow Water on a periodic square with Coriolis (beta plan
% approximation) effect.
global test
global n dx
global X Y
global h0 gp cgrav length ccor

n=63;
mod2;
video='yes';

%% time data
Tmax=5;
cfl=.7;
ddt=cfl*dx/max([cgrav,ccor]);

%% initial function
test = 2;
if test == 1
    r=sqrt((X-.3*length).^2+(Y-.5*length).^2);
    h=h0.*exp(-(r.^2./(.01*length)))+h0;
elseif test == 2
    r1=sqrt((X-.15*length).^2+(Y-.1*length).^2);
    h1=h0.*exp(-(r1.^2./(.01*length)))+h0;
    r2=sqrt((X-.8*length).^2+(Y-.5*length).^2);
    h2=.75*h0.*exp(-(r2.^2./(.01*length)));
    h=h1+h2;
end
ux=zeros(size(X));
uy=zeros(size(X));

mat(1)=1;
energy(1)=1;
time(1)=0;
ref_mat=sum(sum(h));

%% iterations
if strcmp(video,'yes')==1
%     vidObj=VideoWriter(['video_swe.avi']);
%     open(vidObj);
%     set(gca,'nextplot','replacechildren');
%    rmdir('./video/' date '/')
    mkdir(['./swe-video-' date ])
end

t=0;iter=1;
while t<Tmax
    iter=iter+1;
    t=t+ddt;
    clc; disp([t Tmax mat(end)])
    
    %% K1
    tt=t;
    uux=ux;
    uuy=uy;
    hh=h;
    
    hs=relief(X,Y,tt);
    % eq 1
    [div] = div2((hh-hs).*uux,(hh-hs).*uuy);
    kh1=-div;
    
    % eq 2
    gr=0.5*(uux.^2+uuy.^2)+gp*hh;
    [gradx,grady] = grad2(gr);
    [ vort ] = vort2( uux,uuy );
    ku1x=-gradx-vort.*uuy;
    ku1y=-grady+vort.*uux;
    
    %% K2
    tt=t+ddt/2;
    uux=ux+ddt/2*ku1x;
    uuy=uy+ddt/2*ku1y;
    hh=h+ddt/2*kh1;
    
    hs=relief(X,Y,tt);
    % eq 1
    [div] = div2((hh-hs).*uux,(hh-hs).*uuy);
    kh2=-div;
    
    % eq 2
    gr=0.5*(uux.^2+uuy.^2)+gp*hh;
    [gradx,grady] = grad2(gr);
    [ vort ] = vort2( uux,uuy );
    ku2x=-gradx-vort.*uuy;
    ku2y=-grady+vort.*uux;
    
    %% K3
    tt=t+ddt/2;
    uux=ux+ddt/2*ku2x;
    uuy=uy+ddt/2*ku2y;
    hh=h+ddt/2*kh2;
    
    hs=relief(X,Y,tt);
    % eq 1
    [div] = div2((hh-hs).*uux,(hh-hs).*uuy);
    kh3=-div;
    
    % eq 2
    gr=0.5*(uux.^2+uuy.^2)+gp*hh;
    [gradx,grady] = grad2(gr);
    [ vort ] = vort2( uux,uuy );
    ku3x=-gradx-vort.*uuy;
    ku3y=-grady+vort.*uux;
    
    %% K4
    tt=t+ddt;
    uux=ux+ddt*ku3x;
    uuy=uy+ddt*ku3y;
    hh=h+ddt*kh3;
    
    hs=relief(X,Y,tt);
    % eq 1
    [div] = div2((hh-hs).*uux,(hh-hs).*uuy);
    kh4=-div;
    
    % eq 2
    gr=0.5*(uux.^2+uuy.^2)+gp*hh;
    [gradx,grady] = grad2(gr);
    [ vort ] = vort2( uux,uuy );
    ku4x=-gradx-vort.*uuy;
    ku4y=-grady+vort.*uux;
    
    
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
        clf
        
        figure(1)
        hs=relief(X,Y,t);
        plot2curve(h,hs)
        title(['solution at time : ' num2str(t)])
        
        % caxis([-h0/3 1*h0])
        % colormap winter
        
%         currFrame = getframe;
%         writeVideo(vidObj,currFrame);
        print('-dpng', ['./swe-video-' date '/swe' num2str(1000000000+iter) '.png'])
    end
    
    mat(iter)=sum(sum(h))./ref_mat;
    time(iter)=t;
end
% if strcmp(video,'yes') == 1
%     close(vidObj);
% end


figure(1)
plot(time,mat-1)
title('relative mass')
xlabel('time (sec.)')


figure(2)
surf(X,Y,h)
shading interp;
title(['solution at time : ' num2str(t)])