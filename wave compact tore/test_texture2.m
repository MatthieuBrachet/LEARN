%% // clear workspace
clear all; close all; clc;

%% // Default parameters
param.meshsize  = 128 ;     %// main grid size
param.patchsize = 200 ;     
param.windSpeed = 100  ;    %// what unit ? [m/s] ??
param.winddir   = 90   ;    %// Azimuth
param.rng = 13 ;            %// setting seed for random numbers
param.A         = 1e-7 ;    %// Scaling factor
param.g         = 9.81 ;    %// gravitational constant

param.xLim = [-10 10] ;     %// domain limits X
param.yLim = [-10 10] ;     %// domain limits Y
param.zLim = [-1e-4 1e-4]*2 ;

gridSize = param.meshsize * [1 1] ;

%% // Define the grid X-Y domain
x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
[X,Y] = meshgrid(x, y);

%% // get the grid parameters which remain constants (not time dependent)
[H0, W, Grid_Sign] =  initialize_wave( param ) ;

%% // calculate wave at t0
t0 = 0 ;
Z = calc_wave( H0 , W , t0 , Grid_Sign ) ;

%% // populate the display panel
h.fig  = figure('Color','w') ;
h.ax   = handle(axes) ;                 %// create an empty axes that fills the figure
h.surf = handle( surf( NaN(2) ) ) ;     %// create an empty "surface" object

%% // Display the initial wave surface
set( h.surf , 'XData',X , 'YData',Y , 'ZData',Z )
set( h.ax   , 'XLim',param.xLim , 'YLim',param.yLim , 'ZLim',param.zLim )

%% // Change some rendering options
axis off                                %// make the axis grid and border invisible
shading interp                          %// improve shading (remove "faceted" effect)
blue = linspace(0.4, 1.0, 25).' ; cmap = [blue*0, blue*0, blue]; %'// create blue colormap
colormap(cmap)
%// configure lighting
h.light_handle = lightangle(-45,30) ;   %// add a light source
set(h.surf,'FaceLighting','phong','AmbientStrength',.3,'DiffuseStrength',.8,'SpecularStrength',.9,'SpecularExponent',25,'BackFaceLighting','unlit')

%% // Animate
view(75,55) %// no need to reset the view inside the loop ;)

timeStep = 1./25 ;
nSteps = 2000 ;
for time = (1:nSteps)*timeStep    
    %// update wave surface
    Z = calc_wave( H0,W,time,Grid_Sign ) ;
    h.surf.ZData = Z ;
    pause(0.001);
end


%% // This block of code is only if you want to generate a GIF file
%// be carefull on how many frames you put there, the size of the GIF can
%// quickly grow out of proportion ;)

nFrame = 55 ;
gifFileName = 'MyDancingWaves.gif' ;

view(-70,40)
clear im
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;
iframe = 0 ;
for time = (1:nFrame)*.5
    %// update wave surface
    Z = calc_wave( H0,W,time,Grid_Sign ) ;
    h.surf.ZData = Z ;
    pause(0.001);

    f = getframe;
    iframe= iframe+1 ;
    im(:,:,1,iframe) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,gifFileName,'DelayTime',0,'LoopCount',inf)
disp([num2str(nFrame) ' frames written in file: ' gifFileName])