clc; clear all; close all;
global n
global X Y

n=10;
mod2;

u=cos(2*pi*X).*cos(2*pi*Y);

gradxe=-2*pi*sin(2*pi*X).*cos(2*pi*Y);
gradye=-2*pi*cos(2*pi*X).*sin(2*pi*Y);

[gradx,grady] = grad2(u);

figure(1)
surf(X,Y,gradx-gradxe)

figure(2)
surf(X,Y,grady-gradye)