clear all; close all; clc;
rng(11); % setting seed for random numbers

meshSize = 64; % field size
windDir = [1, 0]; % ||windDir|| = 1
patchSize = 64;
A = 1e+4;
g = 9.81; % gravitational constant
windSpeed = 1e+2;

x1 = linspace(-10, 10, meshSize+1); x = x1(1:meshSize);
y1 = linspace(-10, 10, meshSize+1); y = y1(1:meshSize);
[X,Y] = meshgrid(x, y);

H0 = zeros(size(X)); % height field at time t = 0

for i = 1:meshSize
    for j = 1:meshSize
        kx = 2.0 * pi / patchSize * (-meshSize / 2.0 + x(i)); % = 2*pi*n / Lx
        ky = 2.0 * pi / patchSize * (-meshSize / 2.0 + y(j)); % = 2*pi*m / Ly
        P = phillips(kx, ky, windDir, windSpeed, A, g); % phillips spectrum
        H0(i,j) = 1/sqrt(2) * (randn(1) + 1i * randn(1)) * sqrt(P);
    end
end

H0 = H0 + conj(H0);

surf(X,Y,abs(ifft(H0)));
axis([-10 10 -10 10 -10 10]);


