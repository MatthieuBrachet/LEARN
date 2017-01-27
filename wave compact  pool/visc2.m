function [visc] = visc2(h)
global lap n
h=reshape(h,[],1);
visc=lap*h;
visc=reshape(visc,n,n);
end