function [div] = div2(ux,uy)
global Px Py Qx Qy n
ux=reshape(ux,[],1);
uy=reshape(uy,[],1);
wx=Qx*ux;
dx=Px\wx;
wy=Qy*uy;
dy=Py\wy;
dx=reshape(dx,n,n);
dy=reshape(dy,n,n);
div=dx+dy;
end