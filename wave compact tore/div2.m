function [div] = div2(ux,uy)
global Px Py Qx Qy n
ux=reshape(ux,[],1);
uy=reshape(uy,[],1);
wx=Qx*ux;
dx=Px\wx;
wy=Qy*uy;
dy=Py\wy;
dx=reshape(dx,n+1,n+1);
dy=reshape(dy,n+1,n+1);
div=dx+dy;
end