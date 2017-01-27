function [gradx,grady] = grad2(u)
global Px Py Qx Qy n
u=reshape(u,[],1);
wx=Qx*u;
gradx=Px\wx;
wy=Qy*u;
grady=Py\wy;
gradx=reshape(gradx,n+1,n+1);
grady=reshape(grady,n+1,n+1);
end