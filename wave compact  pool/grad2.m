function [gradx,grady] = grad2(u)
global Px Py Qx Qy n
u=reshape(u,[],1);
wx=Qx*u;
gradx=Px\wx;
wy=Qy*u;
grady=Py\wy;
gradx=reshape(gradx,n,n);
grady=reshape(grady,n,n);
end