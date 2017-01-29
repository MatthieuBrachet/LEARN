function [convx,convy] = convec2(ux,uy)
global Px Py Qx Qy n
ux=reshape(ux,[],1);
uy=reshape(uy,[],1);
wx=Qx*ux;
duxdx=Px\wx;
wy=Qy*uy;
duydy=Py\wy;
wy=Qx*uy;
duydx=Px\wy;
wx=Qy*ux;
duxdy=Py\wy;

convx=ux.*duxdx+uy.*duxdy;
convy=ux.*duydx+uy.*duydy;

convx=reshape(convx,n+1,n+1);
convy=reshape(convy,n+1,n+1);
end