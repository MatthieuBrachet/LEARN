function [ vort ] = vort2( ux,uy )
global n
global Qx Qy Px Py
ux=reshape(ux,[],1);
uy=reshape(uy,[],1);
w1=Qx*uy;
w2=Qy*ux;
vort=Px\w1+Py\w2;
vort=reshape(vort,n+1,n+1);
end

