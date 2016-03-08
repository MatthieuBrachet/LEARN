load "msh3"
 
border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2= buildmesh(cc(150));

plot(Th2, wait=1);
 
// eps is there to avoid square roots of negative values
// and thin triangles close to the equator
real eps=1e-3;
 
func zmin=-sqrt(1+eps-x*x-y*y);
func zmax=sqrt(1+eps-x*x-y*y);
 

int[int] rup=[0,2],rdown=[0,1],rmid=[1,3,2,3,3,3,4,3];
 
mesh3 Th=buildlayers(Th2,20,
  zbound=[zmin,zmax],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);

plot(Th, wait=1);
plot(Th,Th2);
