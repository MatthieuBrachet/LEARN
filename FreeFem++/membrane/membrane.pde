real theta=4.*pi/3.;
real a=2.,b=1.;
func z=100*cos(x);
border Gamma1(t=0,theta){ x = a * cos(t); y = b*sin(t); }
border Gamma2(t=theta,2*pi){ x = a * cos(t); y = b*sin(t); }
mesh Th=buildmesh(Gamma1(100)+Gamma2(50));
fespace Vh(Th,P2);
Vh phi,w, f=1;
solve Laplace(phi,w)=int2d(Th)(dx(phi)*dx(w) + dy(phi)*dy(w))
      - int2d(Th)(f*w) + on(Gamma1,phi=z);
plot(Th,wait=true, ps="membraneTh.eps");
plot(phi,fill=true, ps="membrane.eps"); // representation of plotting; wait is an orther option.
savemesh(Th,"Th.msh");