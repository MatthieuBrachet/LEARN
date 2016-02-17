int n=100;
border C(t=0,2*pi){x=cos(t);y=sin(t);}
mesh Th = buildmesh(C(n));
fespace Vh(Th,P2);
Vh u, v;
func f=x*y;
problem Poisson(u,v,solver=LU) =
  int2d(Th)(dx(u)*dx(v)+dy(u)*dy(v))
  - int2d(Th)(f*v)
  + on(C,u=0);
real cpu=clock();
Poisson;
plot(u);
cout << " CPU time = " << clock() - cpu << endl;