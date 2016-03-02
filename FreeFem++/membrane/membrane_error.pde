verbosity =0;
real theta=4.*pi/3.;
real a=1.,b=1.;
int p=2;
border Gamma1(t=0,theta){ x = a * cos(t); y = b*sin(t); }
border Gamma2(t=theta,2*pi) { x = a * cos(t); y = b*sin(t); }
func f=-4*(cos(x^2+y^2-1) -(x^2+y^2)*sin(x^2+y^2-1));
func phiexact=sin(x^2+y^2-1);
real[int] L2error(p);
for(int n=0;n<p;n++)
	{
	mesh Th=buildmesh(Gamma1(20*(n+1))+Gamma2(10*(n+1)));
	fespace Vh(Th,P2);
	Vh phi,w;
	solve laplace(phi,w)=int2d(Th)(dx(phi)*dx(w) + dy(phi)*dy(w))
	- int2d(Th)(f*w) - int1d(Th,Gamma2)(2*w)+ on(Gamma1,phi=0);
	plot(Th,phi,wait=true,ps="membrane.eps");
	L2error[n]= sqrt(int2d(Th)((phi-phiexact)^2));
	}
for(int n=0;n<p;n++)
cout << " L2error " << n << " = "<<
L2error[n] <<endl;
cout <<" convergence rate = "<< log(L2error[0]/L2error[p-1])/log(2.^(p-1)) <<endl;
