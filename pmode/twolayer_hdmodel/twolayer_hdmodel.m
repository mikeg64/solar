%compute energy flux through two layer model

%initial parameters

nx=128;
ny=128;
nz=128;

dx=4.0e6/128;
dz=6.0e6/128;


%warning check carfully the sign for l1 and l2
l2=dz*(124-43);
l1=dz*42;

%state vars
%euleriandisp
%rho
%vx
%vy
%vz
%energy
%pressure
nw=6;
period=180.0;
om=2*pi/period;
nt=884;
dt=0.001;
wd=zeros(nw,nx,ny,nz);

gamma=1.66666667;
m=1;
g=-274.0;

%photosphere/chromosphere
rh0p=1;
pres0p=1;
hp=1.0/(rho0p*g);
csp=sqrt(gamma*pres0p/rho0p);

%corona
rh0c=1;
pres0c=1;
hc=1.0/(rho0*g);
csc=sqrt(gamma*pres0c/rho0c);

consts.csc=csc;
consts.csp=csp;
consts.gamma=gamma;
consts.g=g;




%warning check carfully the sign for l1 and l2
[ d1c, d1p, d2c, d2p ] = dconsts( omega, consts, rho0c, rho0p, pres0c, pres0p, epl1, epl2, l1, l2 );







