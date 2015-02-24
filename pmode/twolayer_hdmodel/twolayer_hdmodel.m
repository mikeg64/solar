%compute energy flux through two layer model

%initial parameters

nx=128;
ny=128;
nz=128;

xmin=1953.1;
ymin=1953.1;
zmin=1.3333e5;

dx=(4.0e6-xmin)/128;
dy=dx;
dz=(6.0e6-zmin)/128;




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
omega=2*pi/period;

%driver mode
n1=0;
n2=0;
aa=500;%driver amplitude m/s


dd=0.5e6  ;%driver depth
nt=884;
dt=0.001;
%wd=zeros(nw,nx,ny,nz);


%compute pressure   
%compute lagrange displacement (differential with respect to time to get velpcit            
%compute loc flux
velocity=zeros(nx,ny,nz);
pressure=zeros(nx,ny,nz);
locflux=zeros(nx,ny,nz);
intflux=zeros(nz);
         
totflux=zeros(nt,nz);        





gamma=1.66666667;
m=1;
g=-274.0;

val3c=xlsread('atmos132.xls');

%set these values using val3c just loaded in
%photosphere/chromosphere
rho0p=val3c(8,3);
pres0p=val3c(8,4);
t0c=val3c(8,2);
hp=1.0/(rho0p*g);
csp=sqrt(gamma*pres0p/rho0p);

%corona
rho0c=val3c(51,3);
pres0c=val3c(51,4);
t0p=val3c(51,2);
hc=1.0/(rho0c*g);
csc=sqrt(gamma*pres0c/rho0c);

consts.rho0c=rho0c;
consts.pres0c=pres0c;
consts.t0c=t0c;
consts.hc=hc;

consts.rho0p=rho0p;
consts.pres0p=pres0p;
consts.t0p=t0p;
consts.hp=hp;


consts.csc=csc;
consts.csp=csp;
consts.gamma=gamma;
consts.g=g;

consts.dt=dt;
consts.omega=omega;

consts.nx=nx;
consts.ny=ny;
consts.nz=nz;

consts.dx=dx;
consts.dy=dy;
consts.dz=dz;

consts.xmin=xmin;
consts.ymin=ymin;
consts.zmin=zmin;

consts.aa=aa;

consts.n1=n1;
consts.n2=n2;

consts.dd=dd;


consts.l1=l1;
consts.l2=l2;


%the constants of integration for the evolution of
%lagrange displacement is defined from wave energy flux computation section
%see energyflux_integral.m
%warning check carfully the sign for l1 and l2
% [ d1c, d1p, d2c, d2p ] = dconsts( omega, consts, rho0c, rho0p, pres0c, pres0p, epl1, epl2, l1, l2 );
% 
% consts.d1c=d1c;
% consts.d2c=d2c;
% consts.d2c=d2c;
% consts.d2p=d2p;

%compute energy fluxes and computation
%of wave energy flux integral
%energyflux_integral;





