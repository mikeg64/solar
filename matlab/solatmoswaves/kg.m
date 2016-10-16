%build and compute solutions for the kg-equation with a stritified
%atmosphere
%See Taroyan and Erdelyi Sol. Phys. (2008) 251:523-531

loadatmos;


tri=1419; %index for position in transition layer (between 1324:1419)
TP=180; %driver period in seconds
om=2*pi/TP;
iom=100;

kgc.TL=temp(tri);
kgc.T0=temp(2048);
kgc.L=height(tri);
kgc.a=(1-(kgc.TL/kgc.T0))/kgc.L; 
kgc.c0=cs(2048);
kgc.alpha=4*pi/(TP*kgc.c0*kgc.a);
kgc.beta=kgc.alpha*sqrt(1-kgc.a * kgc.L);
kgc.nu=0;

kgc.c2=cs(1);
kgc.lam2=kgc.c2.^2/(consts.fgamma*consts.ggg);
kgc.om2=kgc.c2/(2*kgc.lam2);
k=sqrt(om.^2-kgc.om2.^2)./kgc.c2;

kgc.k2=k;

kgc.c1=kgc.c0;
kgc.lam1=kgc.c1.^2/(consts.fgamma*consts.ggg);
kgc.om1=kgc.c1/(2*kgc.lam1);
k=sqrt(om.^2-kgc.om1.^2)./kgc.c1;
kgc.k1=k;

% A2=iom*cos(k*kgc.L)/(1-(kgc.L/(2*kgc.lam2))+((kgc.L/kgc.c2)*sqrt(kgc.om2.^2-om.^2)));
% A1=(iom*bessely(kgc.nu,kgc.beta)-A2*(cos(k*kgc.L)*bessely(kgc.nu,kgc.alpha))/(sqrt(1*kgc.a*kgc.L)))*(besseli(kgc.nu,kgc.alpha)*bessely(kgc.nu,kgc.beta)-besseli(kgc.nu,kgc.beta)*bessely(kgc.nu,kgc.alpha));
% B1=(iom-A1*besseli(kgc.nu,kgc.alpha))/bessely(kgc.nu,kgc.alpha);

A2=iom*kgc.k1*exp(-j*kgc.k2*kgc.L)/(exp(j*kgc.k1*kgc.L)-sin(kgc.k1*kgc.L)*(j*kgc.k1+j*kgc.k2+((1/(2*kgc.lam2))-(1/(2*kgc.lam1)))));








