%build and compute solutions for the kg-equation with a stritified
%atmosphere
%See Taroyan and Erdelyi Sol. Phys. (2008) 251:523-531

loadatmos;

num=size(temp);
nr=num(1);
tri=1419; %index for position in transition layer (between 1324:1419)
TP=180; %driver period in seconds
%om=2*pi/TP;
om=1/TP;
iom=100;

kgc.TL=temp(tri);
kgc.T0=temp(2048);
kgc.L=height(tri);
kgc.a=(1-(kgc.TL/kgc.T0))/kgc.L; 
kgc.c0=cs(2048);
kgc.alpha=4*pi/(TP*kgc.c0*kgc.a);
kgc.beta=kgc.alpha*sqrt(1-kgc.a * kgc.L);
kgc.nu=0;
amp=350;


% A2=iom*cos(k*kgc.L)/(1-(kgc.L/(2*kgc.lam2))+((kgc.L/kgc.c2)*sqrt(kgc.om2.^2-om.^2)));
% A1=(iom*bessely(kgc.nu,kgc.beta)-A2*(cos(k*kgc.L)*bessely(kgc.nu,kgc.alpha))/(sqrt(1*kgc.a*kgc.L)))*(besseli(kgc.nu,kgc.alpha)*bessely(kgc.nu,kgc.beta)-besseli(kgc.nu,kgc.beta)*bessely(kgc.nu,kgc.alpha));
% B1=(iom-A1*besseli(kgc.nu,kgc.alpha))/bessely(kgc.nu,kgc.alpha);

% A2=iom*kgc.k1*exp(-1i*kgc.k2*kgc.L)/(exp(1i*kgc.k1*kgc.L)-sin(kgc.k1*kgc.L)*(1i*kgc.k1+1i*kgc.k2+((1/(2*kgc.lam2))-(1/(2*kgc.lam1)))));
% A1=(iom-A2*(exp(1i*(kgc.k1+kgc.k2)*kgc.L)))/(1-exp(2*1i*kgc.k1*kgc.L));
% B1=iom-A1;
% B2=0;

%compute the Q
for j1=1:100
    
   %TP=180; %driver period in seconds
   TP=j1*5;
%    TP=300;
om=1/TP;
iom=100;
    
%     iom=j*5;
    
kgc.c2=cs(1);
kgc.lam2=kgc.c2.^2/(consts.fgamma*consts.ggg);
kgc.om2=kgc.c2/(2*kgc.lam2);
omn=kgc.om2;
csn=kgc.c2;
%k=sqrt(om.^2-kgc.om2.^2)./kgc.c2;
    if om<omn
        k=(-1i).*sqrt(omn.*omn-om.*om)./csn;
    end
    
    if om>=omn
        k=sqrt(om.*om-omn.*omn)./csn;        
    end   
kgc.k2=k;

kgc.c1=kgc.c0;
kgc.lam1=kgc.c1.^2/(consts.fgamma*consts.ggg);
kgc.om1=kgc.c1/(2*kgc.lam1);
omn=kgc.om1;
csn=kgc.c1;
%k=sqrt(om.^2-kgc.om1.^2)./kgc.c1;

    if om<omn
        k=(-1i).*sqrt(omn.*omn-om.*om)./csn;
    end
    
    if om>=omn
        k=sqrt(om.*om-omn.*omn)./csn;        
    end  


kgc.k1=k;
    
    
    
    
for i1=nr:-1:1
    z=height(i1);
    
    if nr>tri   %region1
        kn=kgc.k1;
        %kgc.k1=kn;
        A2=amp*kgc.k1*exp(-1i*kgc.k2*kgc.L)/(exp(1i*kgc.k1*kgc.L)-sin(kgc.k1*kgc.L)*(1i*kgc.k1+1i*kgc.k2+((1/(2*kgc.lam2))-(1/(2*kgc.lam1)))));
        A1=(amp-A2*(exp(1i*(kgc.k1+kgc.k2)*kgc.L)))/(1-exp(2*1i*kgc.k1*kgc.L));
        B1=amp-A1;
        B2=0;
        s(i1,j1)=-om*kn*((A1.^2-B1.^2)+A1*B1*(exp( -2*1i*kn*z)-exp(2*1i*kn*z)));
    end
    if nr<=tri   %region 2 
        kn=kgc.k2;
        %kgc.k2=kn;
        A2=amp*kgc.k1*exp(-1i*kgc.k2*kgc.L)/(exp(1i*kgc.k1*kgc.L)-sin(kgc.k1*kgc.L)*(1i*kgc.k1+1i*kgc.k2+((1/(2*kgc.lam2))-(1/(2*kgc.lam1)))));
        A1=(amp-A2*(exp(1i*(kgc.k1+kgc.k2)*kgc.L)))/(1-exp(2*1i*kgc.k1*kgc.L));
        B1=amp-A1;
        B2=0;
        s(i1,j1)=-om*kn*((A2.^2-B2.^2));
    end
    
    %s(i1,j1)=sqrt(s(i1,j1).*conj(s(i1,j1)));
    s(i1,j1)=real(s(i1,j1));
end
end
figure;
h=surf(real(s'),'LineStyle','none');
% h=surf(s);
% hold on;

figure;
plot(5*(1:60),s(342,1:60)./s(1880,1:60),5*(1:60),s(684,1:60)./s(1880,1:60),5*(1:60),s(1367,1:60)./s(1880,1:60))

figure;

plot(height(:),s(:,6)./s(1880,6),height(:),s(:,36)./s(1880,36),height(:),s(:,60)./s(1880,60));

figure;

plot(5*(1:100),s(342,1:100))
