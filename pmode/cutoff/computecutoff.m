%%
% $\frac{\partial^2 Q}{\partial t^2} - c_t^2(z) \frac{\partial^2 Q}{\partial z^2} + \Omega^2(z)Q = 0$
% 
% $$\omega_{c}=\frac{\gamma g}{4\pi c_{s}}\sqrt{1+2\frac{d}{dz}\frac{P}{\rho g}}$$
% 
% $$\frac{d}{dz}\frac{P}{\rho g}$$
% 
% 


loadatmos;
nsmoothsteps=4;

loadcutoffdata;

[nr,nc]=size(pres);

plot(height(1420:nr),pres(1420:nr),height(1420:nr),pfit1);
plot(height(1420:nr)./1e6,cutoff_chromos(height(1420:nr)),hc./1e6,cutoff,height(1324:1419)./1e6,cutoff_transition(height(1324:1419)));
plot(height(1420:nr)./1e6,cutoff_chromos(height(1420:nr)),hc./1e6,cutoff,height(1324:1419)./1e6,cutoff_transition(height(1324:1419)),height(1:1325)./1e6,cutoff_corona(height(1:1325)));


[pressmooth, rhosmooth]=smoothav(consts, pres, dens,nsmoothsteps);

csav=sqrt(consts.fgamma.*pressmooth./rhosmooth);
plot(height./1e6,csav./1e3);

lam0=pressmooth./(rhosmooth.*consts.ggg);
lamdashs0=zeros(nr);
lamdash0=zeros(nr);
atc0=zeros(nr);

dh=height(1)-height(2);    

for j=1:nsmoothsteps
   lam0(nr+j)=lam0(nr);   
end

for i=3:nr
    lamdash0(i)=-diff5(lam0,i,dh);
    if lamdash0(i)<-200
        lamdash0(i)=-200;
    end
end

for i=1:nr
    sdashtot=0;
    for j=1:nsmoothsteps-1
        sdashtot=sdashtot+lamdash0(i+j);    
    end
    lamdashs0(i)=sdashtot/nsmoothsteps;
end

for i=1:nr
    atc0(i)=1.0/((consts.fgamma.*consts.ggg/(4*pi.*csav(i)))*sqrt(1+2*lamdashs0(i)));
end

plot(height./1e6,atc0);
