function [t, rho, pres]=solar_profiles_phot(z,rho0, pres0, t0,z0,m)
%profile for photospheric/chromospheric region
    fac=1+(z./z0);
    t=t0.*fac;
    rho=rho0.*(fac.^m);
    pres=pres0.*(fac).^(m+1);