function [t, rho, pres]=solar_profiles_corona(z,rho0, pres0, t0,h)
%profile for coronal region
%h is scale height
    t=t0;
    rho=rho0.*exp(-z./h);
    pres=pres0.*exp(-z./h);

