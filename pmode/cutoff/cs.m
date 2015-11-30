function [ csv ] = cs( consts, pres, rho )
%Compute speed of sound in an ideal gas
%   Pressure=pres, Density=rho, consts.fagamma=adiabatic parameter
    csv=sqrt(consts.fgamma.*pres./rho);
end

