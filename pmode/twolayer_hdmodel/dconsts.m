function [ d1c, d1p, d2c, d2p ] = dconsts( omega, consts, rho0c, rho0p, pres0c, pres0p, epl1, epl2, l1, l2 )
%compute constants for solutions of gravitationally stratified adibatic
%atmosphere
    gamma=consts.gamma;
    g=consts.g;
    csp=consts.csp;
    csc=consts.csc;

    rpc=(gamma*g/(2*csc*csc))+((omega/csc)*sqrt((gamma*g/(2*omega*csc))^2)-1);
    rmc=(gamma*g/(2*csc*csc))-((omega/csc)*sqrt((gamma*g/(2*omega*csc))^2)-1);

    rpp=(gamma*g/(2*csp*csp))+((omega/csp)*sqrt((gamma*g/(2*omega*csp))^2)-1);
    rmp=(gamma*g/(2*csp*csp))-((omega/csp)*sqrt((gamma*g/(2*omega*csp))^2)-1);

    a1p=g*(rho0c-rho0p)-(pres0c*rpc);
    a1m=g*(rho0c-rho0p)-(pres0c*rmc);

    a2p=pres0p*rpp;
    a2m=pres0p*rmp;

    al1=a1p-(a1m*exp(rpc/l2)/exp(rmc/l2));
    al2=a2p-(a2m*exp(-rpp/l1)/exp(-rmp/l1));
    al3=1-(exp(rpc/l2)/exp(rmc/l2));
    al4=-1+(exp(-rpp/l1)/exp(-rmp/l1));

    be1=(a1p*epl1*exp(rmc*l2)+a1m*epl2*exp(-rmp*l1))/(exp(rmc*l2-rmp*l1));
    be2=(epl1*exp(rmc*l2)+epl2*exp(-rmp*l1))/(exp(rmc*l2-rmp*l1));

    d1p=(al3*be1-al1*be2)/(al3*al2-al1*al4);
    d1c=(al4*be1-al2*be2)/(al4*al1-al2*al3);

    d2c=(epl2-d1c*exp(rpc*l2))/(exp(rmc*l2));
    d2p=(epl1-d1p*exp(-rpp*l1))/(exp(-rpp*l1));


end

