function [ ld ] = lagrange_disp(z,consts, l1, l2, d1c, d1p, d2c, d2p )
%Compute Lagrangian displacement
%   
    gamma=consts.gamma;
    g=consts.g;
    
    
    if z>l1  
        csc=consts.csc;
        rpc=(gamma*g/(2*csc*csc))+((omega/csc)*sqrt((gamma*g/(2*omega*csc))^2)-1);
        rmc=(gamma*g/(2*csc*csc))-((omega/csc)*sqrt((gamma*g/(2*omega*csc))^2)-1);        
        ldt=d1c*exp(rpc*z)+d2c*exp(rmc*z);
    else
        csp=consts.csp;
        rpp=(gamma*g/(2*csp*csp))+((omega/csp)*sqrt((gamma*g/(2*omega*csp))^2)-1);
        rmp=(gamma*g/(2*csp*csp))-((omega/csp)*sqrt((gamma*g/(2*omega*csp))^2)-1);
        ldt=d1p*exp(rpp*z)+d2p*exp(rmp*z);        
    end
    ld=ldt;
    
    

end

