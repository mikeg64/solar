function [ y ] = poly(cof, x  )
%Compute polynomial from 1 to 9
%   cof is structure with coefficients from p1 to p9

y=cof.p9+cof.p8.*x+cof.p7.*x.^2+cof.p6.*x.^3+cof.p5.*x.^4+cof.p4.*x.^5+cof.p3.*x.^6+cof.p2.*x.^7+cof.p1.*x.^8;

end

