function [ y ] = poly(cof, x  )
%Compute polynomial from 1 to 9
%   cof is structure with coefficients from p1 to p9
%x=x./meanx;
y=cof.p10+cof.p9.*x+cof.p8.*(x.^2)+cof.p7.*(x.^3)+cof.p6.*(x.^4)+cof.p5.*(x.^5)+cof.p4.*(x.^6)+cof.p3.*(x.^7)+cof.p2.*(x.^8)+cof.p1.*(x.^9);

end

