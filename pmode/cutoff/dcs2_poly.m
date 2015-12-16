function dcs2=dcs2_poly(pc.p, x) 
%polynomial fpr cutoff derived from polynomial expansions for speed of
%sound

dcs2=2*(9*pc.p1.*z.^8+8*pc.p2.*z.^7+7*pc.p3.*z.^6+6*pc.p4.*z.^5+5*pc.p5.*z.^4+4*pc.p6.*z**3+3*pc.p7.*z.^2+2*pc.p8.*z+pc.p9)*(pc.p1.*z.^9+pc.p2.*z.^8+pc.p3.*z.^7+pc.p4.*z.^6+pc.p5.*z.^5+pc.p6.*z.^4+pc.p7.*z.^3+pc.p8.*z.^2+pc.p9.*z+pc.p10);
