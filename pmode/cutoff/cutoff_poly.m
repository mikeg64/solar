function cutoff=cutoff_poly(pc, x) 
%polynomial fpr cutoff derived from polynomial expansions for speed of
%sound
omega=2.2833333378999998e+2/(pc.p1.*x.^9+pc.p2.*x.^8+pc.p3.*x.^7+pc.p4.*x.^6+pc.p5.*x.^5+pc.p6.*x.^4+pc.p7.*x.^3+pc.p8.*x.^2+pc.p9.*x+pc.p10);
    
    cutoff=2.*pi./omega;