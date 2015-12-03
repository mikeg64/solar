function cutoff=cutoff_poly_pres_rho(pr,pp, z) 
%polynomial fpr.p cutoff derived from polynomial expansions for speed of
%sound

omega=sqrt((1.37e+2.*(pr.p1.*z.^9+pr.p2.*z.^8+pr.p3.*z.^7+pr.p4.*z.^6+pr.p5.*z.^5+pr.p6.*z.^4+pr.p7.*z.^3+pr.p8.*z.^2+pr.p9.*z+pr.p10).*((1.66666667e+0.*(9.*pp.p1.*z.^8+8.*pp.p2.*z.^7+7.*pp.p3.*z.^6+6.*pp.p4.*z.^5+5.*pp.p5.*z.^4+4.*pp.p6.*z.^3+3.*pp.p7.*z.^2+2.*pp.p8.*z+pp.p9))/(pr.p1.*z.^9+pr.p2.*z.^8+pr.p3.*z.^7+pr.p4.*z.^6+pr.p5.*z.^5+pr.p6.*z.^4+pr.p7.*z.^3+pr.p8.*z.^2+pr.p9.*z+pr.p10)-(1.66666667e+0.*(9.*pr.p1.*z.^8+8.*pr.p2.*z.^7+7.*pr.p3.*z.^6+6.*pr.p4.*z.^5+5.*pr.p5.*z.^4+4.*pr.p6.*z.^3+3.*pr.p7.*z.^2+2.*pr.p8.*z+pr.p9).*(pp.p1.*z.^9+pp.p2.*z.^8+pp.p3.*z.^7+pp.p4.*z.^6+pp.p5.*z.^5+pp.p6.*z.^4+pp.p7.*z.^3+pp.p8.*z.^2+pp.p9.*z+pp.p10))/(pr.p1.*z.^9+pr.p2.*z.^8+pr.p3.*z.^7+pr.p4.*z.^6+pr.p5.*z.^5+pr.p6.*z.^4+pr.p7.*z.^3+pr.p8.*z.^2+pr.p9.*z+pr.p10).^2))/(pp.p1.*z.^9+pp.p2.*z.^8+pp.p3.*z.^7+pp.p4.*z.^6+pp.p5.*z.^5+pp.p6.*z.^4+pp.p7.*z.^3+pp.p8.*z.^2+pp.p9.*z+pp.p10)+(3.128166672923e+4.*(pr.p1.*z.^9+pr.p2.*z.^8+pr.p3.*z.^7+pr.p4.*z.^6+pr.p5.*z.^5+pr.p6.*z.^4+pr.p7.*z.^3+pr.p8.*z.^2+pr.p9.*z+pr.p10))/(pp.p1.*z.^9+pp.p2.*z.^8+pp.p3.*z.^7+pp.p4.*z.^6+pp.p5.*z.^5+pp.p6.*z.^4+pp.p7.*z.^3+pp.p8.*z.^2+pp.p9.*z+pp.p10));




    cutoff=2.*pi./omega;