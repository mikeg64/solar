function [ lagdisp, pres ] = lowerboundary( ix, iy, iz, it, velz, consts,  rho0 )
%drive the lower boundary of the model an oscillating membrane


    x=consts.xmin+ix*consts.dx;
    y=consts.ymin+iy*consts.dy;
    z=consts.zmin+iz*consts.dz;
    
    dz=consts.dz;

    xxmax=consts.xmin+const.nx*consts.dx;
    yymax=consts.ymin+consts.ny*consts.dy;
    %zzmax=consts.zmin+consts.nz*consts.dz;   
    
    n1=consts.n1;
    n2=consts.n2;
    
    dd=consts.dd;
    r1=(z-dd)*(z-dd);
    
    aa=const.aa;
    s_period=2*pi/(consts.omega);
    qt=consts.dt*it;
    
    exp_z=exp(-r1/(dz*dz));
    exp_xyz=sin(pi*y*(n1+1)/xxmax)*sin(pi*z*(n2+1)/yymax)*exp_z;

    tdep=sin(qt*2.0*pi/s_period);
    vvz=aa*exp_xyz*tdep;

    lagdisp=(velz+consts.dt*vvz)/omega;
    
    %pressure
    pres=lagdisp*consts.g*rho0;
    
    %w[fencode3_MODID(p,ii,mom1)]=w[fencode3_MODID(p,ii,mom1)]+(p->dt)*vvz*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
    %w[fencode3_MODID(p,ii,energy)]=w[fencode3_MODID(p,ii,energy)]+(p->dt)*vvz*vvz*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;


end

