simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;




%The config file 3D_128_spic_bin.ini is in the archive
%pmodeini.tgz available at the following link
%https://drive.google.com/open?id=0B-AKVl-pk6ziVUFqUzMycWJoODg



newfilename='2D_hurlb_bin.ini';
%newfilename='/data/cs1mkg/smaug/configs/3D_128_4Mm_bin.ini';

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7


gamma=1.6666667;                     % -gamma (5/3)
mpoly=1;                             % -mpoly
chi=11.0;                          % -chi
qchand=72.0;                          % -qchand
dvx=1.0e-4;                         % -dvx
dvy=1.0e-4;                         % -dvy
nkx=2.0943951;                     % -nkx (2*pi/3)
nky=6.2831531;                     % -nky (2*pi)
sigma=1.0;                           % -sigma
rhat=1.0e5;                          % -rhat
zeta0=0.25;                          % -zeta0

qmpoly=mpoly;

zz0=1./(qchi.^(1./qmpoly)-1);
%write(*,*)'Top layer dimensionless temperature z0=',zz0
temptop_=zz0;

gamma_=gamma;
grav1_= 0;
grav2_=-(qmpoly+1.0);


eta2=(qmpoly+1.0)*(gamma/(gamma-1)-(qmpoly+1))*((gamma-1.0)/gamma)*((zz0+0.5)**(2*qmpoly-1.0)/zz0.^(2*qmpoly))*(zeta0.^2/(sigma*rhat));




ngx1=2;
ngx2=2;
%ngx3=2;

  it=0;
       time=0;
       ndim=2;
       neqpar=5;
       nw=10;
       nx1=200;
       nx2=40;
     %  nx3=128;


      gamma=1.666667;
      adiab=1.0;
       eta=0;
       g1=-274;
       g2=0;
       %g3=0;
       
       xmin=0.0;
       ymin=0.0;
       %zmin=1953.1;
       xmax=3.0;
       ymax=1.0;
       %zmax=4.0e6;
       
       dx=(xmax-xmin)/(nx1-2*ngx1);
       dy=(ymax-ymin)/(nx2-2*ngx2);
       %dz=(zmax-zmin)/(nx3-2*ngx3);
       
       
       xx=zeros(nx1,nx2);
       yy=zeros(nx1,nx2);
       %zz=zeros(nx1,nx2,nx3);

       rheight=zeros(nx1);
       
       
       for i=1:nx1
           for j=1:nx2
               %for k=1:nx3
                   xx(i,j)=(xmin-ngx1*dx)+i*dx;
                   rheight(i)=(xmin-ngx1*dx)+i*dx;
                   yy(i,j)=(ymin-ngx1*dx)+dy*j;
                 %  zz(i,j,k)=(zmin-ngx1*dx)+dz*k;
               %end
           end
       end
       
       
       
       simdata.w=zeros(nx1,nx2,nw);




        simparams.current_iteration=it;
        simparams.current_time=time;
        simparams.dimensionality=ndim;
        simparams.domain_dimensions=[nx1;nx2];
%         simparams.domain_left_edge=[0;0; 0.0];
%         simparams.domain_right_edge=[0; 0; 0];
        simparams.eta=eta;
        simparams.adiab=adiab;
%        field_ordering=1;
        simparams.gamma=gamma;
        simparams.gravity0=g1;
        simparams.gravity1=g2;
       % simparams.gravity2=g3;
        
        
        
        
        simparams.domain_left_edge(1)=xx(1,1);
        simparams.domain_left_edge(2)=yy(1,1);
        %simparams.domain_left_edge(3)=zz(1,1,1);
        
        simparams.domain_right_edge(1)=xx(nx1,1);
        simparams.domain_right_edge(2)=yy(1,nx2);
        %simparams.domain_right_edge(3)=zz(1,1,nx3);

        
        simgridinfo.grid_dimensions(1)=nx1;
        simgridinfo.grid_dimensions(2)=nx2;
        %simgridinfo.grid_dimensions(3)=nx3;
        
        %load atmosphere

      %% Import the data

cs=sqrt(consts.fgamma.*pres./dens); 



%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy
mu=0.6d0;
R=8.31e3;






qmpoly=mpoly;

zz0=1./(qchi.^(1/qmpoly)-1);

temptop_=zz0;

gamma_=gamma;
grav1_= 0;
grav2_=-(qmpoly+one);


% calculate eqpar-array from input parameters

eta2=(qmpoly+1.0)*(gamma/(gamma-1)-(qmpoly+1))*((gamma-1.0)/gamma)*((zz0+0.5)**(2*qmpoly-1.0)/zz0**(2*qmpoly))*(zeta0.^2/(sigma*rhat));
if(eta2<=eps)then
   write(*,*)'Negative or too small value for eta**2:',eta2
end

eqpar(eta_)=sqrt(eta2)
eqpar(nu_)=eqpar(eta_)*sigma/zeta0
eqpar(kappa_)=(gamma/(gamma-one))*eqpar(eta_)/zeta0

write(*,*)'dimensionless values for dissipative coefficients:'
write(*,*)'resistivity          eta=',eqpar(eta_)
write(*,*)'viscosity             nu=',eqpar(nu_)
write(*,*)'thermal conduction kappa=',eqpar(kappa_)

bstrength=dsqrt(qchand*eqpar(nu_)*eqpar(eta_))
write(*,*)'dimensionless magnetic field strength:',bstrength

!
! set polytropic stratification
! assume stratification gravity along -e_y
!

w(ix^S,rho_)= ((zz0+one-x(ix^S,2))/zz0)**qmpoly
w(ix^S,b1_)= zero
w(ix^S,b2_)=bstrength
! set pressure
w(ix^S,e_)= zz0*(((zz0+one-x(ix^S,2))/zz0)**(qmpoly+one))

!
! small velocity perturbations
!
! amplitudes and wavenumbers from input
!
w(ix^S,m1_)=dvx*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky)
w(ix^S,m2_)=dvy*dsin(x(ix^S,1)*nkx)*dsin(x(ix^S,2)*nky)





















      for i=1:nx1
           for j=1:nx2
               %for k=1:nx3
                   simdata.w(i,j,10)=densg(i);  %density
                   simdata.w(i,j,9)=iniene;  %energy
                                    
               %end
           end
      end




for i=nx1-1:-1:1
    comi=-abs(rheight(i+1)-rheight(i));
    presg(i)=presg(i+1)-densg(i)*comi*consts.ggg;
end


for i=3:nx1-2
     comi=-abs(rheight(i+1)-rheight(i));
     %densg(i)=densg(i)-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
     densg(i)=(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
end

for i=1:nx1
    energg(i)=presg(i)/(consts.fgamma -1);
end







disp('rebuilding array');
      for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   
                   for iw=1:13
                      simdata.w(i,j,k,iw)=0.0; 
                   end
                   simdata.w(i,j,k,10)=densg(i);  %density
%                    presg(i)=tempg(i)*densg(i)*R/mu;
%                    energ(i)=presg(i)./(densg(i)*(consts.fgamma-1));
                   simdata.w(i,j,k,9)=energg(i);
                   %simdata.w(i,j,k,9)=presg(i)./(densg(i)*(consts.fgamma-1));

               end
           end
      end

disp('generate field');
%
%[simparams, simgridinfo, simdata]=generatefield(simparams, simgridinfo, simdata, 'fluxtube');

writesac2D(newfilename, simparams, simgridinfo, simdata, 'binary');
        
