simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;




%The config file 3D_128_spic_bin.ini is in the archive
%pmodeini.tgz available at the following link
%https://drive.google.com/open?id=0B-AKVl-pk6ziVUFqUzMycWJoODg



newfilename='3D_128_4Mm_asc.ini';
%newfilename='/data/cs1mkg/smaug/configs/3D_128_4Mm_bin.ini';

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=-274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7

ngx1=2;
ngx2=2;
ngx3=2;

  it=0;
       time=0;
       ndim=3;
       neqpar=7;
       nw=13;
       nx1=128;
       nx2=128;
       nx3=128;


      gamma=1.666667;
      adiab=1.0;
       eta=0;
       g1=-274;
       g2=0;
       g3=0;
       
       xmin=133333.33;
       ymin=1953.1;
       zmin=1953.1;
       xmax=5955555.6e0;
       ymax=4.0e6;
       zmax=4.0e6;
       
       dx=(xmax-xmin)/(nx1-2*ngx1);
       dy=(ymax-ymin)/(nx2-2*ngx2);
       dz=(zmax-zmin)/(nx3-2*ngx3);
       
       
       xx=zeros(nx1,nx2,nx3);
       yy=zeros(nx1,nx2,nx3);
       zz=zeros(nx1,nx2,nx3);

       rheight=zeros(nx1);
       
       
       for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   xx(i,j,k)=(xmin-ngx1*dx)+i*dx;
                   rheight(i)=(xmin-ngx1*dx)+i*dx;
                   yy(i,j,k)=(ymin-ngx1*dx)+dy*j;
                   zz(i,j,k)=(zmin-ngx1*dx)+dz*k;
               end
           end
       end
       
       
       
       simdata.w=zeros(nx1,nx2,nx3,nw);




        simparams.current_iteration=it;
        simparams.current_time=time;
        simparams.dimensionality=ndim;
        simparams.domain_dimensions=[nx1;nx2;nx3];
%         simparams.domain_left_edge=[0;0; 0.0];
%         simparams.domain_right_edge=[0; 0; 0];
        simparams.eta=eta;
        simparams.adiab=adiab;
%        field_ordering=1;
        simparams.gamma=gamma;
        simparams.gravity0=g1;
        simparams.gravity1=g2;
        simparams.gravity2=g3;
        
        
        
        
        simparams.domain_left_edge(1)=xx(1,1,1);
        simparams.domain_left_edge(2)=yy(1,1,1);
        simparams.domain_left_edge(3)=zz(1,1,1);
        
        simparams.domain_right_edge(1)=xx(nx1,1,1);
        simparams.domain_right_edge(2)=yy(1,nx2,1);
        simparams.domain_right_edge(3)=zz(1,1,nx3);

        
        simgridinfo.grid_dimensions(1)=nx1;
        simgridinfo.grid_dimensions(2)=nx2;
        simgridinfo.grid_dimensions(3)=nx3;
        
        %load atmosphere

      %% Import the data
data = xlsread('atmos.xls','VALMc_rho_2048_test');







%% Allocate imported array to column variable names
height = data(:,1);
temp = data(:,2);
dens = data(:,3);
pres = data(:,4);

cs=sqrt(consts.fgamma.*pres./dens); 

% dens=1e12*dens;
tempg=interp1(height,temp,(xmin-ngx1*dx):dx:(xmax+ngx1*dx));
presg=interp1(height,pres,(xmin-ngx1*dx):dx:(xmax+ngx1*dx));
densg=interp1(height,dens,(xmin-ngx1*dx):dx:(xmax+ngx1*dx));

ninterp=size(tempg);
nint=ninterp(2);

tempg(nint)=tempg(nint-2);
tempg(nint-1)=tempg(nint-2);
presg(nint)=presg(nint-2);
presg(nint-1)=presg(nint-2);
densg(nint)=densg(nint-2);
densg(nint-1)=densg(nint-2);

%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy


%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy
mu_thermal=0.6d0;
R=8.31e3;

%parrVALMc=rhoarrVALMc*TarrVALMc*R/mu_thermal
%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
% p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
% p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
%          (w[*,*,*,0]+w[*,*,*,9])/2.0
% p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
%           +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
% p[*,*,*]=(gamma-1.d0)*p[*,*,*]


%compute correct pressure for gravitationally stratified atmosphere

%compute initial energy (at photosphere or temperature minimum)
%mu_thermal=0.6d0;
%R=8.31e3;

% pressure=temp*R*density/((mu_thermal))
%parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
%iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)

% !iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)
% 
% !iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)
% 
% ! 1.6Mm
% 
iniene=6840.d0*R*(2.3409724e-09)/mu_thermal/(consts.fgamma-1.0);
% 
% !iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)
% 
% !iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)


%inix and inirho read from height versus density dta from valiiic mcwhirter
%data
%  do ix_2=ixGlo2,ixGhi2
%   do ix_3=ixGlo3,ixGhi3 
% 
%    x(ix_1,ix_2,ix_3,1)=inix !*1000.d0
%    w(ix_1,ix_2,ix_3,rho_)=inirho
%    w(ix_1,ix_2,ix_3,e_)=iniene
%    w(ix_1,ix_2,ix_3,m1_)=0.0
%    w(ix_1,ix_2,ix_3,m2_)=0.0
%    w(ix_1,ix_2,ix_3,m3_)=0.0   
% 
%   enddo
%  enddo


      for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   simdata.w(i,j,k,10)=densg(i);  %density
                   simdata.w(i,j,k,9)=iniene;  %energy
                                    
               end
           end
      end



% use energy to get pthermal

% p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=half*( w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,m1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,m2_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,m3_)&
%    **2 )/(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)&
%    +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))
% 
% p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=p(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3)+ half*( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)&
%    **2)+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_)**2)&
%    +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)**2) )&
%    +( (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))&
%    +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,bg3_)) )
% 
% p(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)=(eqpar(gamma_)&
%    -one)*(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)-p(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3))
for i=1:nx1
    presg(i)=(consts.fgamma-1)*iniene;
end

presg1=presg;



%compute correct pressure for gravitationally stratified atmosphere

% do ix_3=ixGlo3,ixGhi3
%  do ix_2=ixGlo2,ixGhi2
%   do ix_1=ixGhi1-1,ixGlo1,-1 
% 
% comi=-abs(x(ix_1+1,ix_2,ix_3,1)-x(ix_1,ix_2,ix_3,1))
% 
% w(ix_1,ix_2,ix_3,p_)=w(ix_1+1,ix_2,ix_3,p_)+w(ix_1,ix_2,ix_3,rho_)*comi*1.d0&
%    *eqpar(grav1_)
% 
% 
% 
%   enddo
%  enddo
% enddo


 


for i=nx1-1:-1:1
    comi=-abs(rheight(i+1)-rheight(i));
    presg(i)=presg(i+1)+densg(i)*comi*consts.ggg;
end


for i=3:nx1-2
     comi=-abs(rheight(i+1)-rheight(i));
     %densg(i)=densg(i)-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
     densg(i)=-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
end
%update density
% do ix_3=ixGlo3,ixGhi3
%  do ix_2=ixGlo2,ixGhi2
%   do ix_1=ixGlo1+2,ixGhi1-2
%        
%        w(ix_1,ix_2,ix_3,rho_)=-(1.D0/eqpar(grav1_))*(1.D0/(12.D0*(x(ix_1&
%           +1,ix_2,ix_3,1)-x(ix_1,ix_2,ix_3,1))))*(w(ix_1+2,ix_2,ix_3,p_)&
%           -8.D0*w(ix_1+1,ix_2,ix_3,p_)+8.D0*w(ix_1-1,ix_2,ix_3,&
%           p_)                         -w(ix_1-2,ix_2,ix_3,p_))
%                
% 
% 
%      enddo
%    enddo
%  enddo 




% !lower boundary
% do ix_1=ixmin1+4,ixmin1+2,-1
%         p_1=w(ix_1+2,ix_2,ix_3,p_)-8.d0*w(ix_1+1,ix_2,ix_3,p_)&
%            +8.d0*w(ix_1-1,ix_2,ix_3,p_)
%         p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
%         w(ix_1-2,ix_2,ix_3,p_) = 12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1&
%            -1,ix_2,ix_3,1))*p_2+p_1
 for i=5:-1:3
   p_1=presg(i+2)+8*presg(i+1)-8*presg(i-1);
   p_2=densg(i)*consts.ggg;
   presg(i-2)= -p_1-12.0*(rheight(i)-rheight(i-1))*p_2;
 end




% !upper boundary
% do ix_1=ixmax1-4,ixmax1-2
%    do ix_2=ixmin2,ixmax2
%       do ix_3=ixmin3,ixmax3
%          
%           p_1=w(ix_1-2,ix_2,ix_3,p_)-8.d0*w(ix_1-1,ix_2,ix_3,p_)+8.d0*w(ix_1&
%              +1,ix_2,ix_3,p_)
%           p_2=w(ix_1,ix_2,ix_3,rho_)*eqpar(grav1_)
%           w(ix_1+2,ix_2,ix_3,p_) = -12.d0*(x(ix_1,ix_2,ix_3,1)-x(ix_1&
%              -1,ix_2,ix_3,1))*p_2+p_1
% 
% !           p_1=w(ix_1-2,ix_2,p_)-8.d0*w(ix_1-1,ix_2,p_)+8.d0*w(ix_1+1,ix_2,p_)
% !           p_2=w(ix_1,ix_2,rho_)*eqpar(grav1_)
% !           w(ix_1+2,ix_2,p_) = -12.d0*(x(ix_1,ix_2,1)-x(ix_1-1,ix_2,1))*p_2+p_1
% 
%       enddo
%    enddo
% enddo
for i=nx1-4:nx1-2
  p_1=presg(i-2)-8*presg(i-1)+8*presg(i+1);
  p_2=densg(i)*consts.ggg;
  presg(i+2)= p_1-12.0*(rheight(i)-rheight(i-1))*p_2;
end











presg2=presg;



%compute enrgy using pressure
% ! Calculate total energy from pressure, kinetic and magnetic energy
% 
% w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,e_)=w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,p_)/(eqpar(gamma_)-1)+half*((w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rhob_))&
%    *(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,v1_)**2+w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,v2_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,v3_)**2)+((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_))&
%    **2+(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b2_))**2&
%    +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_))**2))&
%    +((w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b1_)*w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,bg1_))+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
%    ixmin3:ixmax3,b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,bg2_))&
%    +(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b3_)*w(ixmin1:ixmax1,&
%    ixmin2:ixmax2,ixmin3:ixmax3,bg3_)))
for i=1:nx1
    energg(i)=presg(i)/(consts.fgamma -1);
end



%e=p/(rho*(gamma-1.d0))+0.5d0*(bx*bx+bz*bz)
% for i=4:nx1-3   
%     tempg(i)=(tempg(i-3)+tempg(i-2)+tempg(i-1)+tempg(i+1)+tempg(i+2)+tempg(i+3))/6;
% end







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

% writesac3D(newfilename, simparams, simgridinfo, simdata, 'ascii');
        
