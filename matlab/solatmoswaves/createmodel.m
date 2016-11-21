simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
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
       
       
       for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   xx(i,j,k)=(xmin-ngx1*dx)+i*dx;
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


%e=p/(rho*(gamma-1.d0))+0.5d0*(bx*bx+bz*bz)

      for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   simdata.w(i,j,k,10)=densg(i);  %density
                   simdata.w(i,j,k,9)=presg(i)./(densg(i)*(consts.fgamma-1));
               end
           end
       end

%writesac3D(newfilename, simparams, simgridinfo, simdata, 'ascii');
        
