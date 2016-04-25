simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7

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
       
       dx=(xmax-xmin)/(nx1-1);
       dy=(ymax-ymin)/(nx2-1);
       dz=(zmax-zmin)/(nx3-1);
       
       
       xx=zeros(nx1,nx2,nx3);
       yy=zeros(nx1,nx2,nx3);
       zz=zeros(nx1,nx2,nx3);
       
       
       for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   xx(i,j,k)=xmin+dx*(i-1);
                   yy(i,j,k)=ymin+dy*(j-1);
                   zz(i,j,k)=zmin+dz*(k-1);
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



%use fitting to extend atmos to 12.5Mm and 25Mm

tdens=dens(1:1269);
tpres=pres(1:1269);
theight=height(1:1269);
ttemp=temp(1:1269);


%
deltah=height(1)-height(2);
maxheight=12.5e6;
% nvals=(maxheight-height(2048))/deltah;
nvals=4271;
for i=nvals:-1:1
    nheight(i,1)=height(2048,1)+(nvals-i+1)*deltah;
end

ndens(3493:4271,1)=dens(1270:2048);
ntemp(3493:4271,1)=temp(1270:2048);
npres(3493:4271,1)=pres(1270:2048);

%compute values beyound transition region between 6.5Mm and 25Mm
%using data fitted with power law
for i=1:3492
    newh=nheight(i,1);
    ndens(i,1)=1.817e-7*newh.^(-0.667);
    npres(i,1)=6.717e-10*newh.^(1.219);
    ntemp(i,1)=2.669e-7*newh.^(1.886);
end




tempg=interp1(nheight,ntemp,xmin:dx:xmax);
presg=interp1(nheight,npres,xmin:dx:xmax);
densg=interp1(nheight,ndens,xmin:dx:xmax);



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


        
