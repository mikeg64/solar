newfilename='3D_128_2p5_2p5_12p5_asc.ini';


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
       
       %xmin=133333.33;
       xmin=199219.0;
       %ymin=1953.1;
       ymin=39687.5;
       %zmin=1953.1;
       zmin=39687.5;
       %xmax=5955555.6e0;
       xmax=12.8496e6;
       %ymax=4.0e6;
       %zmax=4.0e6;
       ymax=2.559984e6;
       zmax=2.559984e6;
       
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
maxheight=12.8496e6;
% nvals=(maxheight-height(2048))/deltah;
nvals=4392;
for i=nvals:-1:1
    nheight(i,1)=height(2048,1)+(nvals-i+1)*deltah;
end

ndens(3614:4392,1)=dens(1270:2048);
ntemp(3614:4392,1)=temp(1270:2048);
npres(3614:4392,1)=pres(1270:2048);


%load the results from the fitting
load('dens_corona_fittedmodel.mat');
load('temp_corona_fittedmodel.mat');
load('pres_corona_fittedmodel.mat');
dens_corona=cfit(dens_corona_fittedmodel);
temp_corona=cfit(temp_corona_fittedmodel);
pres_corona=cfit(pres_corona_fittedmodel);


%compute values beyound transition region between 6.5Mm and 25Mm
%using data fitted with power law
for i=1:3613
    newh=nheight(i,1);
%using matlab fitting functions
     ndens(i,1)=dens_corona(newh);
     npres(i,1)=pres_corona(newh);
     ntemp(i,1)=temp_corona(newh);

%old power law    
%     ndens(i,1)=1.817e-7*newh.^(-0.667);
%     npres(i,1)=6.717e-10*newh.^(1.219);
%     ntemp(i,1)=2.669e-7*newh.^(1.886);
end



energg=interp1(nvals,nenerg,xmine:dxe:xmaxe);
tempg=interp1(nheight,ntemp,xmin:dx:xmax);
presg=interp1(nheight,npres,xmin:dx:xmax);
densg=interp1(nheight,ndens,xmin:dx:xmax);
energ=zeros(1,nx1);


%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy
mu=0.6d0;
R=8.31e3;

%parrVALMc=rhoarrVALMc*TarrVALMc*R/mu

% p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
% p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
%          (w[*,*,*,0]+w[*,*,*,9])/2.0
% p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
%           +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
% p[*,*,*]=(gamma-1.d0)*p[*,*,*]


%iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)


%e=p/(rho*(gamma-1.d0))+0.5d0*(bx*bx+bz*bz)
for i=4:nx1-3   
    tempg(i)=(tempg(i-3)+tempg(i-2)+tempg(i-1)+tempg(i+1)+tempg(i+2)+tempg(i+3))/6;
end


      for i=1:nx1
           for j=1:nx2
               for k=1:nx3
                   simdata.w(i,j,k,10)=densg(i);  %density
                   presg(i)=tempg(i)*densg(i)*R/mu;
                   energ(i)=presg(i)./(densg(i)*(consts.fgamma-1));
                   simdata.w(i,j,k,9)=energ(i);
               end
           end
       end

% writesac3D(newfilename, simparams, simgridinfo, simdata, 'ascii');
        
