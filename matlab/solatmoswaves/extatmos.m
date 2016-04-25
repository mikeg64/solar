
%
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
for i=1:3492
    newh=nheight(i,1);
    ndens(i,1)=1.817e-7*newh.^(-0.667);
    npres(i,1)=6.717e-10*newh.^(1.219);
    ntemp(i,1)=2.669e-7*newh.^(1.886);
end


       nx1=128;
       nx2=128;
       nx3=128;


       xmin=133333.33;
       ymin=1953.1;
       zmin=1953.1;
       xmax=5955555.6e0;
       ymax=4.0e6;
       zmax=4.0e6;
       
       dx=(xmax-xmin)/(nx1-1);
       dy=(ymax-ymin)/(nx2-1);
       dz=(zmax-zmin)/(nx3-1);




tempg=interp1(nheight,ntemp,xmin:dx:xmax);
presg=interp1(nheight,npres,xmin:dx:xmax);
densg=interp1(nheight,ndens,xmin:dx:xmax);

%ndens=nheight