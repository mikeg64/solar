
%directory='/storage2/mikeg/results/spic6b0_1_3d/';
%directory='/storage2/mikeg/results/spic5b0_b1G_3d/';
%directory='/storage2/mikeg/results/spic4b0_3_3d/';
%directory='/storage2/mikeg/results/spic5b0_3d/';
%directory='/storage2/mikeg/results/spic3p0a_0_2_3d/';
%directory='/storage2/mikeg/results/spic4p3a_0_1_3d/';
%directory='/storage2/mikeg/results/spic6p7a_0_0_3d/';
%directory='/storage2/mikeg/results/spic2p3a_0_3_3d/';
directory='/data/cs1mkg/smaugtest/smaug/out/';
extension='.out';

%ndirectory='/storage2/mikeg/results/spic5b0_b1G_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic6b0_1_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic4b0_3_3d/images_3d_secs/';
%ndirectory='/storage2/mikeg/results/spic5b0_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic3p0a_0_2_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic4p3a_0_1_3d/images_3d_vsecs/';
%ndirectory='/storage2/mikeg/results/spic6p7a_0_0_3d/images_3d_vsecs/';
ndirectory='/data/cs1mkg/smaugtest/smaug/out/images_3d_vsecs/';
nextension='.jpg';

for i=1:1:60
%for i=1519:2632
%for i=2631:2632
    

id=int2str(1000*i);
filename=[directory,'sedov_',id,extension];
timetext=['time=',int2str(i),'s'];
imfile=[ndirectory,'im1_',id,nextension];
disp([id filename]);
   fid=fopen(trim(filename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=trim(setstr(fread(fid,79,'char')'));
   it=fread(fid,1,'integer*4'); time=fread(fid,1,'float64');
 
   ndim=fread(fid,1,'integer*4');
   neqpar=fread(fid,1,'integer*4'); 
   nw=fread(fid,1,'integer*4');
   nx=fread(fid,3,'integer*4');
   
   nxs=nx(1)*nx(2)*nx(3);
   varbuf=fread(fid,7,'float64');
   
   gamma=varbuf(1);
   eta=varbuf(2);
   g(1)=varbuf(3);
   g(2)=varbuf(4);
   g(3)=varbuf(5);
   
   
   varnames=trim(setstr(fread(fid,79,'char')'));
   
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   
   for iw=1:nw
      %fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      %fread(fid,4);
   end
   
   nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);
   
   xx=reshape(X(:,1),nx1,nx2,nx3);
   yy=reshape(X(:,2),nx1,nx2,nx3);
   zz=reshape(X(:,3),nx1,nx2,nx3);
   
   
 
  % extract variables from w into variables named after the strings in wnames
wd=zeros(nw,nx1,nx2,nx3);
for iw=1:nw
  
     tmp=reshape(w(:,iw),nx1,nx2,nx3);
     wd(iw,:,:,:)=tmp;
end


%w=tmp(iw);
  

clear tmp; 
   
   
   fclose(fid);




%plot sections through 3d array
   %slice=48;
   x=linspace(0,4,128);
   y=linspace(0,4,128);
   z=linspace(0,6,128);
   
   
   
   
   nrange=3:126;
   
   ax=x(nrange);
   ay=y(nrange);
   az=z(nrange);
   [x1,x2,x3] = meshgrid(ax,ay,az);
   %val1=reshape(wd(2,nrange,nrange,nrange),124,124,124);
   %val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   %myval=shiftdim(val1./val2,1);

   val1=reshape(wd(5,nrange,nrange,nrange),124,124,124);
   %val2=reshape(wd(1,nrange,nrange,nrange)+wd(10,nrange,nrange,nrange),124,124,124);
   myval=shiftdim(val1,1);



	R=8.3e+003;
	mu=1.257E-6;
	mu_gas=0.6;
	gamma=1.66667;

%sabx=reshape(wd(11,nrange,nrange,nrange),124,124,124);
%saby=reshape(wd(12,nrange,nrange,nrange),124,124,124);
%sabz=reshape(wd(13,nrange,nrange,nrange),124,124,124);
%TP=reshape(wd(9,nrange,nrange,nrange),124,124,124);
%TP=TP-(sabx.^2.0+saby.^2.0+sabz.^2.0)./2.0;
%TP=(gamma-1.d0).*TP;



   %mval is T
  
  % myval=shiftdim(mu_gas.*TP./R./val2,1);




   %P = [2 1 3];
   %x1 = permute(x1, P);
   %x2 = permute(x2, P);
   %x3 = permute(x3, P);
   %myval = permute(myval, P);
   
   
  %h= slice(myval,64, 64, 4);
  %figure('Visible','off','IntegerHandle','Off');
  %hold on;
  %h=slice(myval,80, 64,8);
  %h=slice(myval,96, 96,[5 49 100]);  %used for 0,1 mode
  %h=slice(myval,108, 96,[5 49 100]);  %used for 0,1 mode
  %h=slice(myval,64, 64,[5 64 100]);  %used for 0,0 mode
  h=slice(myval,63, 63,63);  %used for 0,0 mode
  %h=slice(myval,65, 65,[5 49 100]);
  %h=slice(myval,105, 96,8);
  hold on;
  set(h,'EdgeColor','none','FaceColor','interp');
  %grid off;
  %set(h,'XData',ax,'YData',ay,'ZData',az);
  hax=get(h,'Children');
  set(gca,'CameraPosition',[-606.298 -914.02 280.537]);
  set(gca,'Xlim',[0 124],'Ylim',[0 124],'Zlim',[0 124]);
  
  %set(gca,'XTickLabel',{'0';'1.6';'3.2'})
  %set(gca,'YTickLabel',{'0';'1.6';'3.2'})



  %set(gca,'YTickLabel',{'0';'1.6';'3.2'})
  %set(gca,'XTickLabel',{'0';'0.63';'1.26';'1.89';'2.52';'3.15';'3.78'})
  %set(gca,'ZTickLabel',{'0.09';'0.99';'1.94';'2.88';'3.83';'4.77';'5.72'})
  %cmap=colormap('Jet');
  
  
  
  
   max1=max(myval);
  max2=max(max1);
  max3=max(max2);
  
  min1=min(myval);
  min2=min(min1);
  min3=min(min2);
  
  maxval=5;
  minval=0;
  
  if min3<minval
      minval=min3;
  end
  
  if max3>maxval
      maxval=max3;
  end
  
  
  
  
  cmap=colormap(jet(256));
  caxis([minval maxval]);
  %caxis([-0.6 0.6]);
  %caxis([4*10^5 3*10^6]);
  divmap=diverging_map(linspace(0,1,256),[0.23 0.299 0.754],[0.706 0.016 0.15]);
  colormap(divmap);
  %colormap(cmap);
  
  
 
  
  
  hc=colorbar();
  set(hc,'Ylim',[minval maxval]);
  %set(hc,'Ylim',[-0.6 0.6]);
  %set(hc,'Ylim',[4*10^5 3*10^6]);
  text(-100,0,165,timetext);
  title('Vertical Velocity for Solar Atmosphere with a Sinusoidal (0,0) Mode Driver of Period673.4s, Applied at a Height of 100km');
  xlabel('x-distance (Mm)');
  ylabel('y-distance (Mm)');
  zlabel('Height (Mm)');
  
  ylabel(hc,'Vz [m/s]');
  
  
  
  print('-djpeg', imfile);
  
  hold off
  
end 
  
