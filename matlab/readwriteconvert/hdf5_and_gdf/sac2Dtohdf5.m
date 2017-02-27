%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zero1_ot_bin_256.ini';

oldencoding=slCharacterEncoding;
slCharacterEncoding('ISO-8859-1');
%rootfile='zero1_ot_bin_256';
rootfile='zeroOT_0';
filename=[rootfile,'.mdl'];
   fid=fopen(trim(filename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=trim(setstr(fread(fid,79,'char')'));
   
 
   it=fread(fid,1,'integer*4'); time=fread(fid,1,'float64');
 
   ndim=fread(fid,1,'integer*4');
   neqpar=fread(fid,1,'integer*4'); 
   nw=fread(fid,1,'integer*4');
   nx=fread(fid,2,'integer*4');
   
   nxs=nx(1)*nx(2);%*nx(3);
   varbuf=fread(fid,7,'float64');
   
   gamma=varbuf(1);
   eta=varbuf(2);
   g(1)=varbuf(3);
   g(2)=varbuf(4);
   %g(3)=varbuf(5);
   
   
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
   %nx3=nx(3);
   
   xx=reshape(X(:,1),nx1,nx2);
   yy=reshape(X(:,2),nx1,nx2);
   %zz=reshape(X(:,3),nx1,nx2,nx3);
   
   
 
  % extract variables from w into variables named after the strings in wnames
wd=zeros(nw,nx1,nx2);
for iw=1:nw
    iw
     tmp=reshape(w(:,iw),nx1,nx2);
     wd(iw,:,:)=tmp;
end


%w=tmp(iw);
  

clear tmp; 
   
   
   fclose(fid);
   
   disp('writing h5 file');
   
   %write hdf5 in gdf format
   h5filename=[rootfile,'.gdf'];
   h5create(h5filename, '/dataset1', size(testdata))
   
   
   
   
    
