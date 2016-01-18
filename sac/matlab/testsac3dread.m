%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%http://uk.mathworks.com/matlabcentral/answers/97118-how-do-i-read-a-fortran-unformatted-binary-data-file-into-matlab
filename='/fastdata/cs1mkg/sac/p5b0_0_bv10g_a1250/3D_spic_128_bv10g_12.out';

    fid=fopen(filename,'rb');

   hr1=fread(fid, 1, 'int32');
   headline=setstr(fread(fid,79,'uchar')');
   hr2=fread(fid, 1, 'int32'); 
   
    hr1=fread(fid, 1, 'int32');
    it=fread(fid,1,'int32'); 

    time=fread(fid,1,'float64');

    
    %hr1=fread(fid, 1, 'int32');
        ndim=fread(fid,1,'int32');
      %  hr2=fread(fid, 1, 'int32');
    neqpar=fread(fid,1,'int32');
    
 nw=fread(fid,1,'int32');
 hr1=fread(fid, 1, 'int32');
 hr2=fread(fid, 1, 'int32');
 nx=fread(fid,3,'int32');
 
    nxs=nx(1)*nx(2)*nx(3);
     hr1=fread(fid, 1, 'int32');
 hr2=fread(fid, 1, 'int32');
   varbuf=fread(fid,7,'float64');
   
   
      gamma=varbuf(1);
   eta=varbuf(2);
   g(1)=varbuf(3);
   g(2)=varbuf(4);
   g(3)=varbuf(5);
   
   hr1=fread(fid, 1, 'int32');
   varnames=setstr(fread(fid,79,'uchar')');
   hr2=fread(fid, 1, 'int32');
   
   X=zeros(nxs,ndim);
   tw=zeros(nxs,13);
       hr1=fread(fid, 1, 'int32');
     hr2=fread(fid, 1, 'int32');
   for idim=1:ndim
  
      X(:,idim)=fread(fid,nxs,'float64');
      
   end
   
      nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);
   
   xx=reshape(X(:,1),nx1,nx2,nx3);
   yy=reshape(X(:,2),nx1,nx2,nx3);
   zz=reshape(X(:,3),nx1,nx2,nx3);
   
 
for iw1=1:13
      
      tbuf=fread(fid,nxs,'float64');
       stb=size(tbuf);
       if stb(1)<nxs
         for itb=1:nxs-stb(1)
           tbuf(stb(1)+itb)=0.0;
        end
       end
      tw(:,iw1)=tbuf;
end 
   
     % extract variables from w into variables named after the strings in wnames
wd=zeros(nw,nx1,nx2,nx3);
for iw=1:nw
  
     tmp=reshape(tw(:,iw),nx1,nx2,nx3);
     wd(iw,:,:,:)=tmp;
end



  

clear tmp; 
   
   
   fclose(fid);
   val=reshape(wd(1,:,64,:),128,128);
   surf(val,'LineStyle','none');
   colorbar;
   
   val1=reshape(wd(3,:,64,:),128,128)./(reshape(wd(1,:,64,:),128,128)+reshape(wd(10,:,1,:),128,128));
   val2=reshape(wd(4,:,64,:),128,128)./(reshape(wd(1,:,64,:),128,128)+reshape(wd(10,:,1,:),128,128));
   
