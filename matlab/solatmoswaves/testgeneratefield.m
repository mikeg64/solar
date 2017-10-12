simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

%newfilename='3D_128_4Mm_tube1_asc.ini';


newfilename='3D_128_4Mm_asc.ini';

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

        
        
        
        
        
        %getpicttest  3D version
% Read the npict-th picture from 1 or more files
%http://uk.mathworks.com/matlabcentral/answers/97118-how-do-i-read-a-fortran-unformatted-binary-data-file-into-matlab

%The config file 3D_128_spic_bin.ini is in the archive
%pmodeini.tgz available at the following link
%https://drive.google.com/open?id=0B-AKVl-pk6ziVUFqUzMycWJoODg

filename='3D_128_4Mm_bin.ini';

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


simdata.w=shiftdim(wd,1);
  

clear tmp; 
   
   
   fclose(fid);

   
   
   disp('generate field');
%
[simparams, simgridinfo, simdata]=generatefield(simparams, simgridinfo, simdata, 'fluxtube');

writesac3D(newfilename, simparams, simgridinfo, simdata, 'binary');

        
        
        
        
        
        