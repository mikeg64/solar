%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zero1_ot_bin_256.ini';

%sac3Dgetpic
 
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zerospic1__711000.out';
%filename='spruit.gdf';


gdffilename='myfile.gdf';
sacfilename='zerospic1__711000.out';
simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;






%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';

   fid=fopen((sacfilename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=(setstr(fread(fid,79,'char')'));
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
   
   
   varnames=(setstr(fread(fid,79,'char')'));
   
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   
   xx=reshape(X(:,1),nx1,nx2,nx3);
   yy=reshape(X(:,2),nx1,nx2,nx3);
   zz=reshape(X(:,3),nx1,nx2,nx3);
   
   for iw=1:nw
      %fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      %fread(fid,4);
   end
   
   % extract variables from w into variables named after the strings in wnames
   wd=zeros(nw,nx1,nx2,nx3);
   for iw=1:nw
         tmp=reshape(w(:,iw),nx1,nx2,nx3);
         wd(iw,:,:,:)=tmp;
    end
    clear tmp;
   
   nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);


clear tmp; 
fclose(fid);


%write the sac data into the created structures
simparams.unique_identifier=headline;
simparams.current_iteration=it;         
simparams.current_time = time;         
simparams.dimensionality = ndim;                  
simparams.neqpar=neqpar;          
simparams.nw=nw;
simparams.domain_dimensions=nx;

simparams.gamma=gamma;
simparams.eta=eta;
simparams.gravity0=g1;
simparams.gravity1=g2;
simparams.gravity2=g3;
           
simparams.domain_left_edge(1)=xx(1,1,1);
simparams.domain_left_edge(2)=yy(1,1,1);
simparams.domain_left_edge(3)=zz(1,1,1);

simparams.domain_right_edge(1)=xx(nx1,1,1);
simparams.domain_right_edge(1)=yy(1,nx2,1); 
simparams.domain_right_edge(1)=zz(1,1,nx3);

nw=13;
for iw=1:nw
  simdata.w(:,:,:,iw)=wd(iw,:,:,:);
end




disp('writing h5 file');
writegdf3D(gdffilename,simparams,simgridinfo,simdata);











   
% simgridinfo.create_structures_h5(filename);   
% 
% simparams.write_params_h5(filename);
% simgridinfo.write_gridinfo_h5(filename);
% simdata.write_data_h5(filename); 
   %write hdf5 in gdf format
   %h5filename=[rootfile,'.gdf'];
   %h5create(h5filename, '/dataset1', size(testdata))
   %h5create('test.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5write('template.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5writeatt('test.gdf','/','/simulation_parameters','');
   %simparams.write_params_h5('test.gdf');
   
   
    
