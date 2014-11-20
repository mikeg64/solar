%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zero1_ot_bin_256.ini';

%sac3Dgetpic
 
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zerospic1__711000.out';
filename='spruit.gdf';

disp('Reading gdf file ');

simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

t=simparams.read_params_h5(filename);
simparams=t;

t=simgridinfo.read_gridinfo_h5(filename);
simgridinfo=t;

simdata.setsim_params(simparams);
simdata.setsim_gridinfo(simgridinfo);

%data = h5read('spruit.gdf','/data/grid_0000000000/density_bg');
t=simdata.read_data_h5(filename);

   disp('writing h5 file');
   
   %create goups
   fid = H5F.create('myfile.gdf');
plist = 'H5P_DEFAULT';
gid = H5G.create(fid,'/simulation_parameters',plist,plist,plist);
H5G.close(gid);


%grid_dimensions
%
type_id = H5T.copy('H5T_NATIVE_INT64');
dims = [3];
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'grid_dimensions',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
%H5F.close(fid);
%h5disp('myfile.h5');


gid = H5G.create(fid,'/particle_types',plist,plist,plist);
H5G.close(gid);

gid = H5G.create(fid,'/gridded_data_format',plist,plist,plist);
H5G.close(gid);
   
   
 
H5F.close(fid);
 simparams.write_params_h5('myfile.gdf');  
   %write hdf5 in gdf format
   %h5filename=[rootfile,'.gdf'];
   %h5create(h5filename, '/dataset1', size(testdata))
   %h5create('test.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5write('template.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5writeatt('test.gdf','/','/simulation_parameters','');
   %simparams.write_params_h5('test.gdf');
   
   
    
