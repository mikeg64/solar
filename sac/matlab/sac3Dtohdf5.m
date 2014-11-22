%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zero1_ot_bin_256.ini';

%sac3Dgetpic
 
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zerospic1__711000.out';
filename='spruit.gdf';
filename='myfile.gdf';

simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

   disp('writing h5 file');
   
   %create goups
   fid = H5F.create(filename);
plist = 'H5P_DEFAULT';
gid = H5G.create(fid,'/simulation_parameters',plist,plist,plist);
H5G.close(gid);

%particle types



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

%grid_dimensions
%
type_id = H5T.copy('H5T_NATIVE_INT64');
dims = [3];
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(1,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(fid,'grid_left_index',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

gid = H5G.create(fid,'/particle_types',plist,plist,plist);
H5G.close(gid);

gid = H5G.create(fid,'/gridded_data_format',plist,plist,plist);


acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,'H5T_VARIABLE');
H5T.set_strpad(type_id,'H5T_STR_NULLTERM');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(gid,'data_software',type_id,space_id,acpl);
%H5A.write(attr_id,'H5ML_DEFAULT',10.0)
H5A.close(attr_id);
H5T.close(type_id);

acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,'H5T_VARIABLE');
H5T.set_strpad(type_id,'H5T_STR_NULLTERM');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(gid,'data_software_version',type_id,space_id,acpl);
%H5A.write(attr_id,'H5ML_DEFAULT',10.0)
H5A.close(attr_id);
H5T.close(type_id);


H5G.close(gid);
   

  
%data
%create data group

%i1= int32  simgridinfo.grid_dimensions(1);
%i2=  (int32) (simgridinfo.grid_dimensions(2));
%i3=  (int32) (simgridinfo.grid_dimensions(3));
dim=simgridinfo.grid_dimensions;
i1=double(dim(1));
i2=double(dim(2));
i3=double(dim(3));
dims = [i1 i2 i3];

plist = 'H5P_DEFAULT';
dgid = H5G.create(fid,'data',plist,plist,plist);
gid = H5G.create(dgid,'grid_0000000000',plist,plist,plist);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'density_pert',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'velocity_x',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'velocity_y',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'velocity_z',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'internal_energy_pert',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_x_pert',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_y_pert',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_z_pert',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'internal_energy_bg',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'density_bg',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_x_bg',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_y_bg',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);

type_id = H5T.copy('H5T_NATIVE_DOUBLE');
h5_dims = fliplr(dims);
h5_maxdims = h5_dims;
space_id = H5S.create_simple(3,h5_dims,h5_maxdims);
dcpl = 'H5P_DEFAULT';
dset_id = H5D.create(gid,'mag_field_z_bg',type_id,space_id,dcpl);
H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);



H5G.close(gid);
H5G.close(dgid);

%field_types group
gid = H5G.create(fid,'/field_types',plist,plist,plist);

sgid = H5G.create(gid,'/field_types/density_bg',plist,plist,plist);

acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,'H5T_VARIABLE');
H5T.set_size(type_id,40);
H5T.set_strpad(type_id,'H5T_STR_NULLTERM');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(sgid,'field_name',type_id,space_id,acpl);
%H5A.write(attr_id,'H5ML_DEFAULT','BackgroundDensity')
H5A.close(attr_id);
H5T.close(type_id);

acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_NATIVE_DOUBLE');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(sgid,'field_to_cgs',type_id,space_id,acpl_id);
H5A.write(attr_id,'H5ML_DEFAULT',0.0010000000002);
H5A.close(attr_id);
H5T.close(type_id);

acpl = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_C_S1');
H5T.set_size(type_id,25);
H5T.set_strpad(type_id,'H5T_STR_NULLTERM');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(sgid,'field_units',type_id,space_id,acpl);
H5A.write(attr_id,'H5ML_DEFAULT','\mathrm{\frac{kg}{m^{3}}}');
H5A.close(attr_id);
H5T.close(type_id);

acpl_id = H5P.create('H5P_ATTRIBUTE_CREATE');
type_id = H5T.copy('H5T_NATIVE_INT64');
space_id = H5S.create('H5S_SCALAR');
attr_id = H5A.create(sgid,'staggering',type_id,space_id,acpl_id);
H5A.write(attr_id,'H5ML_DEFAULT',0);
H5A.close(attr_id);
H5T.close(type_id);


H5G.close(sgid);



H5G.close(gid);





 

H5F.close(fid);
simparams.write_params_h5(filename);
simgridinfo.write_gridinfo_h5(filename);
simdata.write_data_h5(filename); 
   %write hdf5 in gdf format
   %h5filename=[rootfile,'.gdf'];
   %h5create(h5filename, '/dataset1', size(testdata))
   %h5create('test.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5write('template.h5','/data/grid_0000000000/density_bg',[64 64 64]);
   %h5writeatt('test.gdf','/','/simulation_parameters','');
   %simparams.write_params_h5('test.gdf');
   
   
    
