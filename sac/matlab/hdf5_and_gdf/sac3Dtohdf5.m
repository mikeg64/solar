%getpicttest  3D version
% Read the npict-th picture from 1 or more files
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zero1_ot_bin_256.ini';

%sac3Dgetpic
 
%filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
%filename='zerospic1__711000.out';
%filename='spruit.gdf';
filename='myfile.gdf';

simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

   disp('writing h5 file');
   
simgridinfo.create_structures_h5(filename);   

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
   
   
    
