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

t=simparams.read_params_h5(filename);
simparams=t;

data = h5read('spruit.gdf','/data/grid_0000000000/density_bg');

   disp('writing h5 file');
   
   %write hdf5 in gdf format
   %h5filename=[rootfile,'.gdf'];
   %h5create(h5filename, '/dataset1', size(testdata))
   
   
   
   
    
