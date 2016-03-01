
%filename=('../../sac/matlab/configs/3D_tubeact_128_128_128_asc_50.ini');
filename='../../sac/matlab/configs/test_asc.ini';
simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

[simparams,simgridinfo,simdata]=readsac3D(filename, simparams, simgridinfo, simdata, 'ascii');