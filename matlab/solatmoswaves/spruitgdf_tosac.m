%filename='Y:\Shared\simulations\sacinit\spruit_const.gdf';
filename='/shared/sp2rc2/Shared/simulations/sacinit/spruit_const.gdf';
ofilename='/shared/sp2rc2/Shared/simulations/sacinit/spruit_const_asc.ini';
simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

simgridinfo=simgridinfo.read_gridinfo_h5(filename);
simparams=simparams.read_params_h5(filename);
simdata=simdata.read_data_h5(filename,simparams,simgridinfo);

writesac3D(ofilename,simparams,simgridinfo,simdata,'ascii');