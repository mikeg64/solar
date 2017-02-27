path(path, './hdf5_and_gdf/')

%h5filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';
gdffilename='spruit.gdf';
sac3dfilename='spruittest1.out'

% gdfinfo=h5info(h5filename);
% 
% cosmological_simulation=gdfinfo.Groups(5).Attributes(1).Value;
% boundary_conditions=gdfinfo.Groups(5).Attributes(2).Value;
% current_iteration=gdfinfo.Groups(5).Attributes(3).Value;
% current_time=gdfinfo.Groups(5).Attributes(4).Value;
% dimensionality=gdfinfo.Groups(5).Attributes(5).Value;
% domain_dimensions=gdfinfo.Groups(5).Attributes(6).Value;
% domain_left_edge=gdfinfo.Groups(5).Attributes(7).Value;
% domain_right_edge=gdfinfo.Groups(5).Attributes(8).Value;
% eta=gdfinfo.Groups(5).Attributes(9).Value;
% field_ordering=gdfinfo.Groups(5).Attributes(10).Value;
% gamma=gdfinfo.Groups(5).Attributes(11).Value;
% gravity0=gdfinfo.Groups(5).Attributes(12).Value;
% gravity1=gdfinfo.Groups(5).Attributes(13).Value;
% gravity2=gdfinfo.Groups(5).Attributes(14).Value;
% nu=gdfinfo.Groups(5).Attributes(15).Value;
% num_ghost_zones=gdfinfo.Groups(5).Attributes(16).Value;
% refine_by=gdfinfo.Groups(5).Attributes(17).Value;
% unique_identifier=gdfinfo.Groups(5).Attributes(18).Value;

disp('Reading gdf file ');

simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

simparams=simparams.read_params_h5(gdffilename);
%simparams=t;

simgridinfo=simgridinfo.read_gridinfo_h5(gdffilename);
%simgridinfo=t;

%simdata.setsim_params(simparams);
%simdata.setsim_gridinfo(simgridinfo);

%data = h5read('spruit.gdf','/data/grid_0000000000/density_bg');
simdata=simdata.read_data_h5(gdffilename, simparams, simgridinfo);

writesac3D(sac3dfilename,simparams,simgridinfo,simdata,'binary');

