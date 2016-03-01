


%% Hold all the constants in consts structure

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7

simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

%% Import the data
data = xlsread('atmos.xls','VALMc_rho_2048_test');







%% Allocate imported array to column variable names
height = data(:,1);
temp = data(:,2);
dens = data(:,3);
pres = data(:,4);

cs=sqrt(consts.fgamma.*pres./dens);





%% Clear temporary variables
clearvars data raw;


%cases
%generate a configuration
%1. uniform horizontal
%2. uniform vertical
%3. inclined vertical
%4. fluxtube

% read a gdf file
% read an ascii file


