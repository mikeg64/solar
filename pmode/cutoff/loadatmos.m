


%% Hold all the constants in consts structure

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7



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




%load('atmos.xls');
h1=height(1420:2048);
h2=height(1324:1419);
h3=height(1:1325);
t1=temp(1420:2048);
t2=temp(1324:1419);
t3=temp(1:1325);
d1=dens(1420:2048);
d2=dens(1324:1419);
d3=dens(1:1325);

cs1=cs(1420:2048);
cs2=cs(1324:1419);
cs3=cs(1:1325);


ra1=0.006757, rb1=-0.2982, rc1=-9.158e-5;
pa1=7.515e5, pb1=-0.3386, pc1=-5681;
%pa1, pb1, pc1, r21, r22, r23, r24, r25, p21, p22, p23, p24, p25, ra3, rb3, rc3, pa3, pb3, pc3;

rhofit1=ra1*(h1.^rb1)+rc1;
pfit1=pa1*(h1.^pb1)+pc1;

