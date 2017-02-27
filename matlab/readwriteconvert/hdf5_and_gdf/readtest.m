% filename='Y:\Shared\configs\3D_128_2p5_2p5_12p5_spruit_1kg_asc.ini';

newfilename='Y:\Shared\configs\3D_128_2p5_2p5_12p5_asc.ini';
% simparams=sim_params;
% simgridinfo=sim_gridinfo;
% simdata=sim_data;
% 
% [simparams, simgridinfo, simdata]=readsac3D(filename, simparams, simgridinfo, simdata, 'ascii');


%re-order the data


    is=1;
    js=1;
    ks=1;
    iif=simparams.domain_dimensions(1);
    jf=simparams.domain_dimensions(2);
    kf=simparams.domain_dimensions(3);

    nw=13;

    newsimdata=simdata;

           for k1=ks:kf
               k1
               for j1=js:jf
                     for i1=is:iif
                             for field=1:nw
                                newsimdata.w(i1,j1,k1,field)=simdata.w(i1,128,128,field);
                                
                                newsimdata.w(i1,j1,k1,6)=0;
                                newsimdata.w(i1,j1,k1,7)=0;
                                newsimdata.w(i1,j1,k1,8)=0;

                                newsimdata.w(i1,j1,k1,11)=0;
                                newsimdata.w(i1,j1,k1,12)=0;
                                newsimdata.w(i1,j1,k1,13)=0;
                             end %loop over field    
                       end %i1
                 end %j1
             end %k1
             
             disp('writing');
             writesac3D(newfilename, simparams, simgridinfo, newsimdata, 'ascii');