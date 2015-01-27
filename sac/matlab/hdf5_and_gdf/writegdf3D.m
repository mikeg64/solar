function writegdf3D(filename, simparams, simgridinfo, simdata)
    simgridinfo.create_structures_h5(filename);   

    simparams.write_params_h5(filename);
    simgridinfo.write_gridinfo_h5(filename);
    simdata.write_data_h5(filename); 
   
 end  %end of function writegdf3D