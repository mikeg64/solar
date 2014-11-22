function writesac3D(filename, simparams, simgridinfo, simdata, mode)

    fid=fopen('filename');
    if mode = 'binary'


    end

    if mode = 'ascii'
        fprintf(fid,'%s\n',simparams.unique_identifier);
        fprintf(fid,'%d %f %d %d\n',simparams.current_iteration, simparams.current_time, simdimensions.dimensionality, simdimensions.nw);
        
        gd=simgridinfo.grid_dimensions;
        fprintf(fid,'%d %d %d\n', gd(1), gd(2), gd(3));
        fprintf(fid,'%f %f %f %f %f %d %d\n',simparams.gamma,simparams.eta,simparams.gravity1,simparams.gravity2,simparams.gravity3,0,0);
        fprintf(fid,'x y z rho mx my mz e bx by bz gamma eta g1 g2 g3\n');

    end
    fclose(fid);


end