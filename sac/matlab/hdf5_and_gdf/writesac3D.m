function writesac3D(filename, simparams, simgridinfo, simdata, mode)

    fid=fopen(filename, 'w');
    
    
    is=1;
    js=1;
    ks=1;
    iif=simparams.domain_dimensions(1);
    jf=simparams.domain_dimensions(2);
    kf=simparams.domain_dimensions(3);
   
   %iif=4;
   %jf=4;
   %kf=4;
    
    p.dx(1)=(simparams.domain_right_edge(1)-simparams.domain_left_edge(1))/(simparams.domain_dimensions(1));
    p.dx(2)=(simparams.domain_right_edge(2)-simparams.domain_left_edge(2))/(simparams.domain_dimensions(2));
    p.dx(3)=(simparams.domain_right_edge(3)-simparams.domain_left_edge(3))/(simparams.domain_dimensions(3));
    
    if strcmp(mode , 'binary')

        display('write binary sac file');
    end

    if strcmp(mode , 'ascii')
        fprintf(fid,'%s\n',simparams.unique_identifier);
        fprintf(fid,'%d %f %d %d\n',simparams.current_iteration, simparams.current_time, simgridinfo.ndimensions, 12);
        
        gd=simgridinfo.grid_dimensions;
        fprintf(fid,'%d %d %d\n', gd(1), gd(2), gd(3));
        fprintf(fid,'%f %f %f %f %f %d %d\n',simparams.gamma,simparams.eta,simparams.gravity0,simparams.gravity1,simparams.gravity2,0,0);
        fprintf(fid,'x y z rho mx my mz e bx by bz gamma eta g1 g2 g3\n');
        
        

for k1=ks:kf

   for j1=js:jf
	     for i1=is:iif
             
                         x=(1+i1)*(p.dx(1));
                      y=(1+j1)*(p.dx(2));
			
		    z=(1+k1)*(p.dx(2));
            
			rho=simdata.w(i1,j1,k1,1);
            mx=simdata.w(i1,j1,k1,2);
            my=simdata.w(i1,j1,k1,3);
            mz=simdata.w(i1,j1,k1,4);
            e=simdata.w(i1,j1,k1,5);
            bx=simdata.w(i1,j1,k1,6);
            by=simdata.w(i1,j1,k1,7);
            bz=simdata.w(i1,j1,k1,8);
            eb=simdata.w(i1,j1,k1,9);
            rhob=simdata.w(i1,j1,k1,10);
            b1b=simdata.w(i1,j1,k1,11);
            b2b=simdata.w(i1,j1,k1,12);
            b3b=simdata.w(i1,j1,k1,13);


                         %shift=(k1*ni*nj+j1*ni+i1);
                         fprintf(fid,'%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n',x,y,z, rho,mx,my,mz,e,bx,by,bz,eb,rhob,b1b,b2b,b3b);



         end  %loop over i values
   end %loop over j values
end  %loop over k values





              
    
        
        

    end
    fclose(fid);


end