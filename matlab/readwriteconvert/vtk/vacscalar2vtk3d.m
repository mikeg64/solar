function [ output_args ] = vacscalar2vtk3d( pict,wd,xx,yy,zz,field,vecsize,vtkfile,directory )
% pict    :: the pict reference
% wd :: the full array of vacdata from getpict 
% xx,yy,zz :: the position data     
% field :: which field
%         ; e.g. 1,2,3
% vecsize :: 1,2,3 how many components field has
%         ; e.g. magnetic field, velocity or momentum
% filename :: is a string of the name of the output file (without.vtk)


        
    
    
     sizexx=size(xx);
     sizeyy=size(yy);
     sizezz=size(zz);
     sizewd=size(wd);

     i=pict;

     if i<9
        filen=[directory,'/vtk/',filename,'00',int2str(i)),'.vtk']
     else
        if i < 99
          filen=[directory,'/vtk/',filename,'0',int2str(i)),'.vtk']
        else
          filen=[directory,'/vtk/',filename,int2str(i)),'.vtk']
        end
     end
       
     fid=fopen(filen, 'w');
     
     %Header
    fprintf(fid,'# vtk DataFile Version 2.0\n');
    fprintf(fid,'Structured Grid\n');
    fprintf(fid,'ASCII\n');
    fprintf(fid,'DATASET RECTILINEAR_GRID\n');
    fprintf(fid,'DIMENSIONS %d %d %d\n',sizewd(1),sizewd(2),sizewd(3));

 

        printf,lu,'X_COORDINATES ',sizew(1),' double'
        for ix=0,sizew(1)-1 do begin
           printf,lu,x(ix,0,0)
        endfor
     
     
     
     
     
     
  
       
       
        fprintf(fid,'%s\n',simparams.unique_identifier);
        fprintf(fid,'%d %f %d %d\n',simparams.current_iteration, simparams.current_time, simgridinfo.ndimensions, 12);
        
        gd=simgridinfo.grid_dimensions;
        fprintf(fid,'%d %d %d\n', gd(1), gd(2), gd(3));
        fprintf(fid,'%f %f %f %f %f %d %d\n',simparams.gamma,simparams.eta,simparams.gravity0,simparams.gravity1,simparams.gravity2,0,0);
        fprintf(fid,'x y z rho mx my mz e bx by bz gamma eta g1 g2 g3\n');
        
        

        for k1=ks:kf

           for j1=js:jf
                 for i1=is:iif

                    x=simparams.domain_left_edge(1)+(i1-1)*p.dx(1);
                    y=simparams.domain_left_edge(2)+(j1-1)*p.dx(2);
                    z=simparams.domain_left_edge(3)+(k1-1)*p.dx(3);

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
        display('write ascii sac file');
 


end

