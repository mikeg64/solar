function writesac3D(filename, simparams, simgridinfo, simdata, mode)

if mode='ascii'
    fid=fopen(filename, 'r');
    
    
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

        
           headline=char(zeros(1,79));
           sz=size(simparams.unique_identifier);
           if(sz(2)>79)
              tstr=simparams.unique_identifier;
              headline=tstr(1:79);
           else
              headline(1:sz(2))=simparams.unique_identifier; 
           end           
           fwrite(fid,headline,'char*1');          
           
           %it
           fwrite(fid,int32(simparams.current_iteration),'integer*4');            
           %time
           fwrite(fid,double(simparams.current_time),'float64');           
           %ndim
           ndim=simparams.dimensionality;
           fwrite(fid,int32(simparams.dimensionality),'integer*4');
           %neqpar
           fwrite(fid,int32(simparams.neqpar),'integer*4');
           %nw
           fwrite(fid,int32(simparams.nw),'integer*4');
           %nx
           fwrite(fid,int32(simparams.domain_dimensions),'integer*4');
           
           
           nx=simparams.domain_dimensions;
           nxs=nx(1)*nx(2)*nx(3);
           
           %varbuf=fread(fid,7,'float64');
           %gamma=varbuf(1);
           fwrite(fid,double(simparams.gamma),'float64');           
           %eta=varbuf(2);
           fwrite(fid,double(simparams.eta),'float64');
           %g(1)=varbuf(3);
           fwrite(fid,double(simparams.gravity0),'float64');
           %g(2)=varbuf(4);
           fwrite(fid,double(simparams.gravity1),'float64');
           %g(3)=varbuf(5);
           fwrite(fid,double(simparams.gravity2),'float64');


           %varnames=(setstr(fread(fid,79,'char')'));
           tvarnames='x y z rho mx my mz e bx by bz gamma eta g1 g2 g3';
           sz=size(tvarnames);
           varnames=char(zeros(1,79));
           varnames(1:sz(2))=tvarnames;
           fwrite(fid,varnames,'char*1');
           
           for k1=ks:kf
               for j1=js:jf
                     for i1=is:iif
                             for idim=1:ndim
                                 if idim==1
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+i1*p.dx(idim);
                                 elseif idim==2
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+j1*p.dx(idim);     
                                 elseif idim==3
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+k1*p.dx(idim);         
                                 end
                               end %loop over idim      
                       end %i1
                 end %j1
             end %k1
           
           
           for idim=1:ndim
              %X(:,idim)=fread(fid,nxs,'float64');
              fwrite(fid,X(:,idim),'float64');
           end

           nw=13;
           for iw=1:nw
              %fread(fid,4);
              %w(:,iw)=fread(fid,nxs,'float64');
              fwrite(fid,simdata.w(:,:,:,iw),'float64');
              %fread(fid,4);
           end
        
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
                    z=(1+k1)*(p.dx(3));

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
   end %strcmp(mode , 'ascii')





              
    
        
        

 %   end
    fclose(fid);
    display('file writer finished (file closed)');

    end  %end of function writesac3D