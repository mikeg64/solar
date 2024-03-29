        

    filename=newfilename;
    
          fid=fopen(filename, 'w');
        
        
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
           disp(varnames);
           fwrite(fid,varnames,'char*1');
           
           
           
           
           
           
           
  
               is=1;
    js=1;
    ks=1;
    iif=simparams.domain_dimensions(1);
    jf=simparams.domain_dimensions(2);
    kf=simparams.domain_dimensions(3);
    
    
       p.dx(1)=(simparams.domain_right_edge(1)-simparams.domain_left_edge(1))/(simparams.domain_dimensions(1)-1);
    p.dx(2)=(simparams.domain_right_edge(2)-simparams.domain_left_edge(2))/(simparams.domain_dimensions(2)-1);
    p.dx(3)=(simparams.domain_right_edge(3)-simparams.domain_left_edge(3))/(simparams.domain_dimensions(3)-1);
 
 
    
    
    
           for k1=ks:kf
               for j1=js:jf
                     for i1=is:iif
                             for idim=1:ndim
                                 if idim==1
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+(i1-1)*p.dx(idim);
                                 elseif idim==2
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+(j1-1)*p.dx(idim);     
                                 elseif idim==3
                                    X(i1,j1,k1,idim)=simparams.domain_left_edge(idim)+(k1-1)*p.dx(idim);         
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
