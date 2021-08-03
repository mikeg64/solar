function writesac2D(filename, simparams, simgridinfo, simdata, mode)
%lotof effort need to complte this
          fid=fopen(filename, 'w');
%           fdatline=zeros(1,7);
%           fdatline(1)=double(simparams.gamma);
%           fdatline(2)=double(simparams.eta);
%           fdatline(3)=double(simparams.adiab);
%           fdatline(4)=double(simparams.gravity0);
%           fdatline(5)=double(simparams.gravity1);
%           fdatline(6)=double(simparams.gravity2);
%           fdatline(7)=;
   
      p.dx(1)=(simparams.domain_right_edge(1)-simparams.domain_left_edge(1))/(simparams.domain_dimensions(1)-1);
    p.dx(2)=(simparams.domain_right_edge(2)-simparams.domain_left_edge(2))/(simparams.domain_dimensions(2)-1);
  %  p.dx(3)=(simparams.domain_right_edge(3)-simparams.domain_left_edge(3))/(simparams.domain_dimensions(3)-1);



    if strcmp(mode , 'binary')

        
        
        
        
        
           headline=char(zeros(1,79));
           sz=size(simparams.unique_identifier);
           if(sz(2)>79)
              tstr=simparams.unique_identifier;
              headline=tstr(1:79);
           else
              headline(1:sz(2))=simparams.unique_identifier; 
           end 
           
%            fwrite(fid, 1, 'int32');
           fwrite(fid,headline,'char*1');          
%            fwrite(fid, 1, 'int32');
           
           
%            fwrite(fid, 1, 'int32');
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
%            fwrite(fid, 1, 'int32');
           
%            fwrite(fid, 1, 'int32');
           %nx
           fwrite(fid,int32(simparams.domain_dimensions),'integer*4');
%            fwrite(fid, 1, 'int32');
           
           nx=simparams.domain_dimensions;
           nxs=nx(1)*nx(2);
%            fwrite(fid, 1, 'int32');
           
           fwrite(fid, 1, 'int32');
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
           fwrite(fid,double(0.0),'float64');
           fwrite(fid,double(0.0),'float64');
           

           %varnames=(setstr(fread(fid,79,'char')'));
           
           tvarnames='x h m1 m2 e b1 b2 eb rhob bg1 bg2 gamma eta grav1 grav2';
           sz=size(tvarnames);
           varnames=char(zeros(1,79));
           varnames(1:sz(2))=tvarnames;
           disp(varnames);
%            fwrite(fid, 1, 'int32');
           fwrite(fid,varnames,'char*1');
%            fwrite(fid, 1, 'int32');
           
           
           
           
           
           
  
               is=1;
    js=1;
    ks=1;
    iif=simparams.domain_dimensions(1);
    jf=simparams.domain_dimensions(2);
   % kf=simparams.domain_dimensions(3);
    
    
%   fwrite(fid, 1, 'int32');
  
%  fwrite(fid, 1, 'int32');
    
    
    
        %   for k1=ks:kf
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
             %end %k1
           
           nw=10;
           %for k1=ks:kf
               for j1=js:jf
                     for i1=is:iif
                             for idim=1:ndim
                                 
                                 fwrite(fid,X(i1,j1,k1,idim),'float64');
                                 
             
                             end %loop over idim 
                             
                             

                               for iw=1:nw
                                  %fread(fid,4);
                                  %w(:,iw)=fread(fid,nxs,'float64');
                                  fwrite(fid,simdata.w(i1,j1,k1,iw),'float64');
                                  %fread(fid,4);
                               end
                             
                             
                             
                       end %i1
                 end %j1
            %end %k1
             
           fwrite(fid, '\n', 'char*1');
             
     
        
        display('write binary sac file');
    end

   if strcmp(mode , 'ascii')
       
    
    
    is=1;
    js=1;
    ks=1;
    iif=simparams.domain_dimensions(1);
    jf=simparams.domain_dimensions(2);
    %kf=simparams.domain_dimensions(3);
   
   %iif=4;
   %jf=4;
   %kf=4;
    
       
       
       
       
        fprintf(fid,'%s\n',simparams.unique_identifier);
        fprintf(fid,'%d %f %d %d %d\n',simparams.current_iteration, simparams.current_time, simgridinfo.ndimensions,6, 10);
        
        gd=simgridinfo.grid_dimensions;
        fprintf(fid,'%d %d \n', gd(1), gd(2));
        fprintf(fid,'%f %f %f %f %f %d %d\n',simparams.gamma,simparams.eta,simparams.adiab,simparams.gravity0,simparams.gravity1,0,0);
        fprintf(fid,'x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2\n');
        
        

        %for k1=ks:kf

           for j1=js:jf
                 for i1=is:iif

                    x=simparams.domain_left_edge(1)+(i1-1)*p.dx(1);
                    y=simparams.domain_left_edge(2)+(j1-1)*p.dx(2);
                    %z=simparams.domain_left_edge(3)+(k1-1)*p.dx(3);

                    rho=simdata.w(i1,j1,1);
                    mx=simdata.w(i1,j1,2);
                    my=simdata.w(i1,j1,3);
                    %mz=simdata.w(i1,j1,k1,4);
                    e=simdata.w(i1,j1,4);
                    bx=simdata.w(i1,j1,5);
                    by=simdata.w(i1,j1,6);
                    %bz=simdata.w(i1,j1,k1,8);
                    eb=simdata.w(i1,j1,7);
                    rhob=simdata.w(i1,j1,8);
                    b1b=simdata.w(i1,j1,9);
                    b2b=simdata.w(i1,j1,10);
                    %b3b=simdata.w(i1,j1,k1,13);


                                 %shift=(k1*ni*nj+j1*ni+i1);
                                 fprintf(fid,'    %G %G %G %G %G %G %G %G %G %G %G %G\n',x,y, rho,mx,my,e,bx,by,eb,rhob,b1b,b2b);



                 end  %loop over i values
           end %loop over j values
        %end  %loop over k values
        display('write ascii sac file');
   end %strcmp(mode , 'ascii')





              
    
        
        

 %   end
    fclose(fid);
    display('file writer finished (file closed)');

    end  %end of function writesac3D
