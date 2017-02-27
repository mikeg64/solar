function [simparams, simgridinfo, simdata]=readsac3D(filename, simparams, simgridinfo, simdata, mode)

    fid=fopen(filename, 'r');
    
    
     
    if strcmp(mode , 'binary')

    
        
       hr1=fread(fid, 1, 'int32');
       headline=setstr(fread(fid,79,'uchar')');
       hr2=fread(fid, 1, 'int32'); 

       hr1=fread(fid, 1, 'int32');
       it=fread(fid,1,'int32'); 

       time=fread(fid,1,'float64');


        %hr1=fread(fid, 1, 'int32');
        ndim=fread(fid,1,'int32');
        %  hr2=fread(fid, 1, 'int32');
        neqpar=fread(fid,1,'int32');

         nw=fread(fid,1,'int32');
         hr1=fread(fid, 1, 'int32');
         hr2=fread(fid, 1, 'int32');
         nx=fread(fid,3,'int32');

         nxs=nx(1)*nx(2)*nx(3);
         hr1=fread(fid, 1, 'int32');
         hr2=fread(fid, 1, 'int32');
         varbuf=fread(fid,7,'float64');


       gamma=varbuf(1);
       eta=varbuf(2);
       g(1)=varbuf(3);
       g(2)=varbuf(4);
       g(3)=varbuf(5);

       hr1=fread(fid, 1, 'int32');
       varnames=setstr(fread(fid,79,'uchar')');
       hr2=fread(fid, 1, 'int32');

       X=zeros(nxs,ndim);
       tw=zeros(nxs,13);
           hr1=fread(fid, 1, 'int32');
         hr2=fread(fid, 1, 'int32');
       for idim=1:ndim

          X(:,idim)=fread(fid,nxs,'float64');

       end

       nx1=nx(1);
       nx2=nx(2);
       nx3=nx(3);

       xx=reshape(X(:,1),nx1,nx2,nx3);
       yy=reshape(X(:,2),nx1,nx2,nx3);
       zz=reshape(X(:,3),nx1,nx2,nx3);


        for iw1=1:13
              tbuf=fread(fid,nxs,'float64');
               stb=size(tbuf);
               if stb(1)<nxs
                 for itb=1:nxs-stb(1)
                   tbuf(stb(1)+itb)=0.0;
                end
               end
              tw(:,iw1)=tbuf;
        end 

        % extract variables from w into variables named after the strings in wnames
        wd=zeros(nw,nx1,nx2,nx3);
        for iw=1:nw

             tmp=reshape(tw(:,iw),nx1,nx2,nx3);
             wd(iw,:,:,:)=tmp;
        end
        clear tmp;
        
        
        
    end
    
    
   if strcmp(mode , 'ascii')
       
         
       headline=fgets(fid);
       buf1=fgets(fid);
       varbuf=sscanf(buf1,'%f');
       it=int32(varbuf(1));
       time=int32(varbuf(2));
       ndim=int32(varbuf(3));
       neqpar=int32(varbuf(4));
       %nw=int32(varbuf(5));
       
       if ndim == 3
           nw=13;
       else
           nw=10;
       end
       buf2=fgets(fid);
       varbuf=sscanf(buf2,'%f');
       nx1=int32(varbuf(1));
       nx2=int32(varbuf(2));
       nx3=int32(varbuf(3));
       
       buf3=fgets(fid);
       varbuf=sscanf(buf3,'%f');
       
       gamma=varbuf(1);
       eta=varbuf(2);
       g1=varbuf(3);
       g2=varbuf(4);
       g3=varbuf(5);
       
       varnames=fgets(fid);
       
       wd=zeros(nw,nx1,nx2,nx3);
       xx=zeros(nx1,nx2,nx3);
       yy=zeros(nx1,nx2,nx3);
       zz=zeros(nx1,nx2,nx3);
       simdata.w=zeros(nx1,nx2,nx3,nw);
       
       
       for i=1:nx1
         for j=1:nx2
            for k=1:nx3
                  buf3=fgets(fid);
                  varbuf=sscanf(buf3,'%f');
                                   
                  xx(i,j,k)=varbuf(1);
                  yy(i,j,k)=varbuf(2);
                  zz(i,j,k)=varbuf(3);

                  for field=1:nw
                      simdata.w(i,j,k,field)=varbuf(3+field);
                  end
           
            end         
          end                  
       end
       
       
   
       display(headline);
       
       
       
       
       
   end
    
   
   
    fclose(fid);
    
    
    
    
        simparams.current_iteration=it;
        simparams.current_time=time;
        simparams.dimensionality=ndim;
        simparams.domain_dimensions=[nx1;nx2;nx3];
%         simparams.domain_left_edge=[0;0; 0.0];
%         simparams.domain_right_edge=[0; 0; 0];
        simparams.eta=eta;
%        field_ordering=1;
        simparams.gamma=gamma;
        simparams.gravity0=g1;
        simparams.gravity1=g2;
        simparams.gravity2=g3;
        
        simgridinfo.grid_dimensions(1)=nx1;
        simgridinfo.grid_dimensions(2)=nx2;
        simgridinfo.grid_dimensions(3)=nx3;
        
        simparams.domain_left_edge(1)=xx(1,1,1);
        simparams.domain_left_edge(2)=yy(1,1,1);
        simparams.domain_left_edge(3)=zz(1,1,1);
        
        simparams.domain_right_edge(1)=xx(nx1,1,1);
        simparams.domain_right_edge(2)=yy(1,nx2,1);
        simparams.domain_right_edge(3)=zz(1,1,nx3);

    
    
    
    display('file writer finished (file closed)');

    %end  %end of function writesac3D