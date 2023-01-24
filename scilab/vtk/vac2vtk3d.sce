function []=vacscalar2vtk3d( pict,wd1,wd2,wd3,xx,yy,zz,field,vecsize,filename,directory )
// pict    :: the pict reference
// wd :: the full array matrix of scalar data
// xx,yy,zz :: the position data     
// field :: which field
//         ; e.g. 1,2,3
// vecsize :: 1,2,3 how many components field has
//         ; e.g. magnetic field, velocity or momentum
// filename :: is a string of the name of the output file (without.vtk)


        
    
    
     sizexx=size(xx);
     sizeyy=size(yy);
     sizezz=size(zz);
     sizewd=size(wd1);

     i=pict;

     if i<9
        filen=msprintf('%s/vtk/%s00%d.vtk',directory,filename,i);
     else
        if i < 99
            filen=msprintf('%s/vtk/%s0%d.vtk',directory,filename,i);         
        else
            filen=msprintf('%s/vtk/%s%d.vtk',directory,filename,i);
        end
     end
      disp(filen);
     fid=mopen(filen, 'w');
     
     disp(sizewd)
     //Header
    mfprintf(fid,'# vtk DataFile Version 2.0\n');
    mfprintf(fid,'Structured Grid\n');
    mfprintf(fid,'ASCII\n');
    mfprintf(fid,'DATASET RECTILINEAR_GRID\n');
    mfprintf(fid,'DIMENSIONS %d %d %d\n',sizewd(1),sizewd(2),sizewd(3));

 
    mfprintf(fid,'X_COORDINATES  %d double \n',sizexx(1));
    for i=1:sizexx(1)
        mfprintf(fid, ' %g ',xx(i));
    end
    mfprintf(fid,'\n');
    
    mfprintf(fid,'Y_COORDINATES  %d double \n',sizeyy(1));
    for i=1:sizeyy(1)
        mfprintf(fid, ' %g ',yy(i));
    end
    mfprintf(fid,'\n');
      
    mfprintf(fid,'Z_COORDINATES  %d double \n',sizezz(1));
    for i=1:sizezz(1)
        mfprintf(fid, ' %g ',zz(i));
    end
    mfprintf(fid,'\n');


      //  printf,lu,'POINT_DATA ',sizew(1)*sizew(2)*sizew(3)
      //  printf,lu,'VECTORS ',filename,' double'


    mfprintf(fid,'POINT_DATA %d \n',sizexx(1)*sizeyy(1)*sizezz(1));
    mfprintf(fid,'VECTORS %s double \n',filename);


//
//        for iz=0,sizew(3)-1 do begin
//            for iy=0,sizew(2)-1 do begin 
//                 for ix=0,sizew(1)-1 do begin
//                    printf,lu,vacdata(ix,iy,iz,field),' ',vacdata(ix,iy,iz,field+1),' ',vacdata(ix,iy,iz,field+2)
//                 endfor
//            endfor
//        endfor



      for i=1:sizexx(1)
              for j=1:sizeyy(1)
                      for k=1:sizezz(1)
                          mfprintf(fid,'%g %g %g',wd1(i,j,k),wd2(i,j,k),wd3(i,j,k));                                                     
                      end
                      mfprintf(fid,'\n');
              end
       end
     
    mclose(fid);
     



endfunction

