function []=vacscalar2vtk3d( pict,wd,xx,yy,zz,field,vecsize,filename,directory )
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
     sizewd=size(wd);

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



//        printf,lu,'POINT_DATA ',sizew(1)*sizew(2)*sizew(3)
//        printf,lu,'SCALARS ',filename,' double 1'
//
    mfprintf(fid,'POINT_DATA %d \n',sizexx(1)*sizeyy(1)*sizezz(1));
    mfprintf(fid,'SCALARS %s double 1\n',filename);











//        printf,lu,'LOOKUP_TABLE TableName '
//        for iz=0,sizew(3)-1 do begin
//           for iy=0,sizew(2)-1 do begin
//              for ix=0,sizew(1)-1 do begin
//                printf,lu,vacdata(ix,iy,iz,field)
//              endfor
//           endfor
//        endfor



      
     mfprintf(fid,'LOOKUP_TABLE TableName \n'); 
      for i=1:sizexx(1)
              for j=1:sizeyy(1)
                      for k=1:sizezz(1)
                          mfprintf(fid,'%g ',wd(i,j,k));                                                     
                      end
                      mfprintf(fid,'\n');

              end
       end
     
    mclose(fid);
     



endfunction

