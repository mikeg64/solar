! types and general functions used for the model builder
module generalmod

    use typesmod

    implicit none

    include 'param.inc'
    public :: writesacfile, writesac3d
    private

contains

subroutine writesacfile(height, dens, press, temp,bruntvas, nitems)
    implicit none

! inputs w matrix with fields

! parameter block


    integer, intent(in) :: nitems
    real, intent(in) :: height(nitems), dens(nitems), press(nitems), temp(nitems),bruntvas(nitems)

    integer :: i

    open(unit=9, file='atmos.txt')

    !starting at point number 4 because lower points had -ve enrgy density
    do i=37,nitems
        write(9,*) height(i),temp(i),dens(i), press(i), bruntvas(i)
200     format(F16.6,2X,F16.6,2X,F16.6,2X,F16.6)
    end do

    close(unit=9)

end subroutine writesacfile

!writesac3D(newfilename, simparams, simgridinfo, simdata, 'ascii');
subroutine writesac3d( newfilename,   ssimparams, ssimgridinfo, ssimdata, smconsts )
        implicit none
       		type(simparams), intent(inout) :: ssimparams
       		type(simgridinfo), intent(inout) :: ssimgridinfo
       		type(simdata), intent(inout)  :: ssimdata
      		type(mconsts), intent(inout)  :: smconsts
            character (len=30), intent(in) :: newfilename
            real :: dx1, dx2, dx3
            real :: x,y,z,rho,mx,my,mz
            real :: e,bx,by,bz
            real :: eb, rhob, b1b, b2b,b3b
            real :: height, logdens, logtt
            integer :: i1, j1,k1
!Use
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/hdf5_and_gdf/writesac3D.m
              ! Open file for unformmated output using the inquired record length
              open(unit = 1, file = newfilename, status='new', action='write')

                dx1=(xmax-xmin)/nx1
                dx2=(ymax-ymin)/nx2
                dx3=(zmax-zmin)/nx3


                !write the header lines
                  print *, 'write out'
                  print *, height, logdens, logtt
                  !       fprintf(fid,'%s\n',simparams.unique_identifier);
                  write(1,'(a100)') ssimparams%uniqueidentifier

                  !it
                !fprintf(fid,'%d %f %d %d %d\n',simparams.current_iteration, simparams.current_time, simgridinfo.ndimensions,7, 13);
                  write(1,'(i6, 2X, f20.14, 2X, i6, 2X, i6, 2X, i6)') &
                   ssimparams%currentiteration, ssimparams%currenttime, ssimgridinfo%ndimensions, 7, 13


!        gd=simgridinfo.grid_dimensions;
!        fprintf(fid,'%d %d %d\n', gd(1), gd(2), gd(3));
                  write(1,'(i6, 2X, i6, 2X, i6)') &
                   ssimgridinfo%grid_dimensions(1), ssimgridinfo%grid_dimensions(2), ssimgridinfo%grid_dimensions(3)

!        fprintf(fid,'%f %f %f %f %f %d %d\n',simparams.gamma,simparams.eta,simparams.adiab,simparams.gravity0,simparams.gravity1,simparams.gravity2,0,0);
                  write(1,'(f20.14, 2X, f20.14, 2X, f20.14, 2X, f20.14, 2X, f20.14, 2X, f20.14, 2X, f20.14, 2X, f20.14)') &
                   ssimparams%fgamma, ssimparams%eta, ssimparams%adiab, &
                   ssimparams%gravity0, ssimparams%gravity1, ssimparams%gravity2, 0.0d0, 0.0d0


!        fprintf(fid,'x y z rho mx my mz e bx by bz gamma eta g1 g2 g3\n');
                    write(1,*) 'x y z rho mx my mz e bx by bz gamma eta g1 g2 g3'

                !write the fields for each location
                 do i1=1,nx1
                    do j1=1,nx2
                        do k1=1,nx3

                            x=ssimdata%w(i1,j1,k1,1)
                            y=ssimdata%w(i1,j1,k1,2)
                            z=ssimdata%w(i1,j1,k1,3)

                            rho=ssimdata%w(i1,j1,k1,4)
                            mx=ssimdata%w(i1,j1,k1,5)
                            my=ssimdata%w(i1,j1,k1,6)
                            mz=ssimdata%w(i1,j1,k1,7)
                            e=ssimdata%w(i1,j1,k1,8)
                            bx=ssimdata%w(i1,j1,k1,9)
                            by=ssimdata%w(i1,j1,k1,10)
                            bz=ssimdata%w(i1,j1,k1,11)
                            eb=ssimdata%w(i1,j1,k1,12)
                            rhob=ssimdata%w(i1,j1,k1,13)
                            b1b=ssimdata%w(i1,j1,k1,14)
                            b2b=ssimdata%w(i1,j1,k1,15)
                            b3b=ssimdata%w(i1,j1,k1,16)
                            write(1,'(3(f20.8, 2X),13(f22.12, 2X))') &
                            x,y,z,rho,mx,my,mz,e,bx,by,bz,eb,rhob,b1b,b2b,b3b

                         end do
                    end do
                 end do


              close(unit = 1)

end subroutine writesac3d




endmodule
