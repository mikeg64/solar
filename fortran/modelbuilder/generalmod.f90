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
            real :: height, logdens, logtt
!Use
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/hdf5_and_gdf/writesac3D.m
              ! Open file for unformmated output using the inquired record length
              open(unit = 1, file = newfilename, status='new', action='write')

                !write the header lines
                  print *, 'write out'
                  print *, height, logdens, logtt
                  write(1,'(a100)') ssimparams%uniqueidentifier

                  !it
                  write(1,'(i6)') ssimparams%currentiteration
                  !write(1,'(f20.14,7X, f20.14, 7X, f20.14)') height, logdens, logtt



                !write the fields for each location



              close(unit = 1)

end subroutine writesac3d




endmodule
