module fieldbuildermod
! subroutines used to generate magnetic field configurations


!   ifieldtype1 is self similar field
!   http://solarwavetheory.blogspot.com/2013/11/solar-atmospheric-mhd-flux-tube.html
    use typesmod

    implicit none

    include 'param.inc'

    public :: genssfield, genvecfield


private


contains

        !generate a field using the self similarity model
       	![simparams, simgridinfo, simdata]=
       	!https://ui.adsabs.harvard.edu/abs/2005A%26A...441..337S/abstract
       	subroutine genssfield(ssimparams, ssimgridinfo, ssimdata, fieldtype)



            integer, intent(in) :: fieldtype
       		type(simparams), intent(inout) :: ssimparams
       		type(simgridinfo), intent(inout) :: ssimgridinfo
       		type(simdata), intent(inout)  :: ssimdata

       		!nx1,nx2,nx3 in param.inc
       		real, dimension nx1 :: x
       		real, dimension nx2 :: y
       		real, dimension nx3 :: z

            real :: Bmax = 0.15  !mag field Tesla
            !Bmin=0.0006d0  ; %mag field Tesla
            real :: Bmin = 0.0002  !mag field Tesla
            real :: d_z = 1.5 !width of Gaussian in Mm
            real :: z_shift = 0.0 !shift in Mm
            real :: A = 0.45 !amplitude
            real :: sscale = 1.0e6
            real :: b0z_top = 0.08
            real :: f0 = 2.0d6 !tube opening factor
            real :: Ab0z = 20.d0 !bz - amplitude
            real :: xr = 0.1d6
            real :: yr = 0.1d6

!Generate magnetic field configuration
!simple fluxtube using self similarity and hydrostatic pressure correction



!Generate field
! calculate hydrostatic pressure update


!cases
!1.  Vertical field
!2. Horizontal field
!3. Inclined (off vertical field)
!4. Flux tube array



       	end subroutine genssfield

        !generate a field using magnetic vector potential
       	![simparams, simgridinfo, simdata]=
       	subroutine genvecfield(ssimparams, ssimgridinfo, ssimdata, fieldtype)



            integer, intent(in) :: fieldtype
       		type(simparams), intent(inout) :: ssimparams
       		type(simgridinfo), intent(inout) :: ssimgridinfo
       		type(simdata), intent(inout)  :: ssimdata





       	end subroutine genvecfield

        real function par4(x,x0,A)
            real, intent(in) :: x, x0, A
            !  [res]=par4(x,x0,A)
            par4=A*exp(-x/(x0))

        end function

        real function par5(x,x0,A)
            real, intent(in) :: x, x0, A
            !res=A.*exp(-1.0*((x+abs(min(x)))/x0)^0.95);
            par5=A*exp(-1.d0*((x+abs(min(x)))/x0)**0.95)
        end function

        function par3(x,x0,A)
            real, intent(in) :: x, x0, A
            real, res
            A=A/(x0**2)

            if (x <= -2*x0) then
                res=0
            else if ((x >= -2*x0) .and. (x <= -x0))  then
                res=A*(x+3*x0)*(x+x0)+A*x0**2
            else if ((x >= -x0) .and. (x <= x0)) then
                res=-A*(x+x0)*(x-x0)+A*x0**2
            else if ((x >= x0) .and. (x <= 2*x0)) then
                res=A*(x-3*x0)*(x-x0)+A*x0**2
            else if (x >= 2*x0) then
                res=0;
            end if
            par3=res

        end function


endmodule
