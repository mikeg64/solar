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
       		real, dimension(nx1) :: x
       		real, dimension(nx2) :: y
       		real, dimension(nx3) :: z, b0z, dbz
       		real, dimension(nx1,nx2,nx3) :: fx

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
            real :: b0zz = 0.001d0

            real :: dx,dy,dz
            real :: bnmin, bnmax,f, fres, tmp, xfmax, tmpx, tmpy, tmpz, sqmumag
            integer :: i,j,k

!Generate magnetic field configuration
!simple fluxtube using self similarity and hydrostatic pressure correction

            fx=0.0d0

            dx=(xmax-xmin)/nx1
            dy=(ymax-ymin)/nx2
            dz=(zmax-zmin)/nx3

            do i=1,nx1
                x(i)=xmin+dx*(i-1)
            end do

            do j=1,nx2
                y(j)=ymin+dy*(j-1)
            end do


            do k=1,nx3
                z(k)=zmin+dz*(k-1)
            end do

            do k=1,nx3
                b0z(k)=par4((z(k)/sscale-z_shift),d_z,A);
            end do

            bnmin=minval(b0z)
            bnmax=maxval(b0z)

            do i=1,nx3
                b0z(i)=((Bmax-Bmin)/(bnmax-bnmin))*(b0z(i)-bnmin)+Bmin
            end do

    b0z=b0z/maxval(b0z)
    b0z=Ab0z*b0z+b0z_top

    dbz=deriv1(b0z,z)


!Generate field

    do i=1,nx1
        do j=1,nx2
            do k=1,nx3
                tmp=b0z(k)*sqrt(x(i)*x(i)+y(j)*y(j))
                tmp=par4(tmp,f0,A)
                fres=tmp*tmp
                fx(i,j,k) = fres

!                %f=b0z(i)*sqrt((y(j)).^2+(z(k)).^2);
!                %xf(i,j,k)=(par4(f,f0,0.5)).^2;

!                f=(y(j).^2+z(k).^2)./R2;
!                xf(i,j,k)=exp(-f);
            end do
        end do
    end do

    xfmax=maxval(fx)
    fx=fx/xfmax

    sqmumag=sqrt(mumag)

    do i=1,nx1
        do j=1,nx2
            do k=1,nx3
!    %use this one for the multiple tubes
!% bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
!% bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j)-ybp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
!% by(i,j,k)=by(i,j,k)-dbz(i)*(y(k)-zbp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k);
                tmp=sqrt(x(i)*x(i)+y(j)*y(j))

                tmpx=b0z(k)/tmp
                tmpy=dbz(k)*x(i)/tmp
                tmpz=dbz(k)*y(j)/tmp

                ssimdata%w(1,j,k,14)=tmpx*fx(i,j,k)/sqmumag
                ssimdata%w(1,j,k,15)=tmpy*fx(i,j,k)/sqmumag
                ssimdata%w(1,j,k,16)=tmpz*fx(i,j,k)/sqmumag



!%use this one for the tube
!%bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
!%bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
!%by(i,j,k)=by(i,j,k)-dbz(i)*(y(k))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k);


!bz(i,j,k)=bz(i,j,k)+b0zz.*xf(i,j,k);
!bx(i,j,k)=bx(i,j,k);
!by(i,j,k)=by(i,j,k);

!  convert to sac units




            end do
        end do
    end do


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

       	function deriv1(f,x)

            real, intent(in), dimension(:) :: f,x
            real, dimension(size(f)) :: dres
            real, dimension(size(f)) :: deriv1
            integer :: i,nel

            nel=size(f)
            dres=0
            deriv1=0
 !           res(i)=(1./12.0./(x(i+1)-x(i)))*(8.0*f1(i+1)-8.0*f1(i-1)-f1(i+2)+f1(i-2));
            do i=3,nel-2
                dres(i)=(1.0/12.0/(x(i+1)-x(i))) *(8.0*f(i+1)-8.0*f(i-1)-f(i+2)+f(i-2))
            end do

            dres(1)=dres(3)
            dres(2)=dres(3)
            dres(nel)=dres(nel-2)
            dres(nel-1)=dres(nel-2)

            deriv1 = dres

        end function

        real function par4(x,x0,A)
            real, intent(in) :: x, x0, A
            !  [res]=par4(x,x0,A)
            par4=A*exp(-x/(x0))

        end function

        real function par5(x,x0, minx ,A)
            real, intent(in) :: x, x0, minx
            real, intent(inout) :: A
            !res=A.*exp(-1.0*((x+abs(min(x)))/x0)^0.95);
            par5=A*exp(-1.d0*((x+abs(minx))/x0)**0.95)
        end function

        real function par3(x,x0,A)
            real, intent(in) :: x, x0
            real, intent(inout) :: A
            real :: res
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
