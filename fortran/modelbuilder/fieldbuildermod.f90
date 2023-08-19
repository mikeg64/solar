module fieldbuildermod
! subroutines used to generate magnetic field configurations
! based on original code at
! https://github.com/mikeg64/sac_working/blob/master/mag_field/B_field_vertical_tube.pro

!   ifieldtype1 is self similar field
!   http://solarwavetheory.blogspot.com/2013/11/solar-atmospheric-mhd-flux-tube.html
    use typesmod

    implicit none

    include 'param.inc'

    public :: genssfield, genvecfield, hsbalancefield


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
       		real, dimension(nx3) :: x
       		real, dimension(nx2) :: y
       		real, dimension(nx1) :: z, b0z, dbz
       		real, dimension(nx1,nx2,nx3) :: fx

            real :: Bmax = 0.01  !0.005 !0.15  !mag field Tesla
            !Bmin=0.0006d0  ; %mag field Tesla
            real :: Bmin = 0.001 !0.0002  !mag field Tesla
            real :: d_z = 0.15 !1.5 !width of Gaussian in Mm
            real :: z_shift = 0.0 !shift in Mm
            real :: A = 0.45 !amplitude
            real :: sscale = 1.0e6
            real :: b0z_top = 0.08
            real :: f0 = 2.0d6 !2.0d6tube opening factor
            real :: Ab0z = 20.d0 !bz - amplitude
            real :: xr = 0.15d6
            real :: yr = 0.15d6
            real :: b0zz = 0.001   !0.001d0

            real :: dx,dy,dz
            real :: bnmin, bnmax,f, fres, tmp, xfmax, tmpx, tmpy, tmpz, sqmumag
            integer :: i,j,k

!Generate magnetic field configuration
!simple fluxtube using self similarity and hydrostatic pressure correction
! check against
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/ ...
! ...  generateinitialconfiguration/generatefield_verttube.m
            fx=0.0d0

            dx=(xmax-xmin)/nx1
            dy=(ymax-ymin)/nx2
            dz=(zmax-zmin)/nx3

            do i=1,nx1
                x(i)=xmin+dx*(i-1)
            end do

            do j=1,nx2
                y(j)=ymin+dy*(j-1)-2.0d6
            end do


            do k=1,nx3
                z(k)=zmin+dz*(k-1)-2.0d6
            end do

            do i=1,nx1
                b0z(i)=par4((x(i)/sscale-z_shift),d_z,A)
            end do

            bnmin=minval(b0z)
            bnmax=maxval(b0z)

            do i=1,nx1
                b0z(i)=((Bmax-Bmin)/(bnmax-bnmin))*(b0z(i)-bnmin)+Bmin
            end do

  !  b0z=b0z/maxval(b0z)
  !  b0z=Ab0z*b0z+b0z_top

    dbz=deriv1(b0z,x)


!Generate field

    do i=1,nx1
        do j=1,nx2
            do k=1,nx3
                tmp=b0z(i)*sqrt(z(k)*z(k)+y(j)*y(j))
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
                tmp=sqrt(z(k)*z(k)+y(j)*y(j))

                tmpx=b0z(i)/tmp
                tmpy=dbz(i)*(y(j))/tmp
                tmpz=dbz(i)*(z(k))/tmp

                ssimdata%w(i,j,k,14)=tmpx*fx(i,j,k)/sqmumag
                ssimdata%w(i,j,k,15)=tmpy*fx(i,j,k)/sqmumag
                ssimdata%w(i,j,k,16)=tmpz*fx(i,j,k)/sqmumag



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
        write(*,*) b0z(i), ssimdata%w(i,16,16,14),tmp,fx(i,16,16)
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

!       For magneto-hydrostatic balance
!        see doi:10.1088/0004-637X/727/1/17
!        https://iopscience.iop.org/article/10.1088/0004-637X/727/1/17
       	subroutine hsbalancefield(ssimparams, ssimgridinfo, ssimdata)


       		type(simparams), intent(inout) :: ssimparams
       		type(simgridinfo), intent(inout) :: ssimgridinfo
       		type(simdata), intent(inout)  :: ssimdata

       		!nx1,nx2,nx3 in param.inc
       		real, dimension(nx1) :: x
       		real, dimension(nx2) :: y
       		real, dimension(nx3) :: z, b0z, dbz
       		real, dimension(nx1,nx2,nx3) :: bx, by, bz
       		real, dimension(nx1,nx2,nx3) :: dbzdz, dbydz, dbydx, dbxdx
     		real, dimension(nx1,nx2,nx3) :: dbxdz
     		real, dimension(nx1,nx2,nx3) :: dbxdy
     		real, dimension(nx1,nx2,nx3) :: dbydy
            real, dimension(nx1,nx2,nx3) :: dbzdy
            real, dimension(nx1,nx2,nx3) :: dbzdx
    		real, dimension(nx1,nx2,nx3) :: bvarix, bvariy, bvari, dpdz, rho1, pres
    		real, dimension(nx1,nx2,nx3) :: bvaridz, bvaridz1, dbvardz, br
    		real, dimension(nx1,nx2,nx3) :: bxby, dbxbydy, bxbz, dbxbzdz
    		real, dimension(nx1,nx2,nx3) :: bybz, dbybzdy, dbxbydx, divb, f,g
    		real, dimension(nx1,nx2,nx3) :: dbybzdz, dbxbybzdz, dbxbzdx
    		real, dimension(nx1,nx2,nx3) :: bxbybz
    		integer :: i,j,k
    		real :: dx, dy,dz, sumi, p_2, emin, rhomin, lowval

!   compute magnetostatics pressure correction
            dbzdz=0.0d0
            dbxdz=0.0d0
            dbydz=0.0d0

            dbzdx=0.0d0
            dbxdx=0.0d0
            dbydx=0.0d0

            dbzdy=0.0d0
            dbxdy=0.0d0
            dbydy=0.0d0


            Bvarix=0.0d0
            Bvariy=0.0d0
            Bvari=0.0d0

            dpdz=0.0d0

            pres=0.0d0
            rho1=0.0d0

            Bvaridz=0.0d0
            Bvaridz1=0.0d0
            dBvardz=0.0d0
            br=0.0d0

            bxby=0.0d0
            dbxbydy=0.0d0
            bxbz=0.0d0
            dbxbzdz=0.0d0
            bxby=0.0d0
            bybz=0.0d0
            dbybzdy=0.0d0

            bxbybz=0.0d0

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

!           do i=1,nx1
!                do j=1,nx2
!                    do k=1,nx3
                        bz(1:nx1,1:nx2,1:nx3)=ssimdata%w(1:nx1,1:nx2,1:nx3,9)+ssimdata%w(1:nx1,1:nx2,1:nx3,14)
                        bx(1:nx1,1:nx2,1:nx3)=ssimdata%w(1:nx1,1:nx2,1:nx3,10)+ssimdata%w(1:nx1,1:nx2,1:nx3,15)
                        by(1:nx1,1:nx2,1:nx3)=ssimdata%w(1:nx1,1:nx2,1:nx3,11)+ssimdata%w(1:nx1,1:nx2,1:nx3,16)
!                    end do
!                end do
!           end do


            do k=1,nx3
            do j=1,nx2
                dbzdz(1:nx1,j,k)=deriv1(bz(1:nx1,j,k),x)
                dbxdz(1:nx1,j,k)=deriv1(bx(1:nx1,j,k),x)
                dbydz(1:nx1,j,k)=deriv1(by(1:nx1,j,k),x)
            enddo
            enddo

            do k=1,nx3
            do i=1,nx1
                dbxdx(i,1:nx2,k)=deriv1(bx(i,1:nx2,k),y)
                dbzdx(i,1:nx2,k)=deriv1(bz(i,1:nx2,k),y)
                dbydx(i,1:nx2,k)=deriv1(by(i,1:nx2,k),y)
            enddo
            enddo

            do j=1,nx2
            do i=1,nx1
             dbydy(i,j,1:nx3)=deriv1(by(i,j,1:nx3),z)
             dbzdy(i,j,1:nx3)=deriv1(bz(i,j,1:nx3),z)
             dbxdy(i,j,1:nx3)=deriv1(bx(i,j,1:nx3),z)
            enddo
            enddo

            divb=dbzdz+dbxdx+dbydy


!           %check b_field_vertical_tube.pro
            bxby=bx*by
            do j=1,nx2
            do i=1,nx1
             dbxbydy(i,j,1:nx3)=deriv1(bxby(i,j,1:nx3),z)
            enddo
            enddo


!        %   print,'dBxBydy'
            bxbz=bx*bz
            do k=1,nx3
            do j=1,nx2
             dbxbzdz(1:nx1,j,k)=deriv1(bxbz(1:nx1,j,k),x)
            enddo
            enddo

!            %;print,'dBxBzdz'
            bxby=bx*by
            do i=1,nx1
            do k=1,nx3
             dbxbydx(i,1:nx2,k)=deriv1(bxby(i,1:nx2,k),y)
            enddo
            enddo

 !           %print,'dBxBydx'
            bybz=by*bz
            do j=1,nx2
            do k=1,nx3
             dbybzdz(1:nx1,j,k)=deriv1(bybz(1:nx1,j,k),x)
            enddo
            enddo

!%************* BEGIN INTEGRATION ****************************
            f=dbxbydy+dbxbzdz
            g=dbxbydx+dbybzdz


            do i=1,nx1

            do k=1,nx3
              do j=1,nx2
               sumi=inte1(reshape(f(i,1:nx2,k),(/nx2/)),y(2)-y(1))
               bvarix(i,j,k)=sumi
             enddo
            enddo

            do j=1,nx2
              do k=1,nx3
               sumi=inte1(reshape(g(i,j,1:nx3),(/nx3/)),z(2)-z(1))
               bvariy(i,j,k)=sumi
             enddo
            enddo
            write(*,*) i

            enddo

            bvari=((bvarix+bvariy)/2)-(bz*bz)/2

            do j=1,nx2
            do k=1,nx3
             dpdz(1:nx1,j,k)=deriv1(bvari(1:nx1,j,k),x)
            enddo
            enddo

            bxbybz=(bx*bx+by*by-bz*bz)/2

            do j=1,nx2
            do k=1,nx3
             dbxbybzdz(1:nx1,j,k)=deriv1(bxbybz(1:nx1,j,k),x)
            enddo
            enddo


            bxbz=bx*bz

            do i=1,nx1
            do k=1,nx3
             dbxbzdx(i,1:nx2,k)=deriv1(bxbz(i,1:nx2,k),y)
            enddo
            enddo

            bybz=by*bz

            do i=1,nx1
            do j=1,nx2
             dbybzdy(i,j,1:nx3)=deriv1(bybz(i,j,1:nx3),z)
            enddo
            enddo

            rho1=(dbxbybzdz-dbxbzdx-  dbybzdy+dpdz)/gs

            rho1(1:nx1,1:nx2,1:nx3)=ssimdata%w(1:nx1,1:nx2,1:nx3,4)+rho1+ssimdata%w(1:nx1,1:nx2,1:nx3,13)
            pres=bvari+ssimdata%w(1:nx1,1:nx2,1:nx3,8)*(fgamma-1.0d0)


!%lower boundary

            do i=4,2,-1
              do j=1,nx2
              do k=1,nx3
                     p_2=rho1(i,j,k)*gs
                     pres(i-1,j,k) = (x(2)-x(1))*p_2+pres(i,j,k)
              enddo
              enddo
             enddo


!%upper boundary

            do i=nx1-2,nx1-1
               do j=1,nx2
               do k=1,nx3
                       p_2=rho1(i,j,k)*gs
                       pres(i+1,j,k) = -(x(2)-x(1))*p_2+pres(i,j,k)
               enddo
               enddo
            enddo




!%update the background energy and magnetic fields
           ssimdata%w(1:nx1,1:nx2,1:nx3,4)=0
           ssimdata%w(1:nx1,1:nx2,1:nx3,8)=0

            ssimdata%w(1:nx1,1:nx2,1:nx3,13)=rho1(1:nx1,1:nx2,1:nx3)
            ssimdata%w(1:nx1,1:nx2,1:nx3,12)=pres(1:nx1,1:nx2,1:nx3)/((fgamma-1.0d0))+ &
                0.5d0*(bx(1:nx1,1:nx2,1:nx3)*bx(1:nx1,1:nx2,1:nx3) &
                +bz(1:nx1,1:nx2,1:nx3)*bz(1:nx1,1:nx2,1:nx3)+by(1:nx1,1:nx2,1:nx3)*by(1:nx1,1:nx2,1:nx3))
            ssimdata%w(1:nx1,1:nx2,1:nx3,14)=bz(1:nx1,1:nx2,1:nx3)
            ssimdata%w(1:nx1,1:nx2,1:nx3,15)=bx(1:nx1,1:nx2,1:nx3)
            ssimdata%w(1:nx1,1:nx2,1:nx3,16)=by(1:nx1,1:nx2,1:nx3)


           emin=1.0d99
           rhomin=1.0d99
           lowval=1.0d-12
!          get min energy which is gt  0
           do i=1,nx1
                do j=1,nx2
                    do k=1,nx3
                        if ( (ssimdata%w(i,j,k,12) .gt. 0.0d0) &
                            .and. (ssimdata%w(i,j,k,12) .lt. emin)) then
                            emin=ssimdata%w(i,j,k,12)
                        end if
                    if ( (ssimdata%w(i,j,k,13) .gt. 0.0d0) &
                        .and. (ssimdata%w(i,j,k,13) .lt. rhomin)) then
                            rhomin=ssimdata%w(i,j,k,13)
                        end if
                    enddo
                enddo
            enddo


!          set negative values to the minimum use floor function
           do i=1,nx1
                do j=1,nx2
                    do k=1,nx3
                        if ( ssimdata%w(i,j,k,12) .lt. 0.0d0 ) then
                            ssimdata%w(i,j,k,12)=emin
                        end if
                        if ( ssimdata%w(i,j,k,13) .lt. 0.0d0 ) then
                            ssimdata%w(i,j,k,13)=rhomin
                        end if
                    enddo
                enddo
            enddo



       	endsubroutine hsbalancefield

       	function inte1(f,dx)

            real, intent(in), dimension(:) :: f
            real :: dx, res, inte1
            integer :: i,nel

            nel=size(f)

!            if nsiz(1)>nsiz(2)
!                nel=nsiz(1)
!            else
!                nel=nsiz(2)
!            end

            res=0.0
            if (nel > 1) then
                if (nel == 2) then
                    res=dx*0.5*(f(2)+f(1))
                end if

                if (nel > 2) then
                  do i=2,nel
                      res=res+0.5*(f(i-1)+f(i))*dx
                  enddo
                endif
            endif

            inte1=res

        end function




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
