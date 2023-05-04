program testarray

    implicit none
            integer, parameter :: nx1=6, nx2=6, nx3=6
            real, dimension(nx3) :: x,z, b0z, dbz, bx, by, bz
       		real, dimension(nx1,nx2,nx3) :: dbzdz, dbydz, dbxdz, dbzdx, dbydx, dbxdx
     		real, dimension(nx1,nx2,nx3) :: dbzdy, dbydy, dbxdy
    		real, dimension(nx1,nx2,nx3) :: bvarix, bvariy, bvari, dpdx, rho1
    		real, dimension(nx1,nx2,nx3) :: bvaridz, bvaridz1, dbvardz, br
    		real, dimension(nx1,nx2,nx3) :: dbxbydy, bxbz, dbxbzdz
    		real, dimension(nx1,nx2,nx3) :: bxby,bybz, dbybzdy, divb, f,g

    		integer, dimension(3) :: shap
            real, dimension(nx2) :: res
            real :: integral

            integer :: i,j,k

            do i=1,nx1
                x(i)=i*0.25d0
            end do

            do i=1,nx1
                do j=1,nx2
                    do k=1,nx3
                        bxby(i,j,k)=cos(3.141d0*(i+j+k)/12.0d0)
                    end do
                end do
             end do

            do k=1,nx3
                do j=1,nx2
                    dbzdz(1:nx1,j,k)=deriv1(bxby(1:nx1,j,k),x)
                enddo
            enddo

            i=1
            k=1
            do j=1, nx2
                write(*,*) bxby(i,j,k), dbzdz(i,j,k)
            end do

            shap=shape(dbzdz)
            write(*,*) shap

            res = reshape(dbzdz(1,1:nx2,1),(/nx2/))
            write(*,*) res

            integral=inte1(res, x(2)-x(1))
            write(*,*) integral

    contains

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

end program
