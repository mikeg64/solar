module modbuildermod

    use typesmod

    implicit none
!mu=0.6d0;
!R=8.31e3;
! parrVALMc=rhoarrVALMc*TarrVALMc*R/mu

!rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
! p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
! p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
!          (w[*,*,*,0]+w[*,*,*,9])/2.0
! p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
!           +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
! p[*,*,*]=(gamma-1.d0)*p[*,*,*]


!compute correct pressure for gravitationally stratified atmosphere

!compute initial energy (at photosphere or temperature minimum)
!mu_thermal=0.6d0;
!R=8.31e3;

! temp*R*density/((mu_thermal))
!parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
!iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)

! !iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)
!
! !iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)
!
! ! 1.6Mm
!
! !iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=6840.d0*R*(2.3409724e-09)/mu/(consts.fgamma-1.0);

!%
!% !iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)
!%
!% !iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)



    real, parameter :: pi=4.*atan(1.0)
    real, parameter :: mu_mass=2.2d0          !0.6d0
    real, parameter :: R=8.31d3
    real, parameter :: mu_thermal = 0.6

! adiabatic gas parameter and the magnetic permeablity
    real, parameter :: fgamma=1.66666667e0,mumag=4*pi/1.0e7
!   density kg/m^3
!    real, parameter :: rho0=2.34d-4, p0=9228.6447 or rho0=5.48d-12
!   pressure and density taken fromVALIIIc data
!    real, parameter :: rho0=2.34d-4,p0=2.7865d-4
!    real, parameter :: rho0=4.68d-4,p0=2.7865d-4
    real, parameter :: rho0=1.17d-4,p0=2.7865d-4
!    real, parameter :: rho0=0.885d-4,p0=2.7865d-4
!    real, parameter :: rho0=0.00017d-4,p0=2.2865d-4
!   solar gravity m/s^2
    real, parameter :: gs=274.0

!   temperature profile parameters (K))
    real, parameter :: Tch=8000, Tc=1.8d6

!   position of the transition zone (ytr) and width of transition zone (wtr) in metres wtr=0.02
    real, parameter :: ytr=2.0d6, wtr=0.1d6


    public :: pi
    public :: mu_mass, R,fgamma,mumag
    public :: rho0, p0, gs
    public :: Tch, Tc, ytr, wtr
    public :: writefile, temp, hydropres, dens, hydrodens, hydropres2, bruntvaisalla


private


contains

!compute temp at height using tanh function
real function temp( height )
    real, intent(in) :: height
    real :: tmptemp

    tmptemp=1+tanh((height-ytr)/wtr)
    temp=Tch+((Tc-Tch)/2.0)*tmptemp
end function



!compute pres
real function hydropres(heights, hindex, npoints, deltah)
    real, intent(in) :: deltah, heights(4096)
    integer, intent(in) :: hindex, npoints
    real :: psum,Hscale,p0c
    integer :: i

!    tmp0=temp(heights(1))
    psum=0.0

    p0c=rho0*R*temp(heights(1))/mu_mass



    if (hindex.eq.1) then
        Hscale=R*temp(heights(hindex))/(mu_mass*gs)
!        tmptemp=p0*tmp0/temp(heights(hindex))
        psum=psum+p0c*exp(-deltah/Hscale)
!         psum=deltah/Hscale
    elseif (hindex.le.npoints) then
 !       tmptemp=p0*tmp0/temp(heights(hindex))
        do i=1,hindex,1
            Hscale=R*temp(heights(i))/(mu_mass*gs)
            psum=psum+deltah/Hscale
        end do
        psum=p0c*exp(-psum)
    endif
    print*,'psum ',hindex,' ',psum
    hydropres=psum


end function

subroutine hydropres2(heights, npoints, deltah, pres, dens)



    real, intent(inout) :: pres(4096), dens(4096)
    real, intent(in) :: deltah, heights(4096)
    integer, intent(in) :: npoints
    integer :: i
    !6840
    real :: iniene = 5840.d0*R*(0.00252e-1)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 6840.d0*R*(2.3409724e-09)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 12840.d0*R*(2.3409724e-09)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 231287.68d0*R*(1.09e-11)/mu_thermal/(fgamma-1.0)
    real  :: p_1, p_2, comi



!iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)
! 1.6Mm
!!!!!!iniene=6840.d0*R*(2.3409724e-09)/mu_thermal/(consts.fgamma-1.0);
!iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)



    do i=1,npoints
!       presg(i)=(consts.fgamma-1)*iniene;
        pres(i)=(fgamma-1)*iniene
    enddo



!presg1=presg;

    do i=npoints-1,1,-1
        comi=-abs(heights(i+1)-heights(i))
        pres(i)=pres(i+1)-dens(i)*comi*gs
    enddo


    do i=3,npoints-2
        comi=-abs(heights(i+1)-heights(i))
!     %densg(i)=densg(i)-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
        dens(i)=(1.0/gs)*(  (1.0/(12*(heights(i+1)-heights(i)))) *(pres(i+2)-8*pres(i+1)+8*pres(i-1)-pres(i-2))     )
    enddo

    do i=5,3,-1
        p_1=pres(i+2)+8*pres(i+1)-8*pres(i-1)
        p_2= -dens(i)*gs
        pres(i-2)= -p_1-12.0*(heights(i)-heights(i-1))*p_2
    enddo




endsubroutine



! this routine is currently the same as hydropres2
! but it will be adapted to compute the parametric pressure and density profile
! using the tanh function parametrisation
subroutine hydropres1(heights, npoints, deltah, pres, dens)



    real, intent(inout) :: pres(4096), dens(4096)
    real, intent(in) :: deltah, heights(4096)
    integer, intent(in) :: npoints
    integer :: i
    !6840
    real :: iniene = 5840.d0*R*(0.00252e-1)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 6840.d0*R*(2.3409724e-09)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 12840.d0*R*(2.3409724e-09)/mu_thermal/(fgamma-1.0)
    !real :: iniene = 231287.68d0*R*(1.09e-11)/mu_thermal/(fgamma-1.0)
    real  :: p_1, p_2, comi



!iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)
! 1.6Mm
!!!!!!iniene=6840.d0*R*(2.3409724e-09)/mu_thermal/(consts.fgamma-1.0);
!iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)



    do i=1,npoints
!       presg(i)=(consts.fgamma-1)*iniene;
        pres(i)=(fgamma-1)*iniene
    enddo



!presg1=presg;

    do i=npoints-1,1,-1
        comi=-abs(heights(i+1)-heights(i))
        pres(i)=pres(i+1)-dens(i)*comi*gs
    enddo


    do i=3,npoints-2
        comi=-abs(heights(i+1)-heights(i))
!     %densg(i)=densg(i)-(1.0/consts.ggg)*(  (1.0/(12*(rheight(i+1)-rheight(i)))) *(presg(i+2)-8*presg(i+1)+8*presg(i-1)-presg(i-2))     );
        dens(i)=(1.0/gs)*(  (1.0/(12*(heights(i+1)-heights(i)))) *(pres(i+2)-8*pres(i+1)+8*pres(i-1)-pres(i-2))     )
    enddo

    do i=5,3,-1
        p_1=pres(i+2)+8*pres(i+1)-8*pres(i-1)
        p_2= -dens(i)*gs
        pres(i-2)= -p_1-12.0*(heights(i)-heights(i-1))*p_2
    enddo




endsubroutine








!compute bruntvaisalla frequency
subroutine bruntvaisalla(heights, npoints, deltah, pres, dens, bruntvas)

    real, intent(inout) :: pres(4096), dens(4096), bruntvas(4096)
    real, intent(in) :: deltah, heights(4096)
    real :: presgrad(4096), densgrad(4096)
    integer, intent(in) :: npoints
    integer :: i

    real  :: p_1, p_2, comi



    !compute density gradient
    do i=3,npoints-2
        comi=-abs(heights(i+1)-heights(i))
        p_1= -dens(i+2)+8*dens(i+1)-8*dens(i-1)+dens(i-2)
        p_2= (1.0d0/dens(i))
        densgrad(i)= p_2*(p_1/comi)
    enddo

    !compute pressure gradient
    do i=3,npoints-2
        comi=-abs(heights(i+1)-heights(i))
        p_1= -pres(i+2)+8*pres(i+1)-8*pres(i-1)+pres(i-2)
        p_2= (1.0d0/(fgamma*pres(i)))
        presgrad(i)= p_2*(p_1/comi)
    enddo

    do i=3,npoints-2
         bruntvas(i)= -gs*(densgrad(i)-presgrad(i))
    enddo


endsubroutine




!using hydrostatic mass balance compute the integrated density

!compute pres
real function hydrodens(heights, hindex, npoints, deltah)
    real, intent(in) :: deltah, heights(4096)
    integer, intent(in) :: hindex, npoints
    real :: rsum,Hscale,T0
    integer :: i

!    tmp0=temp(heights(1))

    rsum=0.0
    T0=temp(heights(1))
    if (hindex.eq.1) then
        Hscale=R*temp(heights(hindex))/(mu_mass*gs)
!        tmptemp=p0*tmp0/temp(heights(hindex))
        rsum=rsum+(deltah/Hscale)
!         psum=deltah/Hscale
    elseif (hindex.le.npoints) then
 !       tmptemp=p0*tmp0/temp(heights(hindex))
        do i=1,hindex,1
            Hscale=R*temp(heights(i))/(mu_mass*gs)
            rsum=rsum+deltah/Hscale
        end do

    endif
    rsum=rho0*(T0/temp(heights(hindex)))*exp(-rsum)
    print*,'rhosum ',R,mu_mass,hindex,' ',rsum
    hydrodens=rsum

end function


! given height and pressure calculate the density
! at that height use the ttemp function to compute the temperature
!returns a single value

!use the equation of state

real function dens( height, spres )
    real, intent(in) :: height, spres
    real :: tmpresult

! parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
    tmpresult= mu_mass*spres/(R*temp(height))
    dens=tmpresult

end function

subroutine writefile(height, dens, press, temp,bruntvas, nitems)
    implicit none
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

end subroutine writefile


end module
