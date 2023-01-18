module statbalancemod
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


    double precision, parameter :: pi=4.*atan(1.0)
    double precision, parameter :: mu_thermal=0.6d0
    double precision, parameter :: R=8.31e3

! adiabatic gas parameter and the magnetic permeablity
    double precision, parameter :: fgamma=1.66666667e0,mumag=4*pi/1.0e7
!   density kg/m^3
    double precision, parameter :: rho0=2.34d-4

!   solar gravity m/s^2
    double precision, parameter :: gs=-274.0

!   temperature profile parameters (K))
    double precision, parameter :: Tch=8000, Tc=1.8d6

!   position of the transition zone (ytr) and width of transition zone (wtr) in metres
    double precision, parameter :: ytr=2.0d6, wtr=0.02d6


    public :: pi
    public :: mu_thermal, R
    public :: rho0, gs
    public :: Tch, Tc, ytr, wtr

    private


contains

!compute temp at height using tanh function
real*8 function temp( height )
    real*8, intent(in) :: height

    temp=height
end function

!do a trapezium integral here
real*8 function presbalanceint(height)
    real*8 , intent(in) :: height

    presbalanceint=height
end function

real*8 function press( temp, dens )
    real*8 , intent(in) :: temp,  dens

    pres=temp*dens
end function



end module
