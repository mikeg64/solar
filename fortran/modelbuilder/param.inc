
    ! parameters for the model creator


	!
	!  pi and its derivatives.
	!
  	real, parameter :: pi=3.14159265358979323846264338327950d0
  	real, parameter :: pi_1=1./pi,pi4_1=pi**(-4),pi5_1=pi**(-5)
  	real, parameter :: sqrtpi=1.77245385090551602729816748334115d0
  	real, parameter :: sqrt2=1.41421356237309504880168872420970d0
  	real, parameter :: sqrt2pi=sqrt2*sqrtpi
  	real, parameter :: four_pi_over_three=4.0/3.0*pi
  	real, parameter :: onethird=1./3., twothird=2./3., fourthird=4./3., onesixth=1./6.
  	real, parameter :: one_over_sqrt3=0.577350269189625764509148780501958d0
  	real, parameter :: twopi = 6.2831853071795864769252867665590d0
  	real, parameter :: dtor = pi/180.d0
	!
  	real, parameter :: lntwo=0.69314718055995d0

    real, parameter :: mu_mass=2.2d0          !0.6d0
    real, parameter :: R=8.31e3
    real, parameter :: mu_thermal = 0.6

    ! adiabatic gas parameter and the magnetic permeablity
    real, parameter :: fgamma=1.66666667e0,mumag=4*pi/1.0e7

!   density kg/m^3
!    real, parameter :: rho0=2.34d-4, p0=9228.6447
!   pressure and density taken fromVALIIIc data
    real, parameter :: rho0=5.48d-12,p0=2.7865d-4
!   solar gravity m/s^2
    real, parameter :: gs=-274.0

!   ifieldtype1 is self similar field
!   http://solarwavetheory.blogspot.com/2013/11/solar-atmospheric-mhd-flux-tube.html
    integer, parameter :: ifieldtype=1

	integer :: ngx1=2
	integer :: ngx2=2
	integer :: ngx3=2

	integer :: it=0
	integer :: rstat

    integer :: ndim=3
    integer :: neqpar=7
    integer, parameter :: nw=13
    integer, parameter :: nx1=256
    integer, parameter :: nx2=64
    integer, parameter :: nx3=64   ! made same as size of input file
!   character (len=30), parameter :: snewfilename='3D_128_4Mm_asc.ini'
   character (len=30), parameter :: snewfilename='3D_256_64_64_4Mm_asc.ini'

    character(100) :: suniqueidentifier = 'parametricmodel'
    real :: rgamma=1.666667
    real :: adiab=1.0
    real :: eta=0
    real :: g1=-274
    real :: g2=0
    real :: g3=0

    real :: xmin=133333.33
    real :: ymin=1953.1
    real :: zmin=1953.1
    real :: xmax=5955555.6e0
    real :: ymax=4.0e6
    real :: zmax=4.0e6
