program createmodel

! program to generate solar atmosphere configuration magnetohydrostatic equilibrium


	implicit none

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

    real :: height, logdens, logtt   !parameters read from file a updated
    real, dimension(305) :: h, ld, ltt

    real :: inene ! initial energy at photosphere obtained from VALIIc
    integer :: i
	integer :: ngx1=2
	integer :: ngx2=2
	integer :: ngx3=2

	integer :: it=0
	integer :: rstat

    integer :: ndim=3
    integer :: neqpar=7
    integer :: nw=13
    integer :: nx1=128
    integer :: nx2=128
    integer :: nx3=305   ! made same as size of input file


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


	character (len=30), parameter :: newfilename='3D_128_4Mm_asc.ini'

	type simgridinfo
		integer ndimensions
        !%grid_dimensions=64*ones(3,1);
        integer grid_dimensions(3)
        integer grid_left_index(3)
	end type simgridinfo

	type simdata
		real w(256,256,256,16)
	end type simdata

	type simparams
	        integer boundary_conditions(3)
        	integer dimensionality

        	real domain_left_edge(3)
         	real domain_right_edge(3)
        	real eta
        	real adiab
        	real fgamma
        	real gravity0
        	real gravity1
        	real gravity2
        	real nu
        	integer num_ghost_zones

        	integer nw
        	integer neqpar
	end type simparams

	type mconsts
        real :: R
		real :: mu
		real :: fgamma
		real :: ggg
		real :: mu_therm
	end type mconsts

	type (mconsts)::consts
    type (simgridinfo) :: gridinfo
    type (simparams) :: params
    type (simdata) :: sdata

	consts%mu_therm=0.6e0
	consts%R=8.31e3
	consts%fgamma=rgamma
	consts%ggg=g1
	consts%mu=4*pi/1.0e7






    !set up simparams, simgridinfo, simdata

    gridinfo%grid_dimensions(1)= nx1
    gridinfo%grid_dimensions(2)= nx2
    gridinfo%grid_dimensions(3)= nx3
    gridinfo%ndimensions = ndim

    params%adiab=adiab
    params%dimensionality=ndim
    params%domain_left_edge(1)=xmin
    params%domain_left_edge(2)=ymin
    params%domain_left_edge(3)=zmin
    params%domain_right_edge(1)=xmax
    params%domain_right_edge(2)=ymax
    params%domain_right_edge(3)=zmax
    params%eta=0
    params%fgamma=rgamma
    params%gravity0=g1
    params%gravity0=g2
    params%gravity0=g3
    params%neqpar=neqpar

    sdata%w=0

    !read the input file with the scale distances

  ! Open file for unformmated output using the inquired record length
  open(unit = 1, file = 'stratification-fromphoto.dat',status='old', action='read', iostat=rstat)
  do i = 1,nx3
      read(1,'(6X,d8.3,7x,d9.7,7x,d9.7)',iostat=rstat) height, logdens, logtt
      h(i)=height
      ld(i)=logdens
      ltt(i)=logtt
      print *, height, logdens, logtt
  end do
  !read(1,3(1X, F10.0)) height, logdens, logtt
  close(unit = 1)

  print *, height, logdens, logtt
    print *, newfilename





!compute correct pressure for gravitationally stratified atmosphere

!compute initial energy (at photosphere or temperature minimum)
!mu_thermal=0.6d0;
!R=8.31e3;

! pressure=temp*R*density/((mu_thermal))
!parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
!iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)

 !iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)

 !iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)

 ! 1.6Mm

iniene=6840.0e0*consts%R*(2.3409724e-09)/consts%mu_therm/(consts%fgamma-1.0)

 !iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)

 !iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)





 ! Open file for unformmated output using the inquired record length
  open(unit = 1, file = 'stratification-fromphoto-hstatic.dat',status='new', action='write', iostat=rstat)
  do i = 1,nx3
      height=h(i)
      logdens=ld(i)
      logtt=ltt(i)
      print *, height, logdens, logtt
      !write(1,'(6X,d8.3,7x,d9.7,7x,d9.7)',iostat=rstat) height, logdens, logtt
      write(1,'(f20.14,7X, f20.14, 7X, f20.14)',iostat=rstat) height, logdens, logtt
!40   format (3d20.14)
  end do
  !read(1,3(1X, F10.0)) height, logdens, logtt
  close(unit = 1)



       	contains

       	![simparams, simgridinfo, simdata]=
       	subroutine genfield(ssimparams, ssimgridinfo, ssimdata, fieldtype)



            integer, intent(in) :: fieldtype
       		type(simparams), intent(inout) :: ssimparams
       		type(simgridinfo), intent(inout) :: ssimgridinfo
       		type(simdata), intent(inout)  :: ssimdata





       	end subroutine genfield



end program createmodel









