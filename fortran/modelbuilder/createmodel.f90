program createmodel

! program to generate solar atmosphere configuration magnetohydrostatic equilibrium

! refer to
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/generatefield_verttube.m

!
!writing routine
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/hdf5_and_gdf/writesac3D.m


!todo  % complete selfsimilarity field generator genssfield%fieldbuildermod
!todo  % complete sac output writer writesac3d%generalmod
!todo  % check hydropres2 integrating pressure correctly hydropres2%modbuildermod
!todo  % construction of initialisation routine check correct calling of typesmod
!todo  % correct dynamic creation of routines - with arrays set to correct size need param.inc types.mod and initialisation routine
!todo  % replace params.inc with correct use of namelists e.g. see
!        https://github.com/mikeg2105/comp-sci/blob/master/gfortran/advanced/f95features/readnamelist.f90


    use modbuildermod, only: writefile, hydropres, hydrodens, dens, hydropres2, bruntvaisalla
    use typesmod
    use generalmod, only: writesac3d
    use fieldbuildermod, only:  genssfield, genvecfield

	implicit none

    include 'param.inc'

    real :: height, ldens, ltt, lpres   !parameters read from file a updated
    real :: h(4096), ld(4096), ltta(4096), lp(4096), lp1(4096),energg(4096)

    real :: iniene ! initial energy at photosphere obtained from VALIIc
    integer :: i,j,k



	character (len=30), parameter :: newfilename='3D_128_4Mm_asc.ini'

	type (mconsts)::consts
    type (simgridinfo) :: gridinfo
    type (simparams) :: params
    type (simdata) :: sdata

	consts%mu_therm=mu_thermal
	consts%R=R
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

    h=0.0d0
    ld=0.0d0
    ltta=0.0d0
    lp=0.0d0
    lp1=0.0d0
    energg=0.0d0

  !read the input file with the scale distances
  ! create an atmospheric profile
  ! Open file for unformmated output using the inquired record length
  !use interpolation https://docs.scipy.org/doc/scipy/tutorial/interpolate/1D.html
  ! use original VALIIIc data set to generate the correct atmosphere profile
  !e.g. see https://github.com/mikeg64/solar/blob/master/matlab/solatmoswaves/atmos.csv
  !
  ! number of rows must be same as nx3
  ! columns must be as follows
  ! height temperature density pressure
  open(unit = 1, file = 'atmos1.txt',status='old', action='read', iostat=rstat)
  do i = 1,nx3
      read(1,'(6X,d8.3,7x,d9.7,7x,d9.7,7x,d9.7)',iostat=rstat) height, ltt, ldens, lpres
      h(i)=height
      ld(i)=ldens
      ltta(i)=ltt
      lp(i)=lpres
      print *, height, ldens, ltt, lpres
  end do
  !read(1,3(1X, F10.0)) height, logdens, logtt
  close(unit = 1)

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


 do i=1,nx1
    do j=1,nx2
        do k=1,nx3
            sdata%w(i,j,k,10)=ld(k)
            sdata%w(i,j,k,9)=iniene
        end do
    end do
 end do

! use energy to get pthermal


!do i=1,nx1
    lp(1:nx3)=(consts%fgamma-1.0d0)*iniene
!end do

lp1=lp;

!heights, npoints, deltah, pres, dens %modelbuildermod
call hydropres2(h,nx3,h(2)-h(1), lp,ld)


!energg(i)=presg(i)/(consts.fgamma -1);
!do i=1,nx1
    energg(1:nx3)=lp(1:nx3)/(consts%fgamma-1.0d0)
!end do

! rebuild the model
 do i=1,nx1
    do j=1,nx2
        do k=1,nx3
            sdata%w(i,j,k,10)=ld(k)
            sdata%w(i,j,k,9)=energg(k)
        end do
    end do
 end do


! generate the field % fieldbuildermod
    call genssfield(params, gridinfo, sdata, ifieldtype)


! write the final input file (ascii or binary) %generalmod
    call writesac3d( newfilename, params, gridinfo, sdata, consts )

       	contains





end program createmodel









