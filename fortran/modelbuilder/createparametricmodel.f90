program createparametricmodel

! program to generate solar atmosphere configuration magnetohydrostatic equilibrium

! Use the following python routine to generate an initial python routine
!https://github.com/mikeg64/solar/blob/master/python/atmos_builder/buildatmos1.ipynb

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


    use modbuildermod, only: writefile, hydropres, hydrodens, dens, hydropres2, bruntvaisalla, temp
    use typesmod
    use generalmod, only: writesac3d
    use fieldbuildermod, only:  genssfield, genvecfield, hsbalancefield

	implicit none

    include 'param.inc'

    real :: height, ldens, ltt, lpres, deltah   !parameters read from file a updated
    real :: h(4096), ld(4096), ltta(4096), lp(4096), lp1(4096),energg(4096)

    real :: iniene ! initial energy at photosphere obtained from VALIIc
    integer :: ii,i,j,k



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

    height=xmin
    deltah=(xmax-xmin)/nx1
    do i=1,nx1
        !calculate height
        h(i)=height
        ! calculate temperature
        ltta(i)=temp(h(i))
        !move to next height
        height=height+deltah
    enddo

    do i=1,nx1
        !sdens(i)=dens(sh(i),spres(i))
        !ii=1+nx1-i
        ii=i
        ld(ii)=hydrodens(h, ii, nx1, deltah)
    enddo





!heights, npoints, deltah, pres, dens %modelbuildermod
! calculate the pressures from the input VALIIIc model data
!call hydropres2(h,nx1,h(2)-h(1), lp,ld)


    do i=1,nx1
        !sdens(i)=dens(sh(i),spres(i))
        !ii=1+nx1-i
        ii=i
        lp(ii)=hydropres(h, ii, nx1, deltah)
    enddo




!energg(i)=presg(i)/(consts.fgamma -1);
!do i=1,nx1
    energg(1:nx1)=lp(1:nx1)/(consts%fgamma-1.0d0)
!end do

! rebuild the model
 do i=1,nx1
    do j=1,nx2
        do k=1,nx3
            sdata%w(i,j,k,13)=ld(i)
            sdata%w(i,j,k,12)=energg(i)
        end do
    end do
 end do

!TODO write coordinates to sdata w(...1) w(...1) w(...1) etc

! generate the field % fieldbuildermod
call genssfield(params, gridinfo, sdata, ifieldtype)

! rebalance the hydrostatic atmosphere pressure
call hsbalancefield(params, gridinfo, sdata)

! write the final input file (ascii or binary) %generalmod
call writesac3d( newfilename, params, gridinfo, sdata, consts )

       	contains





end program createparametricmodel









