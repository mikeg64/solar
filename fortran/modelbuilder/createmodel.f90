program createmodel

! program to generate solar atmosphere configuration magnetohydrostatic equilibrium

! refer to
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/generatefield_verttube.m

!
!writing routine
!https://github.com/mikeg64/smaug_realpmode/blob/master/matlab/generateinitialconfiguration/hdf5_and_gdf/writesac3D.m

    use modbuildermod, only: writefile, hydropres, hydrodens, dens, hydropres2, bruntvaisalla
    use generalmod

	implicit none

    include 'param.inc'

    real :: height, logdens, logtt   !parameters read from file a updated
    real, dimension(305) :: h, ld, ltt

    real :: iniene ! initial energy at photosphere obtained from VALIIc
    integer :: i



	character (len=30), parameter :: newfilename='3D_128_4Mm_asc.ini'

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









