program hydrostatsolatmos


    ! Routine uses tanh based solar temperature profile
    ! and hrdrostatic pressure balance

    !        height(m)        temp(K)     dens(kg)          pres()
    !        5955555.6       1599238.9   5.5055286e-12      0.12194448
    use statbalancemod, only: temp, writefile, hydropres, hydrodens, dens, hydropres2, bruntvaisalla

    implicit none


    real, parameter :: deltah = 22700.0d0
    !real, parameter :: deltah = 45454.545455

    ! value to include photosphere starth is 200.0
    !real, parameter :: starth = 2138469.0
    real, parameter :: starth = 200.0
    !we save the laST 132 points because the first three are -ve energy density
    integer, parameter :: npoints=293
    !real, parameter :: pi=4.0*atan(1.0)
    real,allocatable :: sh(:),sdens(:),spres(:),stemp(:),sbruntvas(:)
    real :: hcurrent
    integer :: i, ii, presmethod=2


    allocate(sh(1:npoints))
    allocate(sdens(1:npoints))
    allocate(spres(1:npoints))
    allocate(stemp(1:npoints))
    allocate(sbruntvas(1:npoints))

    ! initialise array values
    sh(:)=0
    sdens(:)=0
    spres(:)=0
    stemp(:)=0
    sbruntvas(:)=0
    !perform the computation loop
    !hcurrent=deltah
    hcurrent=starth

    do i=1,npoints
        !calculate height
        sh(i)=hcurrent
        hcurrent=sh(i)
        hcurrent=hcurrent+deltah

        ! calculate temperature
        stemp(i)=temp(sh(i))

    enddo

    do i=1,npoints

        !sdens(i)=dens(sh(i),spres(i))
        ii=1+npoints-i
        sdens(ii)=hydrodens(sh, ii, npoints, deltah)

    enddo


    !pres2 calculation
    !hydropres2(heights, hindex, npoints, deltah, pres, dens)
    call hydropres2(sh, npoints, deltah, spres, sdens)


    !!hcurrent=starth
if (presmethod.eq.1) then
    do i=1,npoints
        !calculate height
        !sh(i)=hcurrent+deltah
        ii=1+npoints-i
        hcurrent=sh(1+npoints-i)


        !use hydrostatic balance to calculate  density
        spres(ii)=hydropres(sh, ii, npoints, deltah)
        print *, ii,' height',sh(ii),' dens ',spres(ii)

        !compute pressure

        !write to ascii output file

!        print *, 'result'
    enddo
endif



    call bruntvaisalla(sh, npoints, deltah, spres, sdens, sbruntvas)

    call writefile(sh,sdens,spres,stemp,sbruntvas,npoints)

    print *, 'complete'
end program  hydrostatsolatmos
