! types and general functions used for the model builder
module typesmod
    implicit none

! labels for field items
! rho 0
! mom1 1
! mom2 2
! mom3 3
! e    4
! bx   5
! by   6
! bz   7
! eb   8
! rhob 9
! bxb  10
! byb  11
! bzb  12

! B[Tesla] = B[Gauss] x 10000 * root(mu)
! B[Gauss] is sac units

    include 'param.inc'

    private
    public simgridinfo, simdata, simparams, mconsts

	type simgridinfo
		integer ndimensions(3)
        !%grid_dimensions=64*ones(3,1);
        integer grid_dimensions(3)
        integer grid_left_index(3)
	end type simgridinfo

	type simdata
		real, dimension(128,128,128,16) :: w
	end type simdata

	type simparams
            character(100) :: uniqueidentifier

	        integer boundary_conditions(3)
        	integer dimensionality

        	real domain_left_edge(3)
         	real domain_right_edge(3)
         	integer :: currentiteration
         	real :: currenttime
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

contains





endmodule
