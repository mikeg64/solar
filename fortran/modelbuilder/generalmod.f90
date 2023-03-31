
module generalmod
    implicit none

    private
    public simgridinfo, simdata, simparams, mconsts

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

contains

    !compute temp at height using tanh function
    real function temp( height )
        real, intent(in) :: height
        real :: tmptemp

        tmptemp=1+tanh((height-ytr)/wtr)

        temp=Tch+((Tc-Tch)/2.0)*tmptemp
    end function


endmodule
