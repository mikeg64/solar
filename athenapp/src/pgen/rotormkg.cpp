//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rotor.cpp
//! \brief Sets up 2D rotor test problem.
//!
//! The center of the grid is assumed to have coordinates (x1,x2) = [0,0]; the grid
//! initialization must be consistent with this
//!
//! REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes", JCP, 161,
//!   605 (2000)
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif



namespace {
// made global to share with BC functions
Real grav_acc;
} // namespace



//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Rotor test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gm1 = peos->GetGamma() - 1.0;

  // Read initial conditions from 'athinput'
  Real v0  = pin->GetReal("problem","v0");
  Real p0  = pin->GetReal("problem","p0");
  Real bx0 = pin->GetReal("problem","bx0");
  Real r0  = pin->GetReal("problem","r0");
  Real r1  = pin->GetReal("problem","r1");

  v0 = 0.0;

  grav_acc = phydro->hsrc.GetG2();

  // Initialize the grid.  Note the center is always assumed to have coordinates
  // x1=0, x2=0; the grid range in the input file must be consistent with this
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IDN,k,j,i) = 60.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        /*if(pcoord->x2v(j)< 0.2)
        {
            phydro->u(IDN,k,j,i) = 60.0;
        }
        else
        {
            Real frac = exp(-(SQR((0.3+pcoord->x2v(j))/0.08)));
            phydro->u(IDN,k,j,i) = 1.0+9.0*frac;
        }*/
        // reset density, velocity if cell is inside rotor
        /*Real rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
        if (rad <= r0) {
          phydro->u(IDN,k,j,i) = 10.0;
          phydro->u(IM1,k,j,i) = -100.0*v0*pcoord->x2v(j);
          phydro->u(IM2,k,j,i) = 100.0*v0*pcoord->x1v(i);
        } else {
          // smooth solution between r0 and r1.  For no smoothing, set r1<0 in input
          if (rad <= r1) {
            Real frac = (0.115 - rad)/(0.015);
            phydro->u(IDN,k,j,i) = 1.0 + 9.0*frac;
            phydro->u(IM1,k,j,i) = -frac*100.0*v0*pcoord->x2v(j);
            phydro->u(IM2,k,j,i) =  frac*100.0*v0*pcoord->x1v(i);
          }
        }*/
      }
    }
  }

  // initialize interface B
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie+1; i++) {
        pfield->b.x1f(k,j,i) = bx0;
      }
    }
  }
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x2f(k,j,i) = 0.0;
      }
    }
  }
  for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        pfield->b.x3f(k,j,i) = 0.0;
      }
    }
  }

  // initialize total energy
  //if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*bx0*bx0 +
                                 (SQR(phydro->u(IM1,k,j,i))
                                  + SQR(phydro->u(IM2,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  //}

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  // do nothing
  Real qt,tdep,s_period,AA;
  Real delta_x, delta_y, delta_z, xxmax, yymax, xxmin, yymin;
  Real exp_x,exp_y,exp_z,exp_xyz;
  Real r1,r2,xp, yp,zp;
  Real vvy;
  Real x1,x2,x3;
  Real dt,time;

  Real xcy,xcx;

  int n1,n2;
  int i,j,k;

  n1=2;
  n2=2;

  s_period=50.0; //Driver period
  AA=1.0;       //Driver amplitude
  //AA=1;
  xcy=-0.4;
  xcx=0.0;
  delta_z=0.01;
  delta_x=0.01;
  delta_y=0.01;

  qt=this->pmy_mesh->time;
  dt=this->pmy_mesh->dt;
  tdep=sin(qt*2.0*PI/s_period);

    qt=time;
  // update vy
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {

        x1=pcoord->x1v(i);
	    x2=pcoord->x2v(j);

		r2=(x2-xcy)*(x2-xcy);
        r1=(x1-xcx)*(x1-xcx);

        exp_y=exp(-r2/(delta_y*delta_y));
		//exp_z=exp(-r2/(delta_z*delta_z));
        exp_x=exp(-r1/(delta_x*delta_x));

		//exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
	    exp_xyz=exp_x*exp_y;
	    vvy=AA*exp_xyz*tdep;
        phydro->u(IM2,k,j,i)+= (dt)*vvy*phydro->u(IDN,k,j,i);

      }
    }
  }

//pGrid->U[k][j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[k][j][i].d)/2.0;
   // update total energy
//  if (NON_BAROTROPIC_EOS) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) = dt*SQR(vvy)*phydro->u(IDN,k,j,i)/2.0;
        }
      }
    }
//  }


  return;
}
