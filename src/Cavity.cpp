#include "Cavity.hpp"

#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif



//==================================================
// Constructor
//==================================================

Cavity::Cavity(int nx_,
               int ny_,
               double Reynolds_,
               double lidVelocity_)
:
nx(nx_),
ny(ny_),
Re(Reynolds_),
U_lid(lidVelocity_),

viscosity(U_lid / Re),

omega_p(
    1.0 /
    (0.5 + 3.0*viscosity)
),

omega_m(
    1.0 /
    (
        0.5
        +
        0.25 /
        (
            1.0/omega_p - 0.5
        )
    )
),

lattice(nx,
        ny,
        omega_p,
        omega_m,
        1.0),

boundary(nx,
         ny,
         U_lid)

{

    oldUx.resize(nx*ny,0.0);

    oldUy.resize(nx*ny,0.0);

}



//==================================================
// Initialization
//==================================================

void Cavity::initialize()
{

    lattice.initialize();


    std::fill(oldUx.begin(),
              oldUx.end(),
              0.0);


    std::fill(oldUy.begin(),
              oldUy.end(),
              0.0);

}



//==================================================
// One simulation step
//==================================================

void Cavity::step()
{

    /*
        LBM cycle:

        1) macroscopic fields
        2) equilibrium
        3) collision
        4) streaming
        5) boundary conditions
        6) update macroscopic fields

    */


    lattice.computeMacroscopic();


    lattice.computeEquilibrium();


    lattice.collision();


    lattice.streaming();


    boundary.apply(lattice);



    /*
        Important:

        Boundary conditions modify
        distribution functions.

        Therefore rho,u must be
        updated afterwards.
    */

    lattice.computeMacroscopic();

}



//==================================================
// Apply boundary only
//==================================================

void Cavity::applyBoundary()
{

    boundary.apply(lattice);

}



//==================================================
// Access lattice
//==================================================

Lattice& Cavity::getLattice()
{

    return lattice;

}



//==================================================
// Velocity difference
//==================================================

double Cavity::velocityDifference()
{

    double error = 0.0;


    bool invalid = false;



    #pragma omp parallel for collapse(2) reduction(max:error)
    for(int x=0;x<nx;x++)
    {
        for(int y=0;y<ny;y++)
        {

            int id =
                x*ny+y;



            double ux =
                lattice.getUx(x,y);


            double uy =
                lattice.getUy(x,y);



            /*
                Detect numerical divergence
            */

            if(!std::isfinite(ux) ||
               !std::isfinite(uy))
            {
                invalid = true;
                continue;
            }



            double dux =
                std::abs(
                    ux-oldUx[id]
                );


            double duy =
                std::abs(
                    uy-oldUy[id]
                );



            error =
                std::max(
                    error,
                    std::max(dux,duy)
                );



            oldUx[id]=ux;

            oldUy[id]=uy;

        }
    }



    if(invalid)
        return 1e30;



    return error;

}



//==================================================
// Convergence test
//==================================================

bool Cavity::converged(double tolerance)
{

    return velocityDifference()
           <
           tolerance;

}
