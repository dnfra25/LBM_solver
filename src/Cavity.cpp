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


/*
    For a square cavity:

        Re = U L / nu

    with:

        L = 1

    therefore:

        nu = U/Re
*/

viscosity(U_lid / Re),



/*
    LBM relation:

        nu = cs^2 (tau - 1/2)

    with:

        cs^2 = 1/3

*/

omega_p(
    1.0 /
    (0.5 + 3.0*viscosity)
),



/*
    TRT parameter

        (tau_p-0.5)(tau_m-0.5)=lambda

*/

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

    // Apply initial boundary conditions
    boundary.apply(lattice);


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
    lattice.computeMacroscopic();
    lattice.computeEquilibrium();
    lattice.collision();
    lattice.streaming();
    boundary.apply(lattice);
    // aggiorna campi macroscopici dopo BC
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



    #pragma omp parallel for reduction(max:error)
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



            double dux =
                std::abs(
                    ux-oldUx[id]
                );


            double duy =
                std::abs(
                    uy-oldUy[id]
                );



            error =
                std::max(error,
                         std::max(dux,duy));



            oldUx[id]=ux;

            oldUy[id]=uy;

        }
    }


    return error;

}






//==================================================
// Convergence test
//==================================================

bool Cavity::converged(double tolerance)
{
    double error = velocityDifference();

    if(!std::isfinite(error))
    {
        return false;
    }

    return error < tolerance;
}
