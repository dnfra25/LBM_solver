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

viscosity(U_lid * (nx_ - 1) / Re),

omega_p(
    1.0 /
    (0.5 + 3.0*viscosity)
),

omega_m(
    1.0 /
    (
        0.25 /
        (
            (0.5 + 3.0*viscosity)
            -
            0.5
        )
        +
        0.5
    )
),

lattice(nx_,
        ny_,
        omega_p,
        omega_m,
        1.0),

boundary(nx_,
         ny_,
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

    lattice.computeMacroscopic();

    lattice.computeEquilibrium();

    lattice.collision();

    lattice.streaming();

    boundary.apply(lattice);

   // lattice.computeMacroscopic();

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

    double num = 0.0;
    double den = 0.0;


    #pragma omp parallel for collapse(2) reduction(+:num,den)
    for(int x=0;x<nx;x++)
    {
        for(int y=0;y<ny;y++)
        {

            int id = x*ny+y;


            double ux =
                lattice.getUx(x,y);

            double uy =
                lattice.getUy(x,y);



            double dux =
                ux - oldUx[id];

            double duy =
                uy - oldUy[id];



            num += dux*dux + duy*duy;


            den += ux*ux + uy*uy;



            oldUx[id]=ux;
            oldUy[id]=uy;

        }
    }


    return std::sqrt(num/(den+1e-30));

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
