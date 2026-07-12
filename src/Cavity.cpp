#include "Cavity.hpp"

#include <cmath>
#include <algorithm>



static double computeOmegaM(double viscosity)
{

    double tau_p =
        0.5 + 3.0*viscosity;


    double tau_m =
        0.25/(tau_p-0.5)+0.5;


    return 1.0/tau_m;

}



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

viscosity(
    U_lid*(nx_-1)/Re
),

omega_p(
    1.0/
    (
        0.5+3.0*viscosity
    )
),

omega_m(
    computeOmegaM(viscosity)
),

lattice(
    nx_,
    ny_,
    omega_p,
    omega_m,
    1.0
),

boundary(
    nx_,
    ny_,
    U_lid
)

{

    oldUx.resize(nx*ny,0.0);
    oldUy.resize(nx*ny,0.0);

}



//==================================================

void Cavity::initialize()
{

    lattice.initialize();


    std::fill(
        oldUx.begin(),
        oldUx.end(),
        0.0
    );


    std::fill(
        oldUy.begin(),
        oldUy.end(),
        0.0
    );

}



//==================================================

void Cavity::step()
{

    lattice.computeMacroscopic();

    lattice.computeEquilibrium();

    lattice.collision();

    lattice.streaming();

    boundary.apply(lattice);

    lattice.computeMacroscopic();

}



//==================================================

void Cavity::applyBoundary()
{
    boundary.apply(lattice);
}



//==================================================

Lattice& Cavity::getLattice()
{
    return lattice;
}



//==================================================

double Cavity::velocityDifference()
{

    double num=0.0;

    double den=0.0;



#pragma omp parallel for collapse(2) reduction(+:num,den)
    for(int x=0;x<nx;x++)
    {
        for(int y=0;y<ny;y++)
        {

            int id=x*ny+y;


            double ux =
                lattice.getUx(x,y);

            double uy =
                lattice.getUy(x,y);



            double dux =
                ux-oldUx[id];

            double duy =
                uy-oldUy[id];



            num +=
                dux*dux+
                duy*duy;


            den +=
                ux*ux+
                uy*uy;



            oldUx[id]=ux;

            oldUy[id]=uy;

        }
    }


    return sqrt(num/(den+1e-30));

}



bool Cavity::converged(double tolerance)
{

    return velocityDifference()<tolerance;

}
