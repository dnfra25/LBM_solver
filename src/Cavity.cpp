#include "Cavity.hpp"

#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

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
    (
        0.5 + 3.0*viscosity
    )
),

omega_m(0.0),

lattice(nx_,
        ny_,
        omega_p,
        omega_m,
        1.0),

boundary(nx_,
         ny_,
         U_lid)

{

    /*
        TRT parameters

        (tau_p - 0.5)(tau_m - 0.5)=Lambda

        Lambda = 3/16
    */


    double tau_p =
        1.0 / omega_p;



    double tau_m =
        0.5
        +
        (1.0/4.0)
        /
        (tau_p - 0.5);



    omega_m =
        1.0 / tau_m;


    // IMPORTANT:
    // update TRT antisymmetric relaxation
    // after Lattice construction

    lattice.setOmegaM(omega_m);



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

        1) macroscopic variables
        2) equilibrium
        3) TRT collision
        4) streaming
        5) boundary reconstruction
        6) update macroscopic variables
    */


    lattice.computeMacroscopic();


    lattice.computeEquilibrium();


    lattice.collision();


    lattice.streaming();


    boundary.apply(lattice);


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

    double num = 0.0;

    double den = 0.0;



    #pragma omp parallel for collapse(2) reduction(+:num,den)
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
                ux - oldUx[id];


            double duy =
                uy - oldUy[id];



            num +=
                dux*dux
              +
                duy*duy;



            den +=
                ux*ux
              +
                uy*uy;



            oldUx[id]=ux;

            oldUy[id]=uy;

        }
    }



    return
        std::sqrt(
            num/(den+1e-30)
        );

}



//==================================================
// Convergence test
//==================================================

bool Cavity::converged(double tolerance)
{

    return
        velocityDifference()
        <
        tolerance;

}



//==================================================
// Write horizontal centerline Ux
//==================================================

void Cavity::writeCenterlineUx(
        const std::string& filename)
{

    std::ofstream file(filename);


    file
        << "# x ux\n";


    int y =
        (ny-1)/2;



    file
        << std::scientific
        << std::setprecision(12);



    for(int x=0; x<nx; x++)
    {

        double X =
            static_cast<double>(x)
            /
            static_cast<double>(nx-1);



        double ux =
            lattice.getUx(x,y);



        file
            << X
            << " "
            << ux
            << "\n";

    }


    file.close();

}



//==================================================
// Write vertical centerline Uy
//==================================================

void Cavity::writeCenterlineUy(
        const std::string& filename)
{

    std::ofstream file(filename);


    file
        << "# y uy\n";


    int x =
        (nx-1)/2;



    file
        << std::scientific
        << std::setprecision(12);



    for(int y=0; y<ny; y++)
    {

        double Y =
            static_cast<double>(y)
            /
            static_cast<double>(ny-1);



        double uy =
            lattice.getUy(x,y);



        file
            << Y
            << " "
            << uy
            << "\n";

    }


    file.close();

}
