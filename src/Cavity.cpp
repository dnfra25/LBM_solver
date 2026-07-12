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


viscosity(
    U_lid * (nx_-1) / Re
),


omega_p
(
    1.0 /
    (
        0.5 + 3.0*viscosity
    )
),


omega_m(0.0),


lattice
(
    nx_,
    ny_,
    omega_p,
    omega_m,
    1.0
),


boundary
(
    nx_,
    ny_,
    U_lid
)

{


    /*
        TRT magic parameter

        (tau_p-0.5)(tau_m-0.5)=3/16
    */


    double tau_p =
        1.0/omega_p;



    double tau_m =
        0.5
        +
        (3.0/16.0)
        /
        (tau_p-0.5);



    omega_m =
        1.0/tau_m;



    oldUx.resize(nx*ny,0.0);

    oldUy.resize(nx*ny,0.0);

}



//==================================================
// Initialization
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
// One simulation step
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
                ux-oldUx[id];


            double duy =
                uy-oldUy[id];



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
// Export Ghia horizontal profile
//
// u_x along y=0.5
//
//==================================================

void Cavity::exportCenterVelocityProfile(
        const std::string& filename)
{

    std::ofstream file(filename);



    file
        << "# Ghia benchmark\n";

    file
        << "# x_index   X   Ux/U_lid\n";



    int y =
        ny/2;



    file
        << std::setprecision(12);



    for(int x=0;x<nx;x++)
    {


        double X =
            static_cast<double>(x)
            /
            static_cast<double>(nx-1);



        double ux =
            lattice.getUx(x,y)
            /
            U_lid;



        file
            << x
            << " "
            << X
            << " "
            << ux
            << "\n";

    }


    file.close();

}



//==================================================
// Export Ghia vertical profile
//
// u_y along x=0.5
//
//==================================================

void Cavity::exportVerticalVelocityProfile(
        const std::string& filename)
{

    std::ofstream file(filename);



    file
        << "# Ghia benchmark\n";

    file
        << "# y_index   Y   Uy/U_lid\n";



    int x =
        nx/2;



    file
        << std::setprecision(12);



    for(int y=0;y<ny;y++)
    {


        double Y =
            static_cast<double>(y)
            /
            static_cast<double>(ny-1);



        double uy =
            lattice.getUy(x,y)
            /
            U_lid;



        file
            << y
            << " "
            << Y
            << " "
            << uy
            << "\n";

    }


    file.close();

}
