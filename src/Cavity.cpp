#include "Cavity.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif


//==================================================
// Constructor: gestisco la cavità FISICAMENTE con una classe
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
    U_lid * (nx_ - 1) / Re
),


omega_p(
    1.0 /
    (
        0.5
        +
        3.0 *
        viscosity
    )
),


// TRT omega_minus
// Lambda = 3/16 -> valore standard
omega_m(
    1.0 /
    (
        0.5
        +
        (3.0/16.0)
        /
        (
            (1.0/omega_p)
            -
            0.5
        )
    )
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

  std::cout
    << "Re = "
    << Re
    << "\n";

std::cout
    << "U_lid = "
    << U_lid
    << "\n";

std::cout
    << "nu = "
    << viscosity
    << "\n";

std::cout
    << "omega_plus = "
    << omega_p
    << "\n";

std::cout
    << "omega_minus = "
    << omega_m
    << "\n\n";
    oldUx.resize(
        nx*ny,
        0.0
    );

    oldUy.resize(
        nx*ny,
        0.0
    );

}

//==================================================
// Inizializzo un oggetto cavità
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
// Step standard per un test del modello
//==================================================

void Cavity::step()
{
// FASI IN SEQUENZA per il loro significato sifico
    lattice.computeMacroscopic();
    lattice.computeEquilibrium();
    lattice.collision();
    lattice.streaming();
    boundary.apply(lattice);
    lattice.computeMacroscopic();
}

//==================================================
// Applico le condizioni fisiche al contorno: quelle della cavità
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
// Convergence test -> Mi serve per controllare che l'errore di approx non esploda
//==================================================

bool Cavity::converged(double tolerance)
{

    double error =
        velocityDifference();

    return
        std::isfinite(error)
        &&
        error < tolerance;
}

//==================================================
// Focus su vertical centerline Ux -> Ghia
// u_x(x=0.5,y)
//==================================================

void Cavity::writeCenterlineUx(
        const std::string& filename)
{

    std::ofstream file(filename);


    file << "# y ux\n";


    // centro della cavità in x
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


        file
            << Y
            << " "
            << lattice.getUx(x,y)
            << "\n";

    }


    file.close();

}



//==================================================
// Focus su horizontal centerline Uy -> Ghia
// u_y(x,y=0.5)
//==================================================
void Cavity::writeCenterlineUy(
        const std::string& filename)
{
    std::ofstream file(filename);

    file << "# x uy\n";


    int y =
        (ny-1)/2;


    file
        << std::scientific
        << std::setprecision(12);


    for(int x=0;x<nx;x++)
    {

        double X =
            static_cast<double>(x)
            /
            static_cast<double>(nx-1);


        file
            << X
            << " "
            << lattice.getUy(x,y)
            << "\n";

    }


    file.close();
}
