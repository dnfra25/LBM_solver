#include "Boundary.hpp"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif



//==================================================
// Constructor
//==================================================

Boundary::Boundary(int nx_,
                   int ny_,
                   double lidVelocity_)
:
nx(nx_),
ny(ny_),
lidVelocity(lidVelocity_)
{

}



//==================================================
// Apply all boundary conditions
//==================================================

void Boundary::apply(Lattice& lattice)
{

    bounceBackBottom(lattice);

    bounceBackLeft(lattice);

    bounceBackRight(lattice);

    movingTop(lattice);

}



//==================================================
// Bottom wall - no slip bounce back
//==================================================

void Boundary::bounceBackBottom(Lattice& lattice)
{

    auto& f = lattice.distributions();


    #pragma omp parallel for
    for(int x=0;x<nx;x++)
    {

        int y=0;


        f[lattice.index(3,x,y)] =
            f[lattice.index(4,x,y)];


        f[lattice.index(5,x,y)] =
            f[lattice.index(7,x,y)];


        f[lattice.index(8,x,y)] =
            f[lattice.index(6,x,y)];

    }

}



//==================================================
// Left wall - no slip bounce back
//==================================================

void Boundary::bounceBackLeft(Lattice& lattice)
{

    auto& f = lattice.distributions();


    #pragma omp parallel for
    for(int y=0;y<ny;y++)
    {

        int x=0;


        f[lattice.index(1,x,y)] =
            f[lattice.index(2,x,y)];


        f[lattice.index(5,x,y)] =
            f[lattice.index(6,x,y)];


        f[lattice.index(8,x,y)] =
            f[lattice.index(7,x,y)];

    }

}



//==================================================
// Right wall - no slip bounce back
//==================================================

void Boundary::bounceBackRight(Lattice& lattice)
{

    auto& f = lattice.distributions();


    #pragma omp parallel for
    for(int y=0;y<ny;y++)
    {

        int x=nx-1;


        f[lattice.index(2,x,y)] =
            f[lattice.index(1,x,y)];


        f[lattice.index(6,x,y)] =
            f[lattice.index(5,x,y)];


        f[lattice.index(7,x,y)] =
            f[lattice.index(8,x,y)];

    }

}



//==================================================
// Moving top lid - Zou He velocity boundary
//==================================================

//==================================================
// Moving lid - Zou-He velocity boundary
//==================================================

void Boundary::movingTop(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int y = ny - 1;


    #pragma omp parallel for
    for(int x=1; x<nx-1; ++x)
    {

        /*
            Top wall:

            known:
                f0
                f1
                f2
                f4
                f6
                f8

            unknown:
                f3
                f5
                f7

            ux = lidVelocity
            uy = 0
        */


        double rho =
              f[lattice.index(0,x,y)]
            + f[lattice.index(1,x,y)]
            + f[lattice.index(2,x,y)]
            + 2.0 *
              (
                 f[lattice.index(4,x,y)]
               + f[lattice.index(6,x,y)]
               + f[lattice.index(8,x,y)]
              );


        // north population
        f[lattice.index(3,x,y)]
        =
        f[lattice.index(4,x,y)];


        // north-east
        f[lattice.index(5,x,y)]
        =
        f[lattice.index(6,x,y)]
        +
        0.5*
        (
            f[lattice.index(2,x,y)]
          -
            f[lattice.index(1,x,y)]
        )
        +
        (rho*lidVelocity)/6.0;



        // north-west
        f[lattice.index(7,x,y)]
        =
        f[lattice.index(8,x,y)]
        +
        0.5*
        (
            f[lattice.index(1,x,y)]
          -
            f[lattice.index(2,x,y)]
        )
        -
        (rho*lidVelocity)/6.0;


    }

}
