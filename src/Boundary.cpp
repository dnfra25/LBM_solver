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


    bottomLeftCorner(lattice);

    bottomRightCorner(lattice);

    topLeftCorner(lattice);

    topRightCorner(lattice);

}



//==================================================
// Bottom wall
//==================================================

void Boundary::bounceBackBottom(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int y = 0;


    #pragma omp parallel for
    for(int x=1; x<nx-1; ++x)
    {

        f[lattice.index(3,x,y)] =
            f[lattice.index(4,x,y)];


        f[lattice.index(5,x,y)] =
            f[lattice.index(7,x,y)];


        f[lattice.index(8,x,y)] =
            f[lattice.index(6,x,y)];

    }

}



//==================================================
// Left wall
//==================================================

void Boundary::bounceBackLeft(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x = 0;


    #pragma omp parallel for
    for(int y=1; y<ny-1; ++y)
    {

        f[lattice.index(1,x,y)] =
            f[lattice.index(2,x,y)];


        f[lattice.index(5,x,y)] =
            f[lattice.index(6,x,y)];


        f[lattice.index(8,x,y)] =
            f[lattice.index(7,x,y)];

    }

}



//==================================================
// Right wall
//==================================================

void Boundary::bounceBackRight(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x = nx-1;


    #pragma omp parallel for
    for(int y=1; y<ny-1; ++y)
    {

        f[lattice.index(2,x,y)] =
            f[lattice.index(1,x,y)];


        f[lattice.index(6,x,y)] =
            f[lattice.index(5,x,y)];


        f[lattice.index(7,x,y)] =
            f[lattice.index(8,x,y)];

    }

}



//==================================================
// Moving top wall - Zou He velocity BC
//==================================================

void Boundary::movingTop(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int y = ny-1;


    #pragma omp parallel for
    for(int x=1; x<nx-1; ++x)
    {

        double rho =
              f[lattice.index(0,x,y)]
            + f[lattice.index(1,x,y)]
            + f[lattice.index(2,x,y)]
            + 2.0 *
              (
                f[lattice.index(3,x,y)]
              + f[lattice.index(5,x,y)]
              + f[lattice.index(7,x,y)]
              );



        // south population

        f[lattice.index(4,x,y)]
        =
        f[lattice.index(3,x,y)];



        // south-east

        f[lattice.index(8,x,y)]
        =
        f[lattice.index(7,x,y)]
        -
        0.5 *
        (
            f[lattice.index(1,x,y)]
          -
            f[lattice.index(2,x,y)]
        )
        +
        rho*lidVelocity/6.0;



        // south-west

        f[lattice.index(6,x,y)]
        =
        f[lattice.index(5,x,y)]
        +
        0.5 *
        (
            f[lattice.index(1,x,y)]
          -
            f[lattice.index(2,x,y)]
        )
        -
        rho*lidVelocity/6.0;


    }

}



//==================================================
// Bottom-left corner
//==================================================

void Boundary::bottomLeftCorner(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x=0;
    int y=0;


    f[lattice.index(1,x,y)] =
        f[lattice.index(2,x,y)];


    f[lattice.index(3,x,y)] =
        f[lattice.index(4,x,y)];


    f[lattice.index(5,x,y)] =
        f[lattice.index(7,x,y)];


    f[lattice.index(8,x,y)] =
        f[lattice.index(6,x,y)];

}



//==================================================
// Bottom-right corner
//==================================================

void Boundary::bottomRightCorner(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x=nx-1;
    int y=0;


    f[lattice.index(2,x,y)] =
        f[lattice.index(1,x,y)];


    f[lattice.index(3,x,y)] =
        f[lattice.index(4,x,y)];


    f[lattice.index(6,x,y)] =
        f[lattice.index(8,x,y)];


    f[lattice.index(7,x,y)] =
        f[lattice.index(5,x,y)];

}



//==================================================
// Top-left corner
//==================================================

void Boundary::topLeftCorner(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x=0;
    int y=ny-1;


    f[lattice.index(1,x,y)] =
        f[lattice.index(2,x,y)];


    f[lattice.index(4,x,y)] =
        f[lattice.index(3,x,y)];


    f[lattice.index(8,x,y)] =
        f[lattice.index(7,x,y)];

}



//==================================================
// Top-right corner
//==================================================

void Boundary::topRightCorner(Lattice& lattice)
{

    auto& f = lattice.distributions();


    int x=nx-1;
    int y=ny-1;


    f[lattice.index(2,x,y)] =
        f[lattice.index(1,x,y)];


    f[lattice.index(4,x,y)] =
        f[lattice.index(3,x,y)];


    f[lattice.index(6,x,y)] =
        f[lattice.index(5,x,y)];

}
