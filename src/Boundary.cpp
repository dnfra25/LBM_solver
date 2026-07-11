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


    for(int x = 1; x < nx-1; ++x)
    {

        int y = ny-1;


        //--------------------------------------
        // Density reconstruction
        //--------------------------------------

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


        //--------------------------------------
        // Zou-He moving wall
        //
        // D2Q9 convention:
        //
        // 0 rest
        // 1 east
        // 2 west
        // 3 north
        // 4 south
        // 5 NE
        // 6 SW
        // 7 NW
        // 8 SE
        //
        //--------------------------------------


        // south population

        f[lattice.index(4,x,y)] =
            f[lattice.index(3,x,y)]
            -
            (2.0/3.0)
            *
            rho
            *
            0.0;



        // south-west

        f[lattice.index(6,x,y)] =
            f[lattice.index(8,x,y)]
            +
            0.5 *
            (
                f[lattice.index(2,x,y)]
              -
                f[lattice.index(1,x,y)]
            )
            -
            (1.0/6.0)
            *
            rho
            *
            lidVelocity;



        // south-east

        f[lattice.index(8,x,y)] =
            f[lattice.index(6,x,y)]
            +
            0.5 *
            (
                f[lattice.index(1,x,y)]
              -
                f[lattice.index(2,x,y)]
            )
            +
            (1.0/6.0)
            *
            rho
            *
            lidVelocity;

    }

}
