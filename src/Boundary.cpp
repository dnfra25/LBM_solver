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

void Boundary::movingTop(Lattice& lattice)
{

    auto& f = lattice.distributions();


    #pragma omp parallel for
    for(int x=1;x<nx-1;x++)
    {

        int y = ny-1;



        /*
            Known populations:
            
            q=0
            q=1
            q=2
            q=4
            q=6
            q=8


            Unknown populations:

            q=3
            q=5
            q=7

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



        if(rho < 1e-12)
            rho = 1.0;



        double ux = lidVelocity;



        double uy = 0.0;



        /*
            Zou-He reconstruction
        */


        // q=3 north

        f[lattice.index(3,x,y)] =
            f[lattice.index(4,x,y)]
            +
            (2.0/3.0)
            *
            rho
            *
            uy;



        // q=5 north-east

        f[lattice.index(5,x,y)] =
            f[lattice.index(7,x,y)]
            +
            0.5 *
            (
                f[lattice.index(2,x,y)]
              -
                f[lattice.index(1,x,y)]
            )
            +
            (1.0/6.0)
            *
            rho
            *
            ux;



        // q=7 north-west

        f[lattice.index(7,x,y)] =
            f[lattice.index(5,x,y)]
            -
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
            ux;

    }

}
