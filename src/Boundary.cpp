#include "Boundary.hpp"
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

//==================================================
// Constructor: il contorno fisico della cavity 2d gestito qui
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
// In sequenza applico tutti i limiti al contorno
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
// Parete inferiore
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
// PArete a sx
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
// Parete a DX
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
              f[lattice.index(4,x,y)]
            + f[lattice.index(6,x,y)]
            + f[lattice.index(8,x,y)]
            );


        // North
        f[lattice.index(3,x,y)]
        =
        f[lattice.index(4,x,y)];


        // North-East
        f[lattice.index(5,x,y)]
        =
        f[lattice.index(6,x,y)]
        +
        rho*lidVelocity/6.0;


        // North-West
        f[lattice.index(7,x,y)]
        =
        f[lattice.index(8,x,y)]
        -
        rho*lidVelocity/6.0;

    }

}


//==================================================
// Angolo in basso a sx
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
// Angolo in basso a dx
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
// Angolo in alto a sx
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
// Angolo in alto a dx
//==================================================

void Boundary::topRightCorner(Lattice& lattice)
{

    auto& f = lattice.distributions();

    int x = nx-1;
    int y = ny-1;

    f[lattice.index(1,x,y)]
    =
    f[lattice.index(2,x,y)];

    f[lattice.index(3,x,y)]
    =
    f[lattice.index(4,x,y)];

    f[lattice.index(5,x,y)]
    =
    f[lattice.index(6,x,y)];
}
