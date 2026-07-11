#include "Lattice.hpp"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

//==================================================
// Constructor
//==================================================
Lattice::Lattice(int nx_,
                 int ny_,
                 double omega_p_,
                 double omega_m_,
                 double rho0_)
:
nx(nx_),
ny(ny_),
omega_p(omega_p_),
omega_m(omega_m_),
rho0(rho0_)
{
    const int size = nx * ny;
    //----------------------------------------------
    // Allocate memory
    //----------------------------------------------

    f.resize(Q * size, 0.0);
    f_eq.resize(Q * size, 0.0);
    f_post.resize(Q * size, 0.0);
    rho.resize(size, rho0);
    ux.resize(size, 0.0);
    uy.resize(size, 0.0);

    //----------------------------------------------
    // D2Q9 lattice velocities
    //----------------------------------------------

    c =
    {{
        {{ 0, 0 }},

        {{ 1, 0 }},
        {{-1, 0 }},

        {{ 0, 1 }},
        {{ 0,-1 }},

        {{ 1, 1 }},
        {{-1,-1 }},

        {{-1, 1 }},
        {{ 1,-1 }}
    }};

    //----------------------------------------------
    // D2Q9 weights
    //----------------------------------------------

    w =
    {{
        4.0/9.0,

        1.0/9.0,
        1.0/9.0,

        1.0/9.0,
        1.0/9.0,

        1.0/36.0,
        1.0/36.0,

        1.0/36.0,
        1.0/36.0
    }};

    //----------------------------------------------
    // Opposite directions
    //----------------------------------------------

    opposite =
    {{
        0,
        2,
        1,
        4,
        3,
        6,
        5,
        8,
        7
    }};

}


//==================================================
// Initialize distributions
//==================================================

void Lattice::initialize()
{

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            int id = cell(x,y);
            rho[id] = rho0;
            ux[id] = 0.0;
            uy[id] = 0.0;

            for(int q = 0; q < Q; ++q)
            {

                double equilibrium =
                    rho0 * w[q];


                f[index(q,x,y)] =
                    equilibrium;

                f_eq[index(q,x,y)] =
                    equilibrium;

                f_post[index(q,x,y)] =
                    equilibrium;

            }

        }
    }

}



//==================================================
// Compute macroscopic variables
//==================================================

void Lattice::computeMacroscopic()
{

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            double density = 0.0;

            double vx = 0.0;

            double vy = 0.0;



            for(int q = 0; q < Q; ++q)
            {

                double fq =
                    f[index(q,x,y)];


                density += fq;


                vx += fq * c[q][0];

                vy += fq * c[q][1];

            }



            int id = cell(x,y);


            rho[id] = density;


            ux[id] = vx / density;

            uy[id] = vy / density;

        }
    }

}


//==================================================
// Compute equilibrium distribution
//==================================================

void Lattice::computeEquilibrium()
{

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            int id = cell(x,y);


            double u2 =
                ux[id]*ux[id]
              +
                uy[id]*uy[id];



            for(int q = 0; q < Q; ++q)
            {

                double cu =
                    ux[id]*c[q][0]
                  +
                    uy[id]*c[q][1];


                double feq =
                    1.0
                    + 3.0*cu
                    + 4.5*cu*cu
                    - 1.5*u2;



                f_eq[index(q,x,y)] =
                    rho[id]
                    *
                    w[q]
                    *
                    feq;

            }

        }
    }

}
//==================================================
// TRT collision step
//==================================================

void Lattice::collision()
{

    /*
    
    Equivalent to:

        nb_col_str()

    from the original Python repository.

    Here we compute:

        f_post

    from:

        f
        f_eq

    */

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            for(int q = 0; q < Q; ++q)
            {

                int qb = opposite[q];


                double fq =
                    f[index(q,x,y)];


                double fqb =
                    f[index(qb,x,y)];


                double feq =
                    f_eq[index(q,x,y)];


                double feqb =
                    f_eq[index(qb,x,y)];



                double sum =
                    omega_p + omega_m;


                double diff =
                    omega_p - omega_m;



                f_post[index(q,x,y)] =

                    (1.0 - 0.5*sum)*fq

                    -

                    0.5*diff*fqb

                    +

                    0.5*sum*feq

                    +

                    0.5*diff*feqb;

            }

        }
    }

}






//==================================================
// Streaming step
//==================================================

void Lattice::streaming()
{

    /*
    
    Equivalent to:

        g <- streamed g_up

    in the Python implementation.

    */


    //----------------------------------------------
    // Rest population
    //----------------------------------------------

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            f[index(0,x,y)] =
                f_post[index(0,x,y)];

        }
    }



    //----------------------------------------------
    // Moving populations
    //----------------------------------------------

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {


            // q=1 : east

            if(x+1 < nx)
            {
                f[index(1,x+1,y)] =
                    f_post[index(1,x,y)];
            }



            // q=2 : west

            if(x-1 >= 0)
            {
                f[index(2,x-1,y)] =
                    f_post[index(2,x,y)];
            }



            // q=3 : north

            if(y+1 < ny)
            {
                f[index(3,x,y+1)] =
                    f_post[index(3,x,y)];
            }



            // q=4 : south

            if(y-1 >= 0)
            {
                f[index(4,x,y-1)] =
                    f_post[index(4,x,y)];
            }



            // q=5 : north-east

            if((x+1 < nx) && (y+1 < ny))
            {
                f[index(5,x+1,y+1)] =
                    f_post[index(5,x,y)];
            }



            // q=6 : south-west

            if((x-1 >= 0) && (y-1 >= 0))
            {
                f[index(6,x-1,y-1)] =
                    f_post[index(6,x,y)];
            }



            // q=7 : north-west

            if((x-1 >= 0) && (y+1 < ny))
            {
                f[index(7,x-1,y+1)] =
                    f_post[index(7,x,y)];
            }



            // q=8 : south-east

            if((x+1 < nx) && (y-1 >= 0))
            {
                f[index(8,x+1,y-1)] =
                    f_post[index(8,x,y)];
            }

        }
    }

}






//==================================================
// Perform one LBM iteration
//==================================================

void Lattice::step()
{

    computeMacroscopic();

    computeEquilibrium();

    collision();

    streaming();

}






//==================================================
// Get density
//==================================================

double Lattice::getRho(int x,int y) const
{
    return rho[cell(x,y)];
}



//==================================================
// Get x velocity
//==================================================

double Lattice::getUx(int x,int y) const
{
    return ux[cell(x,y)];
}


//==================================================
// Get y velocity
//==================================================

double Lattice::getUy(int x,int y) const
{
    return uy[cell(x,y)];
}
