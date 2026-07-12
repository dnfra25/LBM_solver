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

    f.resize(Q * size, 0.0);
    f_eq.resize(Q * size, 0.0);
    f_post.resize(Q * size, 0.0);

    rho.resize(size, rho0);
    ux.resize(size, 0.0);
    uy.resize(size, 0.0);

    // D2Q9 velocities

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

    // D2Q9 weights

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

    // Opposite directions

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

            double ux0 = 0.0;
            double uy0 = 0.0;


            ux[id] = ux0;
            uy[id] = uy0;



            double u2 =
                ux0*ux0
              +
                uy0*uy0;



            for(int q = 0; q < Q; ++q)
            {

                double cu =
                    ux0*c[q][0]
                  +
                    uy0*c[q][1];



                double feq =
                    rho0
                    *
                    w[q]
                    *
                    (
                        1.0
                        +
                        3.0*cu
                        +
                        4.5*cu*cu
                        -
                        1.5*u2
                    );



                f[index(q,x,y)] =
                    feq;


                f_eq[index(q,x,y)] =
                    feq;


                f_post[index(q,x,y)] =
                    feq;

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
    for(int x=0;x<nx;x++)
    {
        for(int y=0;y<ny;y++)
        {

            double density=0.0;

            double vx=0.0;

            double vy=0.0;



            for(int q=0;q<Q;q++)
            {

                double fq =
                    f[index(q,x,y)];


                density += fq;


                vx += fq*c[q][0];

                vy += fq*c[q][1];

            }



            int id =
                cell(x,y);



            rho[id]=density;



            if(density>1e-14)
            {
                ux[id]=vx/density;
                uy[id]=vy/density;
            }
            else
            {
                ux[id]=0.0;
                uy[id]=0.0;
            }

        }
    }

}



//==================================================
// Compute equilibrium distribution
//==================================================

void Lattice::computeEquilibrium()
{

    #pragma omp parallel for collapse(2)
    for(int x=0;x<nx;x++)
    {
        for(int y=0;y<ny;y++)
        {

            int id =
                cell(x,y);



            double u2 =
                ux[id]*ux[id]
              +
                uy[id]*uy[id];



            for(int q=0;q<Q;q++)
            {

                double cu =
                    ux[id]*c[q][0]
                  +
                    uy[id]*c[q][1];



                double feq =
                    1.0
                    +
                    3.0*cu
                    +
                    4.5*cu*cu
                    -
                    1.5*u2;



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
// TRT collision
//==================================================

void Lattice::collision()
{

    #pragma omp parallel for collapse(2)
    for(int x = 0; x < nx; ++x)
    {
        for(int y = 0; y < ny; ++y)
        {

            for(int q = 0; q < Q; ++q)
            {

                int qb = opposite[q];


                // evita di processare due volte la stessa coppia
                if(q > qb)
                    continue;


                double fq =
                    f[index(q,x,y)];

                double fqb =
                    f[index(qb,x,y)];


                double feq =
                    f_eq[index(q,x,y)];

                double feqb =
                    f_eq[index(qb,x,y)];



                //--------------------------------------
                // symmetric part
                //--------------------------------------

                double fp =
                    0.5*(fq + fqb);


                double feqp =
                    0.5*(feq + feqb);



                //--------------------------------------
                // antisymmetric part
                //--------------------------------------

                double fm =
                    0.5*(fq - fqb);


                double feqm =
                    0.5*(feq - feqb);



                //--------------------------------------
                // TRT relaxation
                //--------------------------------------

                double fp_new =
                    fp
                    -
                    omega_p*(fp-feqp);



                double fm_new =
                    fm
                    -
                    omega_m*(fm-feqm);



                //--------------------------------------
                // reconstruction
                //--------------------------------------

                double fq_new =
                    fp_new + fm_new;


                double fqb_new =
                    fp_new - fm_new;



                f_post[index(q,x,y)]
                    =
                    fq_new;


                f_post[index(qb,x,y)]
                    =
                    fqb_new;

            }

        }
    }

}



//==================================================
// Streaming - pull scheme
//==================================================
//==================================================
// Streaming - pull scheme
//==================================================

void Lattice::streaming()
{

    std::vector<double> f_new(
        Q * nx * ny,
        0.0
    );


    #pragma omp parallel for collapse(2)
    for(int x=0; x<nx; ++x)
    {
        for(int y=0; y<ny; ++y)
        {

            for(int q=0; q<Q; ++q)
            {

                int xs =
                    x - c[q][0];

                int ys =
                    y - c[q][1];


                if(xs >= 0 &&
                   xs < nx &&
                   ys >=0 &&
                   ys < ny)
                {

                    f_new[index(q,x,y)]
                    =
                    f_post[index(q,xs,ys)];

                }
                else
                {

                    /*
                     lascia vuoto.
                     La boundary ricostruisce
                     le popolazioni mancanti.
                    */

                    f_new[index(q,x,y)]
                    =
                    0.0;

                }

            }

        }
    }


    f.swap(f_new);

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
