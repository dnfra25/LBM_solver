#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "Lattice.hpp"


class Boundary
{

public:

    //=========================================
    // Constructor
    //=========================================

    Boundary(int nx,
             int ny,
             double lidVelocity);


    //=========================================
    // Apply all boundary conditions
    //=========================================

    void apply(Lattice& lattice);



private:

    int nx;
    int ny;


    // velocity of moving lid

    double lidVelocity;



    //=========================================
    // Boundary conditions
    //=========================================

    void bounceBackBottom(Lattice& lattice);

    void bounceBackLeft(Lattice& lattice);

    void bounceBackRight(Lattice& lattice);


    void movingTop(Lattice& lattice);


};

#endif
