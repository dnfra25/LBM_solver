#ifndef CAVITY_HPP
#define CAVITY_HPP


#include "Lattice.hpp"
#include "Boundary.hpp"


class Cavity
{

public:

    //=========================================
    // Constructor
    //=========================================

    Cavity(int nx,
           int ny,
           double Reynolds,
           double lidVelocity);



    //=========================================
    // Initialization
    //=========================================

    void initialize();



    //=========================================
    // Perform one iteration
    //=========================================

    void step();



    //=========================================
    // Boundary application
    //=========================================

    void applyBoundary();



    //=========================================
    // Access lattice
    //=========================================

    Lattice& getLattice();



    //=========================================
    // Convergence
    //=========================================

    double velocityDifference();



    bool converged(double tolerance);



private:


    int nx;
    int ny;


    double Re;

    double U_lid;


    double viscosity;


    double omega_p;

    double omega_m;



    Lattice lattice;

    Boundary boundary;



    std::vector<double> oldUx;

    std::vector<double> oldUy;



};


#endif
