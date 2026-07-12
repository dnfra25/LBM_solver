#ifndef CAVITY_HPP
#define CAVITY_HPP


#include "Lattice.hpp"
#include "Boundary.hpp"

#include <vector>
#include <string>


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



    //=========================================
    // Ghia benchmark output
    //=========================================

    void writeCenterlineUx(
        const std::string& filename
    );


    void writeCenterlineUy(
        const std::string& filename
    );



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
