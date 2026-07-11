#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP


namespace Parameters
{

    // -----------------------------
    // Lattice size
    // -----------------------------

    constexpr int Nx = 200;
    constexpr int Ny = 200;


    // -----------------------------
    // Physical parameters
    // -----------------------------

    // Lid velocity
    constexpr double U_lid = 0.1;


    // Reference length
    // for cavity: L = Nx-1 lattice units
    constexpr double L = Nx - 1;


    // Reynolds number
    constexpr double Re = 100.0;



    // -----------------------------
    // LBM parameters
    // -----------------------------

    // Speed of sound squared
    constexpr double cs2 = 1.0/3.0;


    // Kinematic viscosity
    constexpr double nu =
        U_lid * L / Re;



    // Relaxation time
    // nu = cs2*(tau-0.5)

    constexpr double tau =
        0.5 + nu/cs2;



    // TRT parameters

    constexpr double lambda_trt = 0.25;


    constexpr double tau_p = tau;


    constexpr double tau_m =
        lambda_trt/(tau_p-0.5)+0.5;


    constexpr double omega_p =
        1.0/tau_p;


    constexpr double omega_m =
        1.0/tau_m;



    // -----------------------------
    // Simulation control
    // -----------------------------

    constexpr int max_iterations = 50000;


    // check convergence every N iterations
    constexpr int output_interval = 500;


    // convergence tolerance
    constexpr double tolerance = 1.0e-8;



    // Initial density

    constexpr double rho0 = 1.0;


}



#endif
