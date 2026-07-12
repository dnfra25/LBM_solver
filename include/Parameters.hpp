#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP


namespace Parameters
{

constexpr int Nx = 200;
constexpr int Ny = 200;


constexpr double U_lid = 0.1;


constexpr double L = Nx - 1;


constexpr double Re = 100.0;



// LBM

constexpr double cs2 = 1.0/3.0;


constexpr double nu =
    U_lid * L / Re;



constexpr double tau =
    0.5 + nu/cs2;



// TRT

constexpr double lambda_trt = 0.25;


constexpr double tau_p =
    tau;


constexpr double tau_m =
    lambda_trt/(tau_p-0.5)+0.5;


constexpr double omega_p =
    1.0/tau_p;


constexpr double omega_m =
    1.0/tau_m;



constexpr int max_iterations = 50000;


constexpr int output_interval = 500;


constexpr double tolerance = 1.0e-8;


constexpr double rho0 = 1.0;


}

#endif
