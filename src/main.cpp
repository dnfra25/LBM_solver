#include "Cavity.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif



struct TestCase
{
    int nx;
    int ny;

    double Re;

    double U;

    std::string name;
};




//==================================================
// Main
//==================================================

int main()
{

#ifdef _OPENMP

    std::cout
        << "OpenMP threads: "
        << omp_get_max_threads()
        << "\n\n";

#endif



    //------------------------------------------------
    // Ghia benchmark cases
    //------------------------------------------------

    std::vector<TestCase> tests =
    {
/*
        {
            129,
            129,
            100.0,
            0.05,
            "Re100"
        },


        {
            129,
            129,
            400.0,
            0.05,
            "Re400"
        },


        {
            129,
            129,
            1000.0,
            0.05,
            "Re1000"
        },


        {
            129,
            129,
            3200.0,
            0.05,
            "Re3200"
        },


        {
            129,
            129,
            5000.0,
            0.05,
            "Re5000"
        },


        {
            129,
            129,
            7500.0,
            0.05,
            "Re7500"
        },
*/

        {
            129,
            129,
            10000.0,
            0.05,
            "Re10000"
        }

    };



    //------------------------------------------------
    // Simulation parameters
    //------------------------------------------------

    const int maxIterations =
        100;


    const double tolerance =
        1.0e-8;



    //------------------------------------------------
    // Run cases
    //------------------------------------------------

    for(auto& test : tests)
    {


        std::cout
            << "\n=============================\n";


        std::cout
            << "Running "
            << test.name
            << "\n";



        //------------------------------------------------
        // Create cavity
        //------------------------------------------------

        Cavity cavity(
            test.nx,
            test.ny,
            test.Re,
            test.U
        );



        cavity.initialize();



        //------------------------------------------------
        // Time measurement
        //------------------------------------------------

        auto start =
            std::chrono::high_resolution_clock::now();



        //------------------------------------------------
        // LBM loop
        //------------------------------------------------

        int iteration;


        double error = 1.0;



        for(iteration=0;
            iteration<maxIterations;
            iteration++)
        {


            cavity.step();



            if(iteration % 100 == 0)
            {

                error =
                    cavity.velocityDifference();



                std::cout
                    << "it = "
                    << iteration
                    << " error = "
                    << std::scientific
                    << error
                    << "\n";



                if(error < tolerance)
                {

                    std::cout
                        << "Converged\n";


                    break;

                }

            }


        }



        //------------------------------------------------
        // Timing
        //------------------------------------------------

        auto end =
            std::chrono::high_resolution_clock::now();



        double time =
            std::chrono::duration<double>
            (end-start)
            .count();



        //------------------------------------------------
        // Output Ghia profiles
        //------------------------------------------------


        std::string uxFile =
            "Ux_" + test.name + ".dat";


        std::string uyFile =
            "Uy_" + test.name + ".dat";



        cavity.writeCenterlineUx(
            uxFile
        );


        cavity.writeCenterlineUy(
            uyFile
        );



        //------------------------------------------------
        // Diagnostics
        //------------------------------------------------

        Lattice& lattice =
            cavity.getLattice();



        int cx =
            test.nx/2;


        int cy =
            test.ny/2;



        std::cout
            << "\nFinal diagnostics\n";


        std::cout
            << std::setprecision(12);



        std::cout
            << "Center Ux = "
            << lattice.getUx(cx,cy)/test.U
            << "\n";


        std::cout
            << "Center Uy = "
            << lattice.getUy(cx,cy)/test.U
            << "\n";


        std::cout
            << "Top interior Ux = "
            << lattice.getUx(cx,test.ny-2)/test.U
            << "\n";


        std::cout
            << "Iterations = "
            << iteration
            << "\n";


        std::cout
            << "Time = "
            << time
            << " s\n";



        std::cout
            << "Written:\n";


        std::cout
            << "  "
            << uxFile
            << "\n";


        std::cout
            << "  "
            << uyFile
            << "\n";


    }



    std::cout
        << "\nAll Ghia benchmark cases completed\n";



    return 0;

}
