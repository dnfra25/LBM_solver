#include "Cavity.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

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



int main()
{


#ifdef _OPENMP

    std::cout
        << "OpenMP threads: "
        << omp_get_max_threads()
        << "\n\n";

#endif



    //----------------------------------------------
    // Ghia benchmark test cases
    //----------------------------------------------

    std::vector<TestCase> tests =
    {

        {
            129,
            129,
            100.0,
            0.1,
            "Re100"
        },


        {
            129,
            129,
            400.0,
            0.1,
            "Re400"
        },


        {
            129,
            129,
            1000.0,
            0.1,
            "Re1000"
        },


        {
            129,
            129,
            3200.0,
            0.1,
            "Re3200"
        },


        {
            129,
            129,
            5000.0,
            0.1,
            "Re5000"
        },


        {
            129,
            129,
            7500.0,
            0.1,
            "Re7500"
        },


        {
            129,
            129,
            10000.0,
            0.1,
            "Re10000"
        }

    };




    //----------------------------------------------
    // Simulation parameters
    //----------------------------------------------

    const int maxIterations = 100000;

    const double tolerance = 1e-8;




    //----------------------------------------------
    // Run simulations
    //----------------------------------------------

    int testNumber = 1;



    for(auto& test : tests)
    {


        std::cout
            << "\n=====================\n";


        std::cout
            << "Running "
            << test.name
            << "\n";




        Cavity cavity
        (
            test.nx,
            test.ny,
            test.Re,
            test.U
        );




        //------------------------------------------
        // Initialization
        //------------------------------------------

        cavity.initialize();




        //------------------------------------------
        // Check top boundary
        //------------------------------------------

        cavity.applyBoundary();

        cavity.getLattice()
               .computeMacroscopic();



        std::cout
            << "Initial lid Ux = "
            << cavity.getLattice()
                    .getUx(
                        test.nx/2,
                        test.ny-1
                    )
            << "\n";




        //------------------------------------------
        // Reset
        //------------------------------------------

        cavity.initialize();




        //------------------------------------------
        // Convergence file
        //------------------------------------------

        std::string filename =
            "convergence_Re"
            +
            std::to_string((int)test.Re)
            +
            ".dat";



        std::ofstream file(filename);



        file
            << "# iteration error\n";




        //------------------------------------------
        // Timer
        //------------------------------------------

        auto start =
            std::chrono::high_resolution_clock::now();




        int iteration;




        //------------------------------------------
        // Main LBM loop
        //------------------------------------------

        for(iteration=0;
            iteration<maxIterations;
            iteration++)
        {


            cavity.step();



            if(iteration % 100 == 0 &&
               iteration > 0)
            {


                double error =
                    cavity.velocityDifference();



                file
                    << iteration
                    << " "
                    << std::scientific
                    << error
                    << "\n";



                std::cout
                    << "it = "
                    << iteration
                    << " error = "
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




        auto end =
            std::chrono::high_resolution_clock::now();



        double time =
            std::chrono::duration<double>
            (end-start)
            .count();




        file.close();




        //------------------------------------------
        // Export Ghia profiles
        //------------------------------------------

        std::string ux_file =
            "Ux_center_Re"
            +
            std::to_string((int)test.Re)
            +
            ".dat";



        std::string uy_file =
            "Uy_center_Re"
            +
            std::to_string((int)test.Re)
            +
            ".dat";



        cavity.exportCenterVelocityProfile(
            ux_file
        );



        cavity.exportVerticalVelocityProfile(
            uy_file
        );





        //------------------------------------------
        // Diagnostics
        //------------------------------------------

        Lattice& lattice =
            cavity.getLattice();



        int cx =
            test.nx/2;


        int cy =
            test.ny/2;



        double centerUx =
            lattice.getUx(cx,cy);


        double centerUy =
            lattice.getUy(cx,cy);



        double topUx =
            lattice.getUx(
                cx,
                test.ny-2
            );



        double bottomUx =
            lattice.getUx(
                cx,
                1
            );




        std::cout
            << "\nFinal result\n";



        std::cout
            << "Iterations: "
            << iteration
            << "\n";



        std::cout
            << std::setprecision(12);



        std::cout
            << "Center Ux = "
            << centerUx
            << "\n";



        std::cout
            << "Center Uy = "
            << centerUy
            << "\n";



        std::cout
            << "Top interior Ux = "
            << topUx
            << "\n";



        std::cout
            << "Bottom interior Ux = "
            << bottomUx
            << "\n";



        std::cout
            << "Time = "
            << time
            << " s\n";



        testNumber++;

    }




    std::cout
        << "\nAll Ghia benchmark tests completed\n";



    return 0;

}
