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
    // Test cases
    //----------------------------------------------

    std::vector<TestCase> tests =
    {

        {
            65,
            65,
            100.0,
            0.1,
            "test_1_Re100_65"
        },


        {
            129,
            129,
            100.0,
            0.1,
            "test_2_Re100_129"
        },


        {
            129,
            129,
            400.0,
            0.1,
            "test_3_Re400_129"
        },


        {
            257,
            257,
            100.0,
            0.1,
            "test_4_Re100_257"
        }

    };



    //----------------------------------------------
    // Parameters
    //----------------------------------------------

    const int maxIterations = 20000;

    const double tolerance = 1e-10;



    //----------------------------------------------
    // Loop over tests
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



        Cavity cavity(
            test.nx,
            test.ny,
            test.Re,
            test.U
        );



        cavity.initialize();



        std::string filename =
            "convergence_test_"
            +
            std::to_string(testNumber)
            +
            ".dat";



        std::ofstream file(filename);



        file
            << "# iteration error\n";



        auto start =
            std::chrono::high_resolution_clock::now();



        int iteration;



        for(iteration=0;
            iteration<maxIterations;
            iteration++)
        {


            cavity.step();



            if(iteration % 50 == 0)
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




        //------------------------------------------
        // Center value
        //------------------------------------------

        Lattice& lattice =
            cavity.getLattice();



        double ux =
            lattice.getUx(
                test.nx/2,
                test.ny/2
            );



        double uy =
            lattice.getUy(
                test.nx/2,
                test.ny/2
            );



        std::cout
            << "\nFinal result\n";


        std::cout
            << "Iterations: "
            << iteration
            << "\n";


        std::cout
            << "Center Ux = "
            << std::setprecision(12)
            << ux
            << "\n";


        std::cout
            << "Center Uy = "
            << uy
            << "\n";


        std::cout
            << "Time = "
            << time
            << " s\n";



        file.close();



        testNumber++;

    }



    std::cout
        << "\nAll tests completed\n";


    return 0;

}
