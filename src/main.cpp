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
    // Ghia benchmark cases -> Faccio dei test allo scopo di confrontarmi con Ghia
    // ogni run prende in input i dati per la griglia, RE, U_LID, e un nome da dare al test
    //------------------------------------------------

    std::vector<TestCase> tests =
    {

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

        {
            129,
            129,
            10000.0,
            0.05,
            "Re10000"
        }

    };

    //------------------------------------------------
    // Parameetri di controllo
    //------------------------------------------------

    const int maxIterations =
        20000;


    const double tolerance =
        1.0e-8;



    //------------------------------------------------
    // Run cases: stampo output parlanti
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
        // Create cavity: gestisco FISICAMENTE la cavità 2d
        //------------------------------------------------

        Cavity cavity(
            test.nx,
            test.ny,
            test.Re,
            test.U
        );

        cavity.initialize();

        //------------------------------------------------
        // Time measurement: così ho le tempistiche
        //------------------------------------------------

        auto start =
            std::chrono::high_resolution_clock::now();

        //------------------------------------------------
        // LBM loop: gestisco i loop del modello
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
        // Timing: formatto l'orario
        //------------------------------------------------

        auto end =
            std::chrono::high_resolution_clock::now();
        double time =
            std::chrono::duration<double>
            (end-start)
            .count();


        //------------------------------------------------
        // Output Ghia profiles: stampo i file per i risultati
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
        // Diagnostics: stampo info utili per l'output a video
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
            << lattice.getUx(cx,cy)
            << "\n";

        std::cout
            << "Center Uy = "
            << lattice.getUy(cx,cy)
            << "\n";

        std::cout
            << "Top interior Ux = "
            << lattice.getUx(cx,test.ny-2)
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
