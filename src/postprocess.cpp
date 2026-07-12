#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>


struct Profile
{
    std::vector<double> coord;
    std::vector<double> value;
};


//==================================================
// Read profile file
//==================================================

Profile readProfile(const std::string& filename)
{
    Profile profile;


    std::ifstream file(filename);


    if(!file)
    {
        std::cout
            << "Cannot open "
            << filename
            << std::endl;

        return profile;
    }



    std::string line;


    while(std::getline(file,line))
    {

        if(line.empty())
            continue;


        if(line[0]=='#')
            continue;



        std::stringstream ss(line);


        double x;
        double u;


        ss >> x >> u;


        profile.coord.push_back(x);

        profile.value.push_back(u);

    }


    file.close();


    return profile;
}



//==================================================
// Main
//==================================================

int main()
{

    //--------------------------------------------------
    // Reynolds cases
    //--------------------------------------------------

    std::vector<int> Reynolds =
    {
        100,
        400,
        1000,
        3200,
        5000,
        7500,
        10000
    };



    //--------------------------------------------------
    // Ghia grid points
    //--------------------------------------------------

    std::vector<int> gridPts =
    {
        129,
        125,
        124,
        123,
        122,
        117,
        111,
        104,
        65,
        31,
        30,
        21,
        13,
        11,
        10,
        9,
        1
    };



    //--------------------------------------------------
    // Read Ux and Uy profiles
    //--------------------------------------------------

    std::map<int,Profile> Ux;

    std::map<int,Profile> Uy;



    for(int Re : Reynolds)
    {

        std::string uxFile =
            "Ux_Re"
            +
            std::to_string(Re)
            +
            ".dat";


        std::string uyFile =
            "Uy_Re"
            +
            std::to_string(Re)
            +
            ".dat";



        Ux[Re] = readProfile(uxFile);

        Uy[Re] = readProfile(uyFile);

    }



    //--------------------------------------------------
    // Write Ux table
    //--------------------------------------------------

    std::ofstream fileUx("ghia_u.csv");



    fileUx
        << "gridpt,x";


    for(int Re : Reynolds)
    {
        fileUx
            << ",Re="
            << Re;
    }


    fileUx << "\n";



    for(int gp : gridPts)
    {

        int index =
            gp-1;


        double x =
            double(index)/128.0;



        fileUx
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << x;



        for(int Re : Reynolds)
        {

            if(Ux[Re].value.size() > index)
            {
                fileUx
                    << ","
                    << std::scientific
                    << Ux[Re].value[index];
            }
            else
            {
                fileUx << ",";
            }

        }


        fileUx << "\n";

    }



    fileUx.close();



    //--------------------------------------------------
    // Write Uy table
    //--------------------------------------------------

    std::ofstream fileUy("ghia_v.csv");



    fileUy
        << "gridpt,y";


    for(int Re : Reynolds)
    {
        fileUy
            << ",Re="
            << Re;
    }


    fileUy << "\n";



    for(int gp : gridPts)
    {

        int index =
            gp-1;



        double y =
            double(index)/128.0;



        fileUy
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << y;



        for(int Re : Reynolds)
        {

            if(Uy[Re].value.size() > index)
            {
                fileUy
                    << ","
                    << std::scientific
                    << Uy[Re].value[index];
            }
            else
            {
                fileUy << ",";
            }

        }


        fileUy << "\n";

    }



    fileUy.close();



    std::cout
        << "\nPostprocess completed\n"
        << "Generated:\n"
        << "  ghia_u.csv\n"
        << "  ghia_v.csv\n";



    return 0;
}
