#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
#include <cmath>
#include <limits>


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


        if(ss >> x >> u)
        {
            profile.coord.push_back(x);
            profile.value.push_back(u);
        }

    }


    return profile;
}



//==================================================
// Find nearest grid point
//==================================================

int nearestIndex(const Profile& p,double target)
{

    double minDist =
        std::numeric_limits<double>::max();

    int index = -1;


    for(size_t i=0;i<p.coord.size();i++)
    {

        double dist =
            std::abs(p.coord[i]-target);


        if(dist < minDist)
        {
            minDist = dist;
            index = static_cast<int>(i);
        }

    }


    return index;
}



//==================================================
// Main
//==================================================

int main()
{

    const std::vector<int> Reynolds =
    {
        100,
        400,
        1000,
        3200,
        5000,
        7500,
        10000
    };



    // punti Ghia (129 nodi)
    const std::vector<int> gridPts =
    {
        129,
        126,
        125,
        124,
        123,
        110,
        95,
        80,
        65,
        59,
        31,
        23,
        14,
        10,
        9,
        8,
        1
    };



    std::map<int,Profile> Ux;
    std::map<int,Profile> Uy;



    for(int Re : Reynolds)
    {

        Ux[Re] =
            readProfile(
                "Ux_Re"
                +
                std::to_string(Re)
                +
                ".dat"
            );


        Uy[Re] =
            readProfile(
                "Uy_Re"
                +
                std::to_string(Re)
                +
                ".dat"
            );

    }



//==================================================
// Ux centerline
// Ghia: u(y) con x=0.5
//==================================================


    std::ofstream fileUx("ghia_u.csv");


    fileUx
        << "gridpt,y";


    for(int Re:Reynolds)
        fileUx << ",Re=" << Re;


    fileUx << "\n";



    for(int gp:gridPts)
    {

        double y =
            static_cast<double>(gp-1)/128.0;


        fileUx
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << y;


        for(int Re:Reynolds)
        {

            int id =
                nearestIndex(
                    Ux[Re],
                    y
                );


            if(id>=0)
            {
                fileUx
                    << ","
                    << std::scientific
                    << Ux[Re].value[id];
            }
            else
            {
                fileUx << ",";
            }

        }


        fileUx << "\n";

    }


    fileUx.close();





//==================================================
// Uy centerline
// Ghia: v(x) con y=0.5
//==================================================


    std::ofstream fileUy("ghia_v.csv");


    fileUy
        << "gridpt,x";


    for(int Re:Reynolds)
        fileUy << ",Re=" << Re;


    fileUy << "\n";



    for(int gp:gridPts)
    {

        double x =
            static_cast<double>(gp-1)/128.0;


        fileUy
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << x;



        for(int Re:Reynolds)
        {


            int id =
                nearestIndex(
                    Uy[Re],
                    x
                );


            if(id>=0)
            {
                fileUy
                    << ","
                    << std::scientific
                    << Uy[Re].value[id];
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
        << " ghia_u.csv\n"
        << " ghia_v.csv\n";


    return 0;

}
