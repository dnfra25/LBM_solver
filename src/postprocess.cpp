#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>

struct Node
{
    double rho;
    double ux;
    double uy;
};

int main()
{
    //--------------------------------------------------
    // Reynolds numbers available
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
    // Grid points used by Ghia
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

    const int Nx = 129;
    const int Ny = 129;

    const int centerX = Nx/2;
    const int centerY = Ny/2;

    //--------------------------------------------------
    // Storage
    //--------------------------------------------------

    std::map<int,std::vector<Node>> solutions;

    //--------------------------------------------------
    // Read every field
    //--------------------------------------------------

    for(int Re : Reynolds)
    {
        std::string filename =
            "field_Re" + std::to_string(Re) + ".dat";

        std::ifstream file(filename);

        if(!file)
        {
            std::cout
                << "Cannot open "
                << filename
                << std::endl;

            continue;
        }

        std::vector<Node> field(Nx*Ny);

        std::string line;

        while(std::getline(file,line))
        {
            if(line.empty()) continue;

            if(line[0]=='#') continue;

            std::stringstream ss(line);

            int x,y;

            Node n;

            ss >> x
               >> y
               >> n.rho
               >> n.ux
               >> n.uy;

            field[x*Ny+y]=n;
        }

        solutions[Re]=field;

        std::cout
            << "Loaded "
            << filename
            << std::endl;
    }

    //--------------------------------------------------
    // Horizontal centreline (u)
    //--------------------------------------------------

    std::ofstream uFile("ghia_u.csv");

    uFile
        << "gridpt,x";

    for(int Re : Reynolds)
        uFile << ",Re=" << Re;

    uFile << "\n";

    for(int gp : gridPts)
    {
        int x = gp-1;

        double xcoord =
            double(x)/(Nx-1);

        uFile
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << xcoord;

        for(int Re : Reynolds)
        {
            auto it = solutions.find(Re);

            if(it==solutions.end())
            {
                uFile << ",";
                continue;
            }

            const auto& field = it->second;

            double u =
                field[x*Ny+centerY].ux;

            uFile
                << ","
                << std::scientific
                << u;
        }

        uFile << "\n";
    }

    uFile.close();

    //--------------------------------------------------
    // Vertical centreline (v)
    //--------------------------------------------------

    std::ofstream vFile("ghia_v.csv");

    vFile
        << "gridpt,y";

    for(int Re : Reynolds)
        vFile << ",Re=" << Re;

    vFile << "\n";

    for(int gp : gridPts)
    {
        int y = gp-1;

        double ycoord =
            double(y)/(Ny-1);

        vFile
            << gp
            << ","
            << std::fixed
            << std::setprecision(6)
            << ycoord;

        for(int Re : Reynolds)
        {
            auto it = solutions.find(Re);

            if(it==solutions.end())
            {
                vFile << ",";
                continue;
            }

            const auto& field = it->second;

            double v =
                field[centerX*Ny+y].uy;

            vFile
                << ","
                << std::scientific
                << v;
        }

        vFile << "\n";
    }

    vFile.close();

    std::cout
        << "\nGenerated:\n"
        << "   ghia_u.csv\n"
        << "   ghia_v.csv\n";

    return 0;
}
