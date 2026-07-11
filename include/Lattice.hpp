#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <array>
#include <vector>

class Boundary;

class Lattice
{
    friend class Boundary;

public:

    static constexpr int Q = 9;

    //=========================================
    // Constructor
    //=========================================

    Lattice(int nx,
            int ny,
            double omega_p,
            double omega_m,
            double rho0 = 1.0);

    //=========================================
    // Initialization
    //=========================================

    void initialize();

    //=========================================
    // LBM steps
    //=========================================

    void computeMacroscopic();

    void computeEquilibrium();

    void collision();

    void streaming();

    void step();

    //=========================================
    // Getters
    //=========================================

    int getNx() const
    {
        return nx;
    }

    int getNy() const
    {
        return ny;
    }

    double getRho(int x,int y) const;

    double getUx(int x,int y) const;

    double getUy(int x,int y) const;

    //=========================================
    // Distribution access
    //=========================================

    std::vector<double>& distributions()
    {
        return f;
    }

    const std::vector<double>& distributions() const
    {
        return f;
    }

    //=========================================
    // Indexing
    //=========================================

    inline int index(int q,int x,int y) const
    {
        return q*nx*ny + x*ny + y;
    }

    inline int cell(int x,int y) const
    {
        return x*ny + y;
    }

private:

    //-----------------------------------------
    // Grid
    //-----------------------------------------

    int nx;
    int ny;

    //-----------------------------------------
    // TRT parameters
    //-----------------------------------------

    double omega_p;
    double omega_m;

    double rho0;

    //-----------------------------------------
    // D2Q9
    //-----------------------------------------

    std::array<std::array<int,2>,Q> c;

    std::array<double,Q> w;

    std::array<int,Q> opposite;

    //-----------------------------------------
    // Distribution functions
    //-----------------------------------------

    std::vector<double> f;

    std::vector<double> f_eq;

    std::vector<double> f_post;

    //-----------------------------------------
    // Macroscopic fields
    //-----------------------------------------

    std::vector<double> rho;

    std::vector<double> ux;

    std::vector<double> uy;
};

#endif
