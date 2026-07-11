#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <vector>

class Lattice
{
private:

    //-------------------------
    // Grid dimensions
    //-------------------------

    int nx;
    int ny;

    //-------------------------
    // D2Q9 constants
    //-------------------------

    static constexpr int Q = 9;

    int c[Q][2];

    double w[Q];

    int opposite[Q];

    //-------------------------
    // TRT parameters
    //-------------------------

    double omega_p;
    double omega_m;

    //-------------------------
    // Distribution functions
    //-------------------------

    std::vector<double> f;

    std::vector<double> f_eq;

    std::vector<double> f_post;

    //-------------------------
    // Macroscopic variables
    //-------------------------

    std::vector<double> rho;

    std::vector<double> ux;

    std::vector<double> uy;

    //-------------------------
    // Indexing
    //-------------------------

    inline int index(int q,int x,int y) const
    {
        return q*nx*ny + x*ny + y;
    }

    inline int cell(int x,int y) const
    {
        return x*ny + y;
    }

public:

    //-------------------------
    // Constructor
    //-------------------------

    Lattice(int nx,
            int ny,
            double omega_p,
            double omega_m);

    //-------------------------
    // Initialization
    //-------------------------

    void initialize();

    //-------------------------
    // LBM steps
    //-------------------------

    void computeMacroscopic();

    void computeEquilibrium();

    void collision();

    void streaming();

    //-------------------------
    // Accessors
    //-------------------------

    int getNx() const;

    int getNy() const;

    double getRho(int x,int y) const;

    double getUx(int x,int y) const;

    double getUy(int x,int y) const;

    //-------------------------
    // Access to distributions
    //-------------------------

    double& F(int q,int x,int y);

    const double& F(int q,int x,int y) const;
};

#endif
