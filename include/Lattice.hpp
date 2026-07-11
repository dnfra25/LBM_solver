#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <array>
#include <vector>

class Lattice
{
public:

    // ============================
    // Constructor
    // ============================

    Lattice(int nx,
            int ny,
            double omega_p,
            double omega_m,
            double rho0 = 1.0);

    // ============================
    // Initialization
    // ============================

    void initialize();

    // ============================
    // LBM algorithm
    // ============================

    void computeMacroscopic();

    void computeEquilibrium();

    void collision();

    void streaming();

    // ============================
    // Getters
    // ============================

    int getNx() const { return nx; }

    int getNy() const { return ny; }

    double getRho(int x,int y) const;

    double getUx(int x,int y) const;

    double getUy(int x,int y) const;

    // ============================
    // Access to distributions
    // (needed by Boundary.cpp)
    // ============================

    std::vector<double>& distributions()
    {
        return f;
    }

    std::vector<double>& postCollision()
    {
        return f_post;
    }

    // ============================
    // D2Q9 information
    // ============================

    const std::array<std::array<int,2>,9>& velocities() const
    {
        return c;
    }

    const std::array<int,9>& oppositeDirections() const
    {
        return opposite;
    }

    const std::array<double,9>& weights() const
    {
        return w;
    }

    // ============================
    // Index utilities
    // ============================

    inline int index(int q,int x,int y) const
    {
        return q*nx*ny + x*ny + y;
    }

    inline int cell(int x,int y) const
    {
        return x*ny + y;
    }

private:

    // ============================
    // Grid
    // ============================

    int nx;
    int ny;

    static constexpr int Q = 9;

    // ============================
    // TRT parameters
    // ============================

    double omega_p;
    double omega_m;

    double rho0;

    // ============================
    // D2Q9
    // ============================

    std::array<std::array<int,2>,9> c;

    std::array<double,9> w;

    std::array<int,9> opposite;

    // ============================
    // Distribution functions
    // ============================

    std::vector<double> f;

    std::vector<double> f_eq;

    std::vector<double> f_post;

    // ============================
    // Macroscopic variables
    // ============================

    std::vector<double> rho;

    std::vector<double> ux;

    std::vector<double> uy;

};

#endif
