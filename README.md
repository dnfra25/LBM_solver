# A Fluid Dynamic Solver Based on the Lattice Boltzmann Method

**Author:** Francesco Di Napoli

## Overview

This project implements a two-dimensional fluid dynamic solver based on the **Lattice Boltzmann Method (LBM)**.

The solver is designed to simulate incompressible fluid flows using the **lid-driven cavity benchmark problem**, a classical test case widely used to validate numerical methods for computational fluid dynamics.

The implementation is based on a **D2Q9 lattice model** and uses a **Two-Relaxation-Time (TRT) collision operator** to improve numerical stability and accuracy compared to the standard single-relaxation-time BGK formulation.

The results are generated for different Reynolds numbers and can be compared with the reference benchmark data reported by U. Ghia, K. N. Ghia, and C. T. Shin in their work on the lid-driven cavity flow.

---

## Physical Problem

The simulated problem is a two-dimensional square cavity:

* The bottom wall and the lateral walls are stationary.
* The upper wall (lid) moves with a constant velocity.
* The motion of the lid generates a recirculating flow inside the cavity.

The simulations are performed for the following Reynolds numbers:

```
Re = 100
Re = 400
Re = 1000
Re = 3200
Re = 5000
Re = 7500
Re = 10000
```

For every case, the solver computes the velocity field and extracts the horizontal and vertical centerline velocity profiles.

---

## Numerical Method

### Lattice Boltzmann Model

The solver uses a D2Q9 lattice:

* 2 spatial dimensions
* 9 discrete velocity directions

The particle distribution functions are evolved according to the standard LBM procedure:

1. Computation of macroscopic quantities:

   * density
   * x velocity component
   * y velocity component

2. Equilibrium distribution evaluation

3. Collision step

4. Streaming step

5. Application of physical boundary conditions

The LBM time evolution is therefore:

```
f → macroscopic variables → equilibrium → collision → streaming → boundary treatment
```

---

## Collision Model

The collision operator is based on the **Two-Relaxation-Time (TRT)** formulation.

The distribution functions are decomposed into:

* symmetric component
* antisymmetric component

Each component is relaxed using a different relaxation parameter:

* `omega_plus`
* `omega_minus`

The TRT model provides improved stability, especially for simulations at higher Reynolds numbers.

---

## Boundary Conditions

The cavity walls are implemented through specific boundary treatments:

* **Bottom wall:** bounce-back condition
* **Left wall:** bounce-back condition
* **Right wall:** bounce-back condition
* **Top moving lid:** moving-wall condition with prescribed velocity

The corner nodes are explicitly treated to correctly handle the intersection between different boundaries.

---

## Project Structure

```
.
├── main.cpp
├── cavity.cpp
├── cavity.hpp
├── lattice.cpp
├── lattice.hpp
├── boundary.cpp
├── boundary.hpp
└── postprocess.cpp
```

### Main solver components

### `Lattice`

Handles the LBM numerical model:

* D2Q9 lattice definition
* equilibrium distribution
* macroscopic variables
* TRT collision
* streaming operation

---

### `Cavity`

Represents the physical problem:

* cavity geometry
* Reynolds number definition
* viscosity calculation
* relaxation parameters computation
* convergence monitoring
* velocity profile extraction

---

### `Boundary`

Manages all physical boundary conditions:

* stationary walls
* moving lid
* corner treatment

---

### `postprocess`

Processes the raw solver output.

Its purpose is not to solve the fluid problem but to format the computed profiles into CSV tables with a structure suitable for comparison with the benchmark data of:

* U. Ghia
* K. N. Ghia
* C. T. Shin

---

# Compilation and Execution

The project produces two executables:

```
cavity
postprocess
```

The execution order is important.

## 1. Run the LBM solver

First execute:

```
./cavity
```

This runs all Reynolds number cases and generates the velocity profiles:

```
Ux_Re100.dat
Uy_Re100.dat

Ux_Re400.dat
Uy_Re400.dat

...

Ux_Re10000.dat
Uy_Re10000.dat
```

The output files contain:

* horizontal centerline velocity (`Ux`)
* vertical centerline velocity (`Uy`)

---

## 2. Run the post-processing tool

After the solver has completed, execute:

```
./postprocess
```

This reads the generated profiles and creates:

```
ghia_u.csv
ghia_v.csv
```

These files are formatted to reproduce the structure commonly used for comparison with the Ghia benchmark data.

---

# Convergence Criterion

The solver monitors the relative velocity variation between consecutive iterations.

The simulation stops when:

```
error < 1.0e-8
```

or when the maximum number of iterations is reached:

```
maxIterations = 20000
```

---

# Parallelization

The implementation supports OpenMP parallelization.

The main computational loops are parallelized using OpenMP directives, including:

* initialization
* macroscopic variable computation
* equilibrium calculation
* collision
* streaming
* convergence evaluation

If compiled with OpenMP support, the number of available threads is printed at runtime.

---

# Output

During execution, the solver prints:

* Reynolds number
* lid velocity
* viscosity
* relaxation parameters
* convergence history
* final velocity diagnostics
* execution time

Example generated profiles:

```
Ux_Re100.dat
Uy_Re100.dat
```

contain normalized coordinates and velocity values suitable for visualization and validation.

---

# Validation

The solver is intended to be validated through comparison with the classical lid-driven cavity benchmark:

**U. Ghia, K. N. Ghia, and C. T. Shin**
*High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method*

The generated CSV files allow direct comparison between the LBM results and the reference centerline velocity profiles.

---

# Author

**Francesco Di Napoli**

Project:
**A Fluid Dynamic Solver Based on the Lattice Boltzmann Method**
