#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "Lattice.hpp"

class Boundary
{

public:

```
Boundary(int nx,
         int ny,
         double lidVelocity);


void apply(Lattice& lattice);
```

private:

```
int nx;
int ny;

double lidVelocity;


// walls

void bounceBackBottom(Lattice& lattice);

void bounceBackLeft(Lattice& lattice);

void bounceBackRight(Lattice& lattice);


// moving lid

void movingTop(Lattice& lattice);


// corners

void bottomLeftCorner(Lattice& lattice);

void bottomRightCorner(Lattice& lattice);

void topLeftCorner(Lattice& lattice);

void topRightCorner(Lattice& lattice);
```

};

#endif
