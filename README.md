# LibVMMC

Copyright &copy; 2015 Lester Hedges.

Released under the [GPL](http://www.gnu.org/copyleft/gpl.html).

## About
A simple C++ library to implement the "virtual-move" Monte Carlo (VMMC)
algorithm of [Steve Whitelam](http://nanotheory.lbl.gov/people/SteveWhitelam.html)
and [Phill Geissler](http://www.cchem.berkeley.edu/plggrp/index.html), see:

* Avoiding unphysical kinetic traps in Monte Carlo simulations of strongly
attractive particles, S. Whitelam and P.L. Geissler,
[Journal of Chemical Physics, 127, 154101 (2007)](http://dx.doi.org/10.1063/1.2790421)

* Approximating the dynamical evolution of systems of strongly interacting
overdamped particles, S. Whitelam,
[Molecular Simulation, 37 (7) (2011)](http://dx.doi.org/10.1080/08927022.2011.565758).
(Preprint version available [here](http://arxiv.org/abs/1009.2008).)

Our primary goal is to make VMMC accessible to a wider audience, for whom the
time required to code the algorithm poses a significant barrier to using the
method. This allows the user to focus on model development.

## Installation
A `Makefile` is included for building and installing LibVMMC.

To compile LibVMMC, then install the library, documentation, and demos:

```
$ make build
$ make install
```

By default, the library installs to `/usr/local`. Therefore, you may need admin
priveleges for the final `make install` step above. An alternative is to change
the install location:

```bash
$ PREFIX=MY_INSTALL_DIR make install
```

Further details on using the Makefile can be found by running make without
a target, i.e.

```bash
$ make
```

## Compiling and linking
To use LibVMMC with a C/C++ code first include the LibVMMC header file somewhere
in the code.

```
//example.cpp
#include <vmmc/VMMC.h>
```

Then to compile, we can use something like the following:

```bash
$ g++ example.cpp -lvmmc
```

This assumes that we have used the default install location /usr/local/. If
we specify an install location, we would use a command more like the following:

```bash
$ g++ example.cpp -I/my/path/include -L/my/path/lib -lvmmc
```

## Dependencies
LibVMMC uses the [Mersenne Twister](http://en.wikipedia.org/wiki/Mersenne_Twister)
psuedorandom number generator. A C++ implementation is included as a bundled
header file and is compiled into the libvmmc library. The included
`MersenneTwister.h` is taken from
[here](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/c-lang.html).

## Callback functions
LibVMMC works via four user-defined callback functions that abstract model specific
details, such as the pair potential. These callbacks have the following
prototypes:

### Particle energy
Calculate the total pair interaction energy felt by a particle.
```
typedef double (*energyCallback) (unsigned int index, double position[], double orientation[]);
```
`index` = The particle index.

`position` = An x, y, z (or x, y in 2D) coordinate vector for the particle.

`orientation` = The particle's orientation unit vector.

### Pair energy
Calculate the pair interaction between two particles.
```
typedef double (*pairEnergyCallback) (unsigned int index1, double position1[], double orientation1[], unsigned int index2, double position2[], double orientation2[]);
```
`index1` = The index of the first particle.

`position1` = The coordinate vector of the first particle.

`orientation1` = The orientation unit vector of the first particle.

`index2` = The index of the second particle.

`position2` = The coordinate vector of the second particle.

`orientation2` = The orientation unit vector of the second particle.

### Interactions
Determine the interactions for a given particle.
```
typedef unsigned int (*interactionsCallback) (unsigned int index, double position[], double orientation[], unsigned int interactions[]);
```
`index` = The index of the  particle.

`position` = The coordinate vector of the particle.

`orientation` = The orientation unit vector of the particle.

`interactions` = An array to store the indices of the interactions.

### Post-move
Apply any post-move updates, e.g. update cell lists, or neighbour lists.
```
typedef void (*postMoveCallback) (unsigned int index, double position[], double orientation[]);
```
`index` = The index of the  particle.

`position` = The coordinate vector of the particle following the move.

`orientation` = The orientation unit vector of the particle following the move.

## The VMMC object
To use LibVMMC you will want to create an instance of the VMMC object. This has the following
constructor:
```
VMMC(unsigned int nParticles, unsigned int dimension, double coordinates[], double orientations[], double maxTrialTranslation, double maxTrialRotation, double probTranslate, double referenceRadius, unsigned int maxInteractions, double boxSize[], bool isRepulsive, energyCallback computeEnergy, pairEnergyCallback computePairEnergy, interactionsCallback computeInteractions, postMoveCallback, applyPostMoveUpdates);
```
`nParticles` = The number of particles in the simulation box.

`dimension` = The dimension of the simulation box (either 2 or 3).

`coordinates` = An array containing coordinates for all of the particles in the
system, i.e. `x1, y1, z1, x2, y2, z2, ... , xN, yN, zN.`
Coordinates should run from 0 to the box size in each dimension.

`orientations` = An array containing orientations (unit vectors) for all of the
particles in the system, i.e. `nx1, ny1, nz1, nx2, ny2, nz2, ... , nxN, nyN, nzN.`

`maxTrialTranslation` = The maximum trial translation, in units of the particle
diameter (or typical particle size).

`maxTrialTranslation` = The maximum trial rotation in radians.

`probTranslate` = The probability of attempting a translation move (relative to rotations).

`referenceRadius` = A reference radius for computing the approximate hydrodynamic
damping factor, e.g. the radius of a typical particle in the system.

`maxInteractions` = The maximum number of pair interactions that an individual
particle can make. This will be used to resize LibVMMC's internal data
structures and the user should assert that this limit isn't exceed in the
`interactionsCallback` function. The number can be chosen from the symmetry
of the system, e.g. if particles can only make a certain number of patchy
interactions, or by estimating the average number of neigbours within the
interaction volume around a particle.

`boxSize` = The base length of the simulation box in each dimension.

`isRepulsive` = Whether the potential has finite energy repulsions.

`computeEnergy` = The callback function to calculate the total pair interaction
for a particle.

`computePairEnergy` = The callback function to calculate the pair interaction
between two particles.

`computeInteractions` = The callback function to determine the neighbours with
which a particle interacts.

`applyPostMoveUpdates` = The callback function to perform any required updates
following the move. Here you should copy the updated particle position and
orientation back into your own data structures and implement any additional
updates, e.g. cell lists.

## Executing a virtual-move
Once an instance of the VMMC object is created, e.g.
```
VMMC(...) vmmc;
```
then a single trial move can be executed as follows:
```
vmmc.step();
```
To perform 1000 trial moves:
```
vmmc.step(1000);
```
The same can be achieved by using the overloaded `++` and `+=` operators,
i.e. `vmmc++` for a single step, and `vmmc += 1000` for 1000 steps.

## Demos
The following example codes showing how to interface with LibVMMC are included
in the `demos` directory.

* `square_wellium.cpp`: A simulation of a square-well fluid in two- or three-dimensions.
* `lennard_jones.cpp`: A simulation of a Lennard-Jones fluid in two- or three-dimensions.

The demo code also illustrates how to implement efficient, dynamically
updated cell lists. See `demos/src/CellList.h` and `demos/src/CellList.cpp`
for implementation details.

## Limitations
* The use of simple C-style callback functions means that the user will likely need
to use globals for several variables and data structures. Since this library is
intended to be used for simulation of relatively simple models, likely with a
small code base, this was deemed as a reasonable trade-off in terms of preserving
generality. The [example](#demos) code illustrates some simple examples.
* For spherical particles bearing isotropic interactions, e.g. the square-well
fluid, single particle rotations will always be accepted. While not a problem
from a thermodynamic perspective, this may cause issues if the user wishes to
enforce a strict Stoke's scaling of translational and rotational diffusion.
* The recursive manner in which the trial cluster is built can lead to a stack
overflow if the cluster contains many particles. Typically, thousands, or tens
of thousands of particles should be perfectly manageable. The typically memory
footprint for a simulation of 1000 particles is around 2.5MB for hard particle
simulations. This is roughly doubled if the potential has finite energy repulsions.

## Efficiency
In aid of generality there are several sources of redundancy that impact the
efficiency of the VMMC implementation. As written, LibVMMC performs around 3-4
times worse than a fully optimised VMMC code for square-well fluids. A few
efficiency considerations are listed below in case the user wishes to modify
the VMMC source code in order to improve performance.

* When calculating a list of neighbours with which a given particle interacts
it's likely that you'll need to calculate the pair interaction energy. For
certain models it may be more efficient to return a list of pair energies
along with the interactions, rather than having to recalculate them.
* For models with an isotropic interaction of fixed energy scale the pair
energy is simple a constant. As such, the pair energy calculation is entirely
redundant, i.e. if we knowing that two particles interact is enough to know
the pair energy.
* If using cell lists, the typical size of a trial displacement will be small
enough such that a particle stays within the same neigbourhood of cells
following the trial move. As such, there is often no need to update cell lists
until confirming that the post-move configuration is valid, e.g. no overlaps.
At present the same `postMoveCallback` function is called twice: once in order
to apply the move; again if the move is subsequently rejected. This means that
the cell lists will be updated twice if a move is rejected.

## Tips
* LibVMMC currently assumes that the simulation box is periodic in all dimensions.
To impose non-periodic boundaries simply check whether the move leads to a particle
being displaced by more than half the box width along the restricted dimension and
return an appropriately large energy so that the move will be rejected.
* It is not a requirement that all particles in the simulation box be of the same
type. Make use of the particle indices that are passed to callback functions in
order to distinguish different species.
