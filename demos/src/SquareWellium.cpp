/*
  Copyright (c) 2015-2016 Lester Hedges <lester.hedges+vmmc@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <limits>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "SquareWellium.h"

double INF = std::numeric_limits<double>::infinity();

SquareWellium::SquareWellium(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :

    box(box_),
    particles(particles_),
    cells(cells_),
    maxInteractions(maxInteractions_),
    interactionEnergy(interactionEnergy_),
    interactionRange(interactionRange_),
    squaredCutOffDistance(interactionRange*interactionRange)
{
}

#ifndef ISOTROPIC
double SquareWellium::energyCallback(unsigned int particle, const double* position, const double* orientation)
#else
double SquareWellium::energyCallback(unsigned int particle, const double* position)
#endif
{
    // N.B. This method is somewhat redundant since the same functionality
    // could be achieved by using a combination of the interactionsCallback
    // and model specific pairEnergyCallback methods.

    // Energy counter.
    double energy = 0;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                // Calculate model specific pair energy.
#ifndef ISOTROPIC
                energy += pairEnergyCallback(particle, position, orientation,
                          neighbour, &particles[neighbour].position[0],
                          &particles[neighbour].orientation[0]);
#else
                energy += pairEnergyCallback(particle, position,
                          neighbour, &particles[neighbour].position[0]);
#endif

                // Early exit test for hard core overlaps and large finite energy repulsions.
                if (energy > 1e6) return INF;
            }
        }
    }

    return energy;
}

#ifndef ISOTROPIC
double SquareWellium::pairEnergyCallback(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
#else
double SquareWellium::pairEnergyCallback(unsigned int particle1,
    const double* position1, unsigned int particle2, const double* position2)
#endif
{
    // Separation vector.
    std::vector<double> sep(box.dimension);

    // Calculate separation.
    for (unsigned int i=0;i<box.dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    box.minimumImage(sep);

    double normSqd = 0;

    // Calculate squared norm of vector.
    for (unsigned int i=0;i<box.dimension;i++)
        normSqd += sep[i]*sep[i];

    if (normSqd < 1) return INF;
    if (normSqd < squaredCutOffDistance) return -interactionEnergy;
    return 0;
}

#ifndef ISOTROPIC
unsigned int SquareWellium::interactionsCallback(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
#else
unsigned int SquareWellium::interactionsCallback(unsigned int particle,
    const double* position, unsigned int* interactions)
#endif
{
    // Interaction counter.
    unsigned int nInteractions = 0;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        // Cell index.
        unsigned int cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            // Index of neighbouring particle.
            unsigned int neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = position[k] - particles[neighbour].position[k];

                // Enforce minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Particles interact.
                if (normSqd < squaredCutOffDistance)
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] SquareWellium: Maximum number of interactions exceeded!\n";
                        exit(EXIT_FAILURE);
                    }

                    interactions[nInteractions] = neighbour;
                    nInteractions++;
                }
            }
        }
    }

    return nInteractions;
}

#ifndef ISOTROPIC
void SquareWellium::postMoveCallback(unsigned int particle, const double* position, const double* orientation)
#else
void SquareWellium::postMoveCallback(unsigned int particle, const double* position)
#endif
{
    // Copy coordinates/orientations.
    for (unsigned int i=0;i<box.dimension;i++)
    {
        particles[particle].position[i] = position[i];
#ifndef ISOTROPIC
        particles[particle].orientation[i] = orientation[i];
#endif
    }

    // Calculate the particle's cell index.
    unsigned int newCell = cells.getCell(particles[particle]);

    // Update cell lists if necessary.
    if (particles[particle].cell != newCell)
        cells.updateCell(newCell, particles[particle], particles);
}

#ifndef ISOTROPIC
double SquareWellium::nonPairwiseCallback(unsigned int particle, const double* position, const double* orientation)
#else
double SquareWellium::nonPairwiseCallback(unsigned int particle, const double* position)
#endif
{
    return 0;
}

#ifndef ISOTROPIC
bool SquareWellium::boundaryCallback(unsigned int particle, const double* position, const double* orientation)
#else
bool SquareWellium::boundaryCallback(unsigned int particle, const double* position)
#endif
{
    return false;
}

double SquareWellium::getEnergy()
{
    double energy = 0;

    for (unsigned int i=0;i<particles.size();i++)
#ifndef ISOTROPIC
        energy += energyCallback(i, &particles[i].position[0], &particles[i].orientation[0]);
#else
        energy += energyCallback(i, &particles[i].position[0]);
#endif

    return energy/(2*particles.size());
}
