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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "PatchyDisc.h"

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif

PatchyDisc::PatchyDisc(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_)
{
#ifdef ISOTROPIC
    std::cerr << "[ERROR] PatchyDisc: Cannot be used with isotropic VMMC library!\n";
    exit(EXIT_FAILURE);
#endif

    // Check dimensionality.
    if (box.dimension != 2)
    {
        std::cerr << "[ERROR] PatchyDisc: Model only valid in two dimensions!\n";
        exit(EXIT_FAILURE);
    }

    // Work out the angle between patches.
    patchSeparation = 2.0*M_PI/maxInteractions;

    // Resize rotation matrix lookup tables.
    cosTheta.resize(maxInteractions);
    sinTheta.resize(maxInteractions);

    // Populate lookup tables.
    for (unsigned int i=0;i<maxInteractions;i++)
    {
        cosTheta[i] = cos(i*patchSeparation);
        sinTheta[i] = sin(i*patchSeparation);
    }
}

double PatchyDisc::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
{
    // Separation vector.
    std::vector<double> sep(2);

    // Calculate disc separation.
    sep[0] = position1[0] - position2[0];
    sep[1] = position1[1] - position2[1];

    // Enforce minimum image.
    box.minimumImage(sep);

    // Calculate squared norm of vector.
    double normSqd = sep[0]*sep[0] + sep[1]*sep[1];

    // Discs overlap.
    if (normSqd < 1) return INF;

    // Total interaction energy sum.
    double energy = 0;

    // Test interactions between all patch pairs.
    for (unsigned int i=0;i<maxInteractions;i++)
    {
        // Compute position of patch i on first disc.
        std::vector<double> coord1(2);
        coord1[0] = position1[0] + 0.5*(orientation1[0]*cosTheta[i] - orientation1[1]*sinTheta[i]);
        coord1[1] = position1[1] + 0.5*(orientation1[0]*sinTheta[i] + orientation1[1]*cosTheta[i]);

        // Enforce periodic boundaries.
        box.periodicBoundaries(coord1);

        for (unsigned int j=0;j<maxInteractions;j++)
        {
            // Compute position of patch j on second disc.
            std::vector<double> coord2(2);
            coord2[0] = position2[0] + 0.5*(orientation2[0]*cosTheta[j] - orientation2[1]*sinTheta[j]);
            coord2[1] = position2[1] + 0.5*(orientation2[0]*sinTheta[j] + orientation2[1]*cosTheta[j]);

            // Enforce periodic boundaries.
            box.periodicBoundaries(coord2);

            // Calculate patch separation.
            sep[0] = coord1[0] - coord2[0];
            sep[1] = coord1[1] - coord2[1];

            // Enforce minimum image.
            box.minimumImage(sep);

            // Calculate squared norm of vector.
            normSqd = sep[0]*sep[0] + sep[1]*sep[1];

            // Patches interact.
            if (normSqd < squaredCutOffDistance)
                energy -= interactionEnergy;
        }
    }

    return energy;
}

unsigned int PatchyDisc::computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
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
                // Calculate pair energy.
                double energy = computePairEnergy(particle, position, orientation,
                                neighbour, &particles[neighbour].position[0],
                                &particles[neighbour].orientation[0]);

                // Particles interact.
                if (energy < 0)
                {
                    if (nInteractions == maxInteractions)
                    {
                        std::cerr << "[ERROR] PatchyDisc: Maximum number of interactions exceeded!\n";
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
