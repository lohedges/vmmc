/*
  Copyright (c) 2015 Lester Hedges <lester.hedges+vmmc@gmail.com>

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

#include "LennardJonesium.h"

LennardJonesium::LennardJonesium(Box& box_,
                             std::vector <Particle>& particles_,
                             CellList& cells_,
                             unsigned int maxInteractions_,
                             double interactionEnergy_,
                             double interactionRange_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_)
{
    double cutOffDistance = 1.0 + interactionRange;

    // Work out the potential shift.
    potentialShift = std::pow(1.0/cutOffDistance, 12) - std::pow(1/cutOffDistance, 6);
}

double LennardJonesium::computeEnergy(unsigned int particle, double position[], double orientation[])
{
    unsigned int cell;      // cell index
    unsigned int neighbour; // index of neighbouring particle
    double energy = 0;      // energy counter

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particles[particle].cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure the particles are different.
            if (neighbour != particle)
            {
                std::vector <double> sep(box.dimension);

                // Calculate separation.
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
                    double r2Inv = 1.0 / normSqd;
                    double r6Inv = r2Inv*r2Inv*r2Inv;
                    energy += 4.0*interactionEnergy*((r6Inv*r6Inv) - r6Inv - potentialShift);
                }
            }
        }
    }

    return energy;
}

double LennardJonesium::computePairEnergy(unsigned int particle1, double position1[],
    double orientation1[], unsigned int particle2, double position2[], double orientation2[])
{
    // Separation vector.
    std::vector <double> sep(box.dimension);

    // Calculate separation.
    for (unsigned int i=0;i<box.dimension;i++)
        sep[i] = position1[i] - position2[i];

    // Enforce minimum image.
    box.minimumImage(sep);

    double normSqd = 0;

    // Calculate squared norm of vector.
    for (unsigned int i=0;i<box.dimension;i++)
        normSqd += sep[i]*sep[i];

    // Particles interact.
    if (normSqd < squaredCutOffDistance)
    {
        double r2Inv = 1.0 / normSqd;
        double r6Inv = r2Inv*r2Inv*r2Inv;
        return 4.0*interactionEnergy*((r6Inv*r6Inv) - r6Inv - potentialShift);
    }
    else return 0;
}
