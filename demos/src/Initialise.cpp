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

#include <iostream>

#include "Box.h"
#include "CellList.h"
#include "Particle.h"
#include "Initialise.h"
#include "MersenneTwister.h"

Initialise::Initialise()
{
}

void Initialise::random(std::vector<Particle>& particles, CellList& cells, Box& box, MersenneTwister& rng, bool isSpherocylinder)
{
    if (isSpherocylinder && (box.dimension != 3))
    {
        std::cerr << "[ERROR] Initialise: Spherocylindrical boundary only valid for three dimensional simulation box!\n";
        exit(EXIT_FAILURE);
    }

    // Copy box dimensions.
    boxSize = box.boxSize;

    for (unsigned i=0;i<particles.size();i++)
    {
        // Current number of attempted particle insertions.
        unsigned int nTrials = 0;

        // Whether particle overlaps.
        bool isOverlap = true;

        // Temporary vector.
        std::vector<double> vec(box.dimension);

        // Set particle index.
        particles[i].index = i;

        // Keep trying to insert particle until there is no overlap.
        while (isOverlap)
        {
            nTrials++;

            // Generate a random position.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng()*box.boxSize[j];

            particles[i].position = vec;

            // Generate a random orientation.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] = rng.normal();

            // Calculate vector norm.
            double norm = 0;
            for (unsigned int j=0;j<box.dimension;j++)
                norm += vec[j]*vec[j];
            norm = sqrt(norm);

            // Convert orientation to a unit vector.
            for (unsigned int j=0;j<box.dimension;j++)
                vec[j] /= norm;

            particles[i].orientation = vec;

            // Calculate the particle's cell index.
            particles[i].cell = cells.getCell(particles[i]);

            // Enforce spherocylindrical boundary.
            if (isSpherocylinder)
            {
                // Make sure particle lies within the spherocylinder.
#ifndef ISOTROPIC
                if (!outsideSpherocylinder(i, &particles[i].position[0], &particles[i].orientation[0]))
#else
                if (!outsideSpherocylinder(i, &particles[i].position[0]))
#endif
                {
                    // See if there is any overlap between particles.
                    isOverlap = checkOverlap(particles[i], particles, cells, box);
                }
                else isOverlap = true;
            }
            else
            {
                // See if there is any overlap between particles.
                isOverlap = checkOverlap(particles[i], particles, cells, box);
            }

            // Check trial limit isn't exceeded.
            if (nTrials == MAX_TRIALS)
            {
                std::cerr << "[ERROR] Initialise: Maximum number of trial insertions reached.\n";
                exit(EXIT_FAILURE);
            }
        }

        // Update cell list.
        cells.initCell(particles[i].cell, particles[i]);
    }
}

#ifndef ISOTROPIC
bool Initialise::outsideSpherocylinder(unsigned int particle, const double* position, const double* orientation)
#else
bool Initialise::outsideSpherocylinder(unsigned int particle, const double* position)
#endif
{
    // Centre of sphere or circle.
    std::vector<double> centre(3);

    // Separation vector.
    std::vector<double> sep(3);

    // Squared radius of spherical cap (minus squared radius of particle).
    double radiusSqd = 0.25*(boxSize[0] - 1)*(boxSize[0] - 1);

    // Squared norm of separation vector.
    double normSqd = 0;

    // Initialise x and y coordinates of sphere centres.
    centre[0] = 0.5*boxSize[0];
    centre[1] = 0.5*boxSize[0];

    // Check whether particle lies in lower cap.
    if (position[2] < 0.5*boxSize[0])
    {
        centre[2] = 0.5*boxSize[0];

        // Calculate separation.
        for (unsigned int i=0;i<3;i++)
        {
            sep[i] = position[i] - centre[i];
            normSqd += sep[i]*sep[i];
        }

        // Particle lies outside of cap.
        if (normSqd > radiusSqd) return true;
    }
    else
    {
        // Check whether particle lies in upper cap.
        if (position[2] > (boxSize[2] - 0.5*boxSize[0]))
        {
            centre[2] = boxSize[2] - 0.5*boxSize[0];

            // Calculate separation.
            for (unsigned int i=0;i<3;i++)
            {
                sep[i] = position[i] - centre[i];
                normSqd += sep[i]*sep[i];
            }

            // Particle lies outside of cap.
            if (normSqd > radiusSqd) return true;
        }
        else
        {
            // Calculate separation.
            for (unsigned int i=0;i<2;i++)
            {
                sep[i] = position[i] - centre[i];
                normSqd += sep[i]*sep[i];
            }

            // Particle lies outside of cylinder.
            if (normSqd > radiusSqd) return true;
        }
    }

    // Inside spherocylinder.
    return false;
}

bool Initialise::checkOverlap(Particle& particle, std::vector<Particle>& particles, CellList& cells, Box& box)
{
    unsigned int cell, neighbour;

    // Check all neighbouring cells including same cell.
    for (unsigned int i=0;i<cells.getNeighbours();i++)
    {
        cell = cells[particle.cell].neighbours[i];

        // Check all particles within cell.
        for (unsigned int j=0;j<cells[cell].tally;j++)
        {
            neighbour = cells[cell].particles[j];

            // Make sure particles are different.
            if (neighbour != particle.index)
            {
                // Particle separtion vector.
                std::vector<double> sep(box.dimension);

                // Compute separation.
                for (unsigned int k=0;k<box.dimension;k++)
                    sep[k] = particle.position[k] - particles[neighbour].position[k];

                // Compute minimum image.
                box.minimumImage(sep);

                double normSqd = 0;

                // Calculate squared norm of vector.
                for (unsigned int k=0;k<box.dimension;k++)
                    normSqd += sep[k]*sep[k];

                // Overlap if normSqd is less than particle diameter (box is scaled in diameter units).
                if (normSqd < 1) return true;
            }
        }
    }

    // If we get this far, no overlaps.
    return false;
}
