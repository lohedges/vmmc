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

#include "SquareWellium.h"

SquareWellium::SquareWellium(
    Box& box_,
    std::vector<Particle>& particles_,
    CellList& cells_,
    unsigned int maxInteractions_,
    double interactionEnergy_,
    double interactionRange_) :
    Model(box_, particles_, cells_, maxInteractions_, interactionEnergy_, interactionRange_)
{
}

#ifndef ISOTROPIC
double SquareWellium::computePairEnergy(unsigned int particle1, const double* position1,
    const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
#else
double SquareWellium::computePairEnergy(unsigned int particle1,
    const double position1, unsigned int particle2, const double* position2)
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
