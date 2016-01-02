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

#include "Model.h"

namespace vmmc
{
    Model::Model() {};

    Model::~Model() {};

#ifndef ISOTROPIC
    double Model::energyCallback(unsigned int particle, const double* position, const double* orientation)
#else
    double Model::energyCallback(unsigned int particle, const double* position)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::energyCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }

#ifndef ISOTROPIC
    double Model::pairEnergyCallback(unsigned int particle1, const double* position1,
        const double* orientation1, unsigned int particle2, const double* position2, const double* orientation2)
#else
    double Model::pairEnergyCallback(unsigned int particle1,
        const double* position1, unsigned int particle2, const double* position2)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::pairEnergyCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }

#ifndef ISOTROPIC
    unsigned int Model::interactionsCallback(unsigned int particle,
        const double* position, const double* orientation, unsigned int* interactions)
#else
    unsigned int Model::interactionsCallback(unsigned int particle,
        const double* position, unsigned int* interactions)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::interactionsCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }

#ifndef ISOTROPIC
    void Model::postMoveCallback(unsigned int particle, const double* position, const double* orientation)
#else
    void Model::postMoveCallback(unsigned int particle, const double* position)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::postMoveCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }

#ifndef ISOTROPIC
    double Model::nonPairwiseCallback(unsigned int particle, const double* position, const double* orientation)
#else
    double Model::nonPairwiseCallback(unsigned int particle, const double* position)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::nonPairwiseCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }

#ifndef ISOTROPIC
    bool Model::boundaryCallback(unsigned int particle, const double* position, const double* orientation)
#else
    bool Model::boundaryCallback(unsigned int particle, const double* position)
#endif
    {
        std::cerr << "[ERROR] Model: Virtual function Model::boundaryCallback() must be defined.\n";
        exit(EXIT_FAILURE);
    }
}
