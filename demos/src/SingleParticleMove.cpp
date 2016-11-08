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

#include "Box.h"
#include "CellList.h"
#include "Model.h"
#include "Particle.h"
#include "SingleParticleMove.h"

SingleParticleMove::SingleParticleMove(
    Model* model_,
    double maxTrialTranslation_,
    double maxTrialRotation_,
    double probTranslate_,
    bool isIsotropic_) :

    model(model_),
    maxTrialTranslation(maxTrialTranslation_),
    maxTrialRotation(maxTrialRotation_),
    probTranslate(probTranslate_),
    isIsotropic(isIsotropic_)
{
    // Check dimensionality.
    if (model->box.dimension == 3) is3D = true;
    else is3D = false;

    // Allocate memory.
    moveParams.trialVector.resize(model->box.dimension);

    // Ignore rotations if potential is isotropic.
    if (isIsotropic) probTranslate = 1.0;
}

void SingleParticleMove::step(const int nSteps)
{
    for (int i=0;i<nSteps;i++)
        step();
}

void SingleParticleMove::operator ++ (const int)
{
    step();
}

void SingleParticleMove::operator += (const int nSteps)
{
    step(nSteps);
}

void SingleParticleMove::step()
{
    // Increment number of attempted moves.
    nAttempts++;

    // Propose a trial move.
    proposeMove();

    // Check whether move was accepted.
    if (accept())
    {
        if (!moveParams.isRotation)
        {
            // Increment number of accepted moves.
            nAccepts++;

            // Increment number of rotations.
            nRotations += moveParams.isRotation;

            unsigned int oldCell = moveParams.preMoveParticle.cell;
            unsigned int newCell = model->particles[moveParams.seed].cell;

            // update cell list
            if (oldCell != newCell)
            {
                model->particles[moveParams.seed].cell = oldCell;
                model->cells.updateCell(newCell, model->particles[moveParams.seed], model->particles);
            }
        }
    }
    else
    {
        // Revert particle to pre-move state.
        model->particles[moveParams.seed] = moveParams.preMoveParticle;
    }
}

unsigned long long SingleParticleMove::getAttempts() const
{
    return nAttempts;
}

unsigned long long SingleParticleMove::getAccepts() const
{
    return nAccepts;
}

unsigned long long SingleParticleMove::getRotations() const
{
    return nRotations;
}

void SingleParticleMove::reset()
{
    nAttempts = nAccepts = nRotations = 0;
}

void SingleParticleMove::proposeMove()
{
    // Choose a seed particle.
    moveParams.seed = rng.integer(0, model->particles.size()-1);

    // Choose a random point on the surface of the unit sphere/circle.
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] = rng.normal();

    // Normalise the trial vector.
    double norm = computeNorm(moveParams.trialVector);
    for (unsigned int i=0;i<model->box.dimension;i++)
        moveParams.trialVector[i] /= norm;

    // Choose the move type.
    if (rng() < probTranslate)
    {
        // Translation.
        moveParams.isRotation = false;

        // Scale step-size to uniformly sample unit sphere/circle.
        if (is3D) moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/3.0);
        else moveParams.stepSize = maxTrialTranslation*std::pow(rng(), 1.0/2.0);
    }
    else
    {
        // Rotation.
        moveParams.isRotation = true;
        moveParams.stepSize = maxTrialRotation*(2.0*rng()-1.0);
    }

    // Calculate pre-move energy.
#ifndef ISOTROPIC
    double initialEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0]);
#else
    double initialEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0]);
#endif

    // Store initial coordinates/orientation.
    moveParams.preMoveParticle = model->particles[moveParams.seed];

    // Execute the move.
    if (!moveParams.isRotation) // Translation.
    {
        for (unsigned int i=0;i<model->box.dimension;i++)
            model->particles[moveParams.seed].position[i] += moveParams.stepSize*moveParams.trialVector[i];

        // Apply periodic boundary conditions.
        model->box.periodicBoundaries(model->particles[moveParams.seed].position);

        // Work out new cell index.
        model->particles[moveParams.seed].cell = model->cells.getCell(model->particles[moveParams.seed]);
    }
    else                        // Rotation.
    {
        std::vector<double> vec(model->box.dimension);

        // Calculate orientation rotation vector.
        if (is3D) rotate3D(model->particles[moveParams.seed].orientation, moveParams.trialVector, vec, moveParams.stepSize);
        else rotate2D(model->particles[moveParams.seed].orientation, vec, moveParams.stepSize);

        // Update orientation.
        for (unsigned int i=0;i<model->box.dimension;i++)
            model->particles[moveParams.seed].orientation[i] += vec[i];
    }

    // Calculate post-move energy.
#ifndef ISOTROPIC
    double finalEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0],
        &model->particles[moveParams.seed].orientation[0]);
#else
    double finalEnergy = model->computeEnergy(moveParams.seed,
        &model->particles[moveParams.seed].position[0]);
#endif

    energyChange = finalEnergy - initialEnergy;
}

bool SingleParticleMove::accept()
{
    if (energyChange == 0) return true;
    if (energyChange == INF) return false;
    if (rng() < exp(-energyChange)) return true;
    else return false;
}

void SingleParticleMove::rotate3D(std::vector<double>& v1, std::vector<double>& v2, std::vector<double>& v3, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    double v1Dotv2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

    v3[0] = ((v1[0] - v2[0]*v1Dotv2))*(c - 1) + (v2[2]*v1[1] - v2[1]*v1[2])*s;
    v3[1] = ((v1[1] - v2[1]*v1Dotv2))*(c - 1) + (v2[0]*v1[2] - v2[2]*v1[0])*s;
    v3[2] = ((v1[2] - v2[2]*v1Dotv2))*(c - 1) + (v2[1]*v1[0] - v2[0]*v1[1])*s;
}

void SingleParticleMove::rotate2D(std::vector<double>& v1, std::vector<double>& v2, double angle)
{
    double c = cos(angle);
    double s = sin(angle);

    v2[0] = (v1[0]*c - v1[1]*s) - v1[0];
    v2[1] = (v1[0]*s + v1[1]*c) - v1[1];
}

double SingleParticleMove::computeNorm(std::vector<double>& vec)
{
    double normSquared = 0;

    for (unsigned int i=0;i<vec.size();i++)
        normSquared += vec[i]*vec[i];

    return sqrt(normSquared);
}
