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

#ifdef ISOTROPIC
#error square_wellium.cpp cannot be linked to isotropic VMMC library!
#endif

#include <Python.h>
#include <cstdlib>
#include <iostream>

#include "VMMC.h"

// Global Python objects.
PyObject *box, *particles, *cells, *io, *squareWellium;

// Globals.
unsigned int dimension;
unsigned int maxInteractions;

// Callback function prototypes.
double computeEnergy(unsigned int, const double*, const double*);
double computePairEnergy(unsigned int, const double*, const double*, unsigned int, const double*, const double*);
unsigned int computeInteractions(unsigned int, const double*, const double*, unsigned int*);
void applyPostMoveUpdates(unsigned int, const double*, const double*);

// Function prototypes.
double getEnergy();
void appendXyzTrajectory(bool);

int main(int argc, char** argv)
{
    // Simulation parameters.
    dimension = 3;                          // dimension of simulation box
    unsigned int nParticles = 1000;         // number of particles
    double interactionEnergy = 2.6;         // pair interaction energy scale (in units of kBT)
    double interactionRange = 1.1;          // size of interaction range (in units of particle diameter)
    double density = 0.05;                  // particle density
    double baseLength;                      // base length of simulation box
    maxInteractions = 15;                   // maximum number of interactions per particle

    // Data structures.
    bool isIsotropic[nParticles];           // whether the potential of each particle is isotropic
    std::vector<double> boxSize;            // size of simulation box in each dimension

    // Work out base length of simulation box (particle diameter is one).
    if (dimension == 2) baseLength = std::pow((nParticles*M_PI)/(4.0*density), 1.0/2.0);
    else baseLength = std::pow((nParticles*M_PI)/(6.0*density), 1.0/3.0);

    for (unsigned int i=0;i<dimension;i++)
        boxSize.push_back(baseLength);

    // Add current directory to path.
    setenv("PYTHONPATH", ".", 1);

    // Initialise the Python interpreter.
    Py_Initialize();

    // Common and reuseable Python objects.
    PyObject *pArgs, *pClass, *pString, *pResult;

    // Load the module object.
    PyObject *pModule = PyImport_ImportModule("demo");

    // Create box size list.
    PyObject *pBoxSize = PyList_New(dimension);
    for (unsigned int i=0;i<dimension;i++)
        PyList_SetItem(pBoxSize, i, PyFloat_FromDouble(baseLength));

    // Reference to box object.
    pClass = PyObject_GetAttrString(pModule, "Box");

    // Create argument object.
    pArgs = Py_BuildValue("(O)", pBoxSize);

    // Instatiate box object.
    box = PyObject_CallObject(pClass, pArgs);

    // Reference to particle object.
    pClass = PyObject_GetAttrString(pModule, "Particle");

    // Create particle list.
    particles = PyList_New(nParticles);
    for (unsigned int i=0;i<nParticles;i++)
        PyList_SetItem(particles, i, PyObject_CallObject(pClass, NULL));

    // Reference to cell list object.
    pClass = PyObject_GetAttrString(pModule, "CellList");

    // Create argument object.
    pArgs = Py_BuildValue("(IOd)", dimension, pBoxSize, interactionRange);

    // Instatiate cell list object.
    cells = PyObject_CallObject(pClass, pArgs);

    // Reference to initalise object.
    pClass = PyObject_GetAttrString(pModule, "Initialise");

    // Instatiate initialise object.
    PyObject *initialise = PyObject_CallObject(pClass, NULL);

    // Create method name string.
    pString = PyUnicode_FromString("random");

    // Generate a random particle configuration.
    pResult = PyObject_CallMethodObjArgs(initialise, pString, particles, cells, box, NULL);

    // Reference to input/output object.
    pClass = PyObject_GetAttrString(pModule, "InputOutput");

    // Instatiate input/output object.
    io = PyObject_CallObject(pClass, NULL);

    // Create method name string.
    pString = PyUnicode_FromString("vmd_script");

    // Create VMD TcL script.
    pResult = PyObject_CallMethodObjArgs(io, pString, pBoxSize, NULL);

    // Reference to square-wellium object.
    pClass = PyObject_GetAttrString(pModule, "SquareWellium");

    // Create argument object.
    pArgs = Py_BuildValue("(OOOIdd)", box, particles, cells,
        maxInteractions, interactionEnergy, interactionRange);

    // Instatiate square-wellium object.
    squareWellium = PyObject_CallObject(pClass, pArgs);

    // C-style coordinate and orientation arrays.
    double coordinates[dimension*nParticles];
    double orientations[dimension*nParticles];

    // Initialise Python objects to store particle position and orientation.
    PyObject *pPosition, *pOrientation;

    // Copy coordinates and orientations into arrays.
    for (unsigned int i=0;i<nParticles;i++)
    {
        pPosition = PyObject_GetAttrString(PyList_GetItem(particles, i), "position");
        pOrientation = PyObject_GetAttrString(PyList_GetItem(particles, i), "orientation");

        for (unsigned int j=0;j<dimension;j++)
        {
            coordinates[dimension*i + j] = PyFloat_AsDouble(PyList_GetItem(pPosition, j));
            orientations[dimension*i + j] = PyFloat_AsDouble(PyList_GetItem(pOrientation, j));
        }

        // Set all particles as isotropic.
        isIsotropic[i] = true;
    }

    // Deallocate memory.
    Py_DECREF(pPosition);
    Py_DECREF(pOrientation);

    // Assign VMMC callback functions.
    vmmc::CallbackFunctions callbacks;
    callbacks.energyCallback = computeEnergy;
    callbacks.pairEnergyCallback = computePairEnergy;
    callbacks.interactionsCallback = computeInteractions;
    callbacks.postMoveCallback = applyPostMoveUpdates;

    // Initialise VMMC object.
    vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
        0.15, 0.2, 0.5, 0.5, maxInteractions, &boxSize[0], isIsotropic, false, callbacks);

    // Execute the simulation.
    for (unsigned int i=0;i<1000;i++)
    {
        // Increment simulation by 1000 Monte Carlo Sweeps.
        vmmc += 1000*nParticles;

        // Append particle coordinates to an xyz trajectory.
        if (i == 0) appendXyzTrajectory(true);
        else appendXyzTrajectory(false);

        // Report.
        printf("sweeps = %9.4e, energy = %5.4f\n", ((double) (i+1)*1000), getEnergy());
    }

    std::cout << "\nComplete!\n";

    // Deallocate memory.
    Py_DECREF(box);
    Py_DECREF(particles);
    Py_DECREF(cells);
    Py_DECREF(io);
    Py_DECREF(squareWellium);
    Py_DECREF(pArgs);
    Py_DECREF(pClass);
    Py_DECREF(pString);
    Py_DECREF(pResult);

    // Finalise the Python interpreter.
    Py_Finalize();

    // We're done!
    return (EXIT_SUCCESS);
}

// Function definitions.

double computeEnergy(unsigned int particle, const double* position, const double* orientation)
{
    // Create position and orientation objects.
    PyObject *pPosition = PyList_New(dimension);
    PyObject *pOrientation = PyList_New(dimension);

    // Copy coordinates and orientations into Python lists.
    for (unsigned int i=0;i<dimension;i++)
    {
        PyList_SetItem(pPosition, i, PyFloat_FromDouble(position[i]));
        PyList_SetItem(pOrientation, i, PyFloat_FromDouble(orientation[i]));
    }

    // Create method name string.
    PyObject *pString = PyUnicode_FromString("compute_energy");

    // Particle index.
    PyObject *pParticle = PyLong_FromLong(particle);

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(squareWellium,
        pString, pParticle, pPosition, pOrientation, NULL);

    // Store result.
    double energy = PyFloat_AsDouble(pResult);

    // Dellocate memory.
    Py_DECREF(pPosition);
    Py_DECREF(pOrientation);
    Py_DECREF(pParticle);
    Py_DECREF(pString);
    Py_DECREF(pResult);

    return energy;
}

double computePairEnergy(unsigned int particle1, const double* position1, const double* orientation1,
    unsigned int particle2, const double* position2, const double* orientation2)
{
    // Create position and orientation objects.
    PyObject *pPosition1 = PyList_New(dimension);
    PyObject *pPosition2 = PyList_New(dimension);
    PyObject *pOrientation1 = PyList_New(dimension);
    PyObject *pOrientation2 = PyList_New(dimension);

    // Copy coordinates and orientations into Python lists.
    for (unsigned int i=0;i<dimension;i++)
    {
        PyList_SetItem(pPosition1, i, PyFloat_FromDouble(position1[i]));
        PyList_SetItem(pPosition2, i, PyFloat_FromDouble(position2[i]));
        PyList_SetItem(pOrientation1, i, PyFloat_FromDouble(orientation1[i]));
        PyList_SetItem(pOrientation2, i, PyFloat_FromDouble(orientation2[i]));
    }

    // Create method name string.
    PyObject *pString = PyUnicode_FromString("compute_pair_energy");

    // Particle indices.
    PyObject *pParticle1 = PyLong_FromLong(particle1);
    PyObject *pParticle2 = PyLong_FromLong(particle2);

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(squareWellium, pString, pParticle1,
        pPosition1, pOrientation1, pParticle2, pPosition2, pOrientation2, NULL);

    // Store result.
    double energy = PyFloat_AsDouble(pResult);

    // Dellocate memory.
    Py_DECREF(pPosition1);
    Py_DECREF(pPosition2);
    Py_DECREF(pOrientation1);
    Py_DECREF(pOrientation2);
    Py_DECREF(pParticle1);
    Py_DECREF(pParticle2);
    Py_DECREF(pString);
    Py_DECREF(pResult);

    return energy;
}

unsigned int computeInteractions(unsigned int particle,
    const double* position, const double* orientation, unsigned int* interactions)
{
    // Create position, orientation, and interactions objects.
    PyObject *pPosition = PyList_New(dimension);
    PyObject *pOrientation = PyList_New(dimension);
    PyObject *pInteractions = PyList_New(maxInteractions);

    // Copy coordinates and orientations into Python lists.
    for (unsigned int i=0;i<dimension;i++)
    {
        PyList_SetItem(pPosition, i, PyFloat_FromDouble(position[i]));
        PyList_SetItem(pOrientation, i, PyFloat_FromDouble(orientation[i]));
    }

    // Initialise interaction list.
    for (unsigned int i=0;i<maxInteractions;i++)
        PyList_SetItem(pInteractions, i, PyLong_FromLong(0));

    // Create method name string.
    PyObject *pString = PyUnicode_FromString("compute_interactions");

    // Particle index.
    PyObject *pParticle = PyLong_FromLong(particle);

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(squareWellium, pString,
        pParticle, pPosition, pOrientation, pInteractions, NULL);

    // Store result.
    unsigned int nInteractions = PyLong_AsLong(pResult);

    // Copy interactions back into C array.
    for (unsigned int i=0;i<nInteractions;i++)
        interactions[i] = PyLong_AsLong(PyList_GetItem(pInteractions, i));

    // Dellocate memory.
    Py_DECREF(pPosition);
    Py_DECREF(pOrientation);
    Py_DECREF(pInteractions);
    Py_DECREF(pParticle);
    Py_DECREF(pString);
    Py_DECREF(pResult);

    return nInteractions;
}

void applyPostMoveUpdates(unsigned int particle, const double* position, const double* orientation)
{
    // Create position and orientation objects.
    PyObject *pPosition = PyList_New(dimension);
    PyObject *pOrientation = PyList_New(dimension);

    // Copy coordinates and orientations into Python lists.
    for (unsigned int i=0;i<dimension;i++)
    {
        PyList_SetItem(pPosition, i, PyFloat_FromDouble(position[i]));
        PyList_SetItem(pOrientation, i, PyFloat_FromDouble(orientation[i]));
    }

    // Create method name string.
    PyObject *pString = PyUnicode_FromString("apply_post_move_updates");

    // Particle index.
    PyObject *pParticle = PyLong_FromLong(particle);

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(squareWellium,
        pString, pParticle, pPosition, pOrientation, NULL);

    // Dellocate memory.
    Py_DECREF(pPosition);
    Py_DECREF(pOrientation);
    Py_DECREF(pParticle);
    Py_DECREF(pString);
    Py_DECREF(pResult);
}

double getEnergy()
{
    // Create method name string.
    PyObject *pString = PyUnicode_FromString("get_energy");

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(squareWellium, pString, NULL);

    // Store result.
    double energy = PyFloat_AsDouble(pResult);

    // Dellocate memory.
    Py_DECREF(pString);
    Py_DECREF(pResult);

    return energy;
}

void appendXyzTrajectory(bool clearFile)
{
    // Create method name string.
    PyObject *pString = PyUnicode_FromString("append_xyz_trajectory");

    // Box dimension.
    PyObject *pDimension = PyLong_FromLong(dimension);

    // Whether to wipe file.
    PyObject *pClearFile = PyBool_FromLong(clearFile);

    // Execute Python callback.
    PyObject *pResult = PyObject_CallMethodObjArgs(io, pString,
        pDimension, particles, pClearFile, NULL);

    // Dellocate memory.
    Py_DECREF(pString);
    Py_DECREF(pDimension);
    Py_DECREF(pClearFile);
    Py_DECREF(pResult);
}
