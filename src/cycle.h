/**
 * @file		cycle.h
 * @brief		PINC cycle.
 * @author		Michal Jan Odorczuk <michaljo@uio.no>
 *
 * PINC cycle implementation.
 */

#ifndef CYCLE_H
#define CYCLE_H

#include "core.h"
#include "collisions.h"
#include "pusher.h"
#include "object.h"
#include "population.h"
#include "grid.h"
#include "units.h"
#include "multigrid.h"
#include "spectral.h"

typedef struct {
    void (*acc)();
    void (*distr)();
    void (*solverInterface)();
    void (*extractEmigrants)();
    void (*collide)();
    void (*solve)();
    void *(*solverAlloc)();
    void (*solverFree)();
} Methods;

typedef struct {
    bool isCollisional;
    bool isBoris;
    bool isObject;
    bool usesAcc;
    bool usesDistr;
    bool usesSolver;
    bool usesExtractEmigrants;
    bool usesCollide;
} PINCflags;

typedef struct {
    Units *units;
    MpiInfo *mpiInfo;
    Population *pop;
    Grid *E;
    Grid *rho;
    Grid *rho_e;
    Grid *rho_i;
    Grid *rhoObj;
    Grid *phi;
    PincObject *obj;
    hid_t history;
    void *solver;
    MccVars *mccVars;
    Grid *rhoNeutral;
    double *S;
    double *T;
    gsl_rng *rngSync;
    gsl_rng *rng;
} PINCvariables;

// TODO: Add the description
void cycle(dictionary *ini);

#endif