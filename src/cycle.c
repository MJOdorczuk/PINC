/**
 * @file		cycle.c
 * @brief		PINC cycle.
 * @author		Michal Jan Odorczuk <michaljo@uio.no>
 *
 * PINC cycle implementation.
 */


#include "core.h"
#include "cycle.h"

// TODO: Implement the general cycle and decompose all the modes into smaller
// functions passed to the cycle.
// Modes to consider:
// main.c: regular, BorisTestMode
// multigrid.c: mgMode, mgModeErrorScaling
// spectral.c: sMode
// object.c: oMode
// collisions.c: mccMode, oCollMode, neutTest
//
// Run consists of the pre-cycle initialisation, the cycle, and the post-cycle
// cleanup.

/*************************************************
 *		Auxiliary functions
 ************************************************/

static void selectMethods(dictionary *ini, struct methods *methods){
    methods->acc = select(ini, "methods:acc",
                                            puAcc3D1_set,
                                            puAcc3D1KE_set,
                                            puAccND1_set,
                                            puAccND1KE_set,
                                            puAccND0_set,
                                            puAccND0KE_set,
                                            puBoris3D1_set,
                                            puBoris3D1KE_set,
                                            puBoris3D1KETEST_set);
    methods->distr = select(ini, "methods:distr",
                                            puDistr3D1split_set,
                                            puDistr3D1_set,
                                            puDistrND1_set,
                                            puDistrND0_set);
    methods->solverInterface = select(ini, "methods:poisson",
                                            mgSolver_set,
                                            sSolver_set);
    methods->extractEmigrants = select(ini,	"methods:migrate",
                                            puExtractEmigrants3D_set,
                                            puExtractEmigrantsND_set,
                                            puExtractEmigrants3DOpen_set);
    methods->collide = select(ini, "methods:mcc",
                                            mccCollissionsOff_set,
                                            mccConstCrossect_set,
                                            mccConstFreq_set,
                                            mccFunctionalCrossect_set);
    methods->solve = NULL;
    methods->solverAlloc = NULL;
    methods->solverFree = NULL;
    solverInterface(&solve, &solverAlloc, &solverFree);
}

/*************************************************
 *		Cycle
 ************************************************/

static void cycle(dictionary *ini){
	// TODO: Implement the cycle.

    Methods methods;
    // Not present for mgMode, mgModeErrorScaling and sMode
    selectMethods(ini, &methods);

    // PINC variables
    // Set grid and ghost layers
    // Set random number seeds
    // Prepare files for writing
    // Open files
    // Write initial data
    // Run the cycle
    // Close files
    // Free memory
}