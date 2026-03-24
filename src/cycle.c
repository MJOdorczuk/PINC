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

static void defaultPINCFlags(PINCflags *flags){
    flags->isCollisional = false;
    flags->isBoris = false;
    flags->isObject = true;
    flags->usesAcc = true;
    flags->usesDistr = true;
    flags->usesSolver = true;
    flags->usesExtractEmigrants = true;
    flags->usesCollide = true;
}

static void selectMethods(dictionary *ini, Methods *methods, PINCflags *flags){
    if(flags->usesAcc) {
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
    }
    if(flags->usesDistr) {
        methods->distr = select(ini, "methods:distr",
                                            puDistr3D1split_set,
                                            puDistr3D1_set,
                                            puDistrND1_set,
                                            puDistrND0_set);
    }
    if(flags->usesExtractEmigrants) {
        methods->extractEmigrants = select(ini,	"methods:migrate",
                                            puExtractEmigrants3D_set,
                                            puExtractEmigrantsND_set,
                                            puExtractEmigrants3DOpen_set);
    }
    if(flags->isCollisional){
        methods->collide = select(ini, "methods:mcc",
                                            mccCollissionsOff_set,
                                            mccConstCrossect_set,
                                            mccConstFreq_set,
                                            mccFunctionalCrossect_set);
    }
    methods->solve = NULL;
    methods->solverAlloc = NULL;
    methods->solverFree = NULL;
    if(flags->usesSolver) {
        // TODO: does it need to be a part of methods?
        methods->solverInterface = select(ini, "methods:poisson",
                                            mgSolver_set,
                                            sSolver_set);
        methods->solverInterface(&methods->solve, &methods->solverAlloc, &methods->solverFree);
    }
}

// TODO: is it the best way? Maybe it should be passed from outside?
static bool isCollisional(dictionary *ini){
    char *mode = iniGetStr(ini, "methods:mode");
    if(strcmp(mode, "mccMode") == 0 
    || strcmp(mode, "oCollMode") == 0 
    || strcmp(mode, "neutTest") == 0){
        free(mode);
        return true;
    }
    free(mode);
    return false;
}

static void selectPINCvariables(dictionary *ini, PINCvariables *vars,
                                PINCflags *flags, Methods *methods){
    vars->units = uAlloc(ini);
    vars->mpiInfo = gAllocMpi(ini);
    vars->pop = pAlloc(ini, vars->mpiInfo);
    vars->E = gAlloc(ini, VECTOR, vars->mpiInfo);
    vars->rho = gAlloc(ini, SCALAR, vars->mpiInfo);
    vars->rho_e = gAlloc(ini, SCALAR, vars->mpiInfo);
    vars->rho_i = gAlloc(ini, SCALAR, vars->mpiInfo);
    vars->phi = gAlloc(ini, SCALAR, vars->mpiInfo);
    vars->solver = methods->solverAlloc(ini, vars->rho, vars->phi, vars->mpiInfo);
    if (flags->isObject){
        vars->rhoObj = gAlloc(ini, SCALAR, vars->mpiInfo);
        vars->obj = objoAlloc(ini, vars->mpiInfo, vars->units);
    }
    if (flags->isCollisional){
        vars->mccVars = mccAlloc(ini, vars->units);
        vars->rhoNeutral = gAlloc(ini, SCALAR, vars->mpiInfo);
        gZero(vars->rhoNeutral);
        gAdd(vars->rhoNeutral, vars->mccVars->nt);
    }
    if (flags->isBoris){
        vars->S = (double*)malloc((3)*(vars->pop->nSpecies)*sizeof(double));
        vars->T = (double*)malloc((3)*(vars->pop->nSpecies)*sizeof(double));
    }
    vars->rngSync = gsl_rng_alloc(gsl_rng_mt19937);
    vars->rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(vars->rng, vars->mpiInfo->mpiRank+1); // Seed needs to be >=1

}

static void freePINCvariables(PINCvariables *vars, PINCflags *flags, Methods *methods){
    gsl_rng_free(vars->rngSync);
    gsl_rng_free(vars->rng);
    if (flags->isBoris){
        free(vars->S);
        free(vars->T);
    }
    if (flags->isCollisional){
        gFree(vars->rhoNeutral);
        mccFreeVars(vars->mccVars);
    }
    if (flags->isObject){
        oFree(vars->obj);
        gFree(vars->rhoObj);
    }
    methods->solverFree(vars->solver);
    gFree(vars->phi);
    gFree(vars->rho_i);
    gFree(vars->rho_e);
    gFree(vars->rho);
    gFree(vars->E);
    pFree(vars->pop);
    gFreeMpi(vars->mpiInfo);
    uFree(vars->units);
}

// TODO: check if this conforms to all the modes
static void setGridAndGhostLayers(dictionary *ini, PINCvariables *vars, Methods *methods){
    gCreateNeighborhood(ini, vars->mpiInfo, vars->rho);
    gSetBndSlices(ini, vars->phi, vars->mpiInfo);
    gSetBndSlicesE(ini, vars->E, vars->mpiInfo);
    // Initalize particles
    pPosUniformCell(ini, vars->rho, vars->pop, vars->rng);
    // TODO: do we need it?
    // double maxVel = iniGetDouble(ini, "population:maxVel");
    // Perturb particles
    //pPosPerturb(ini, vars->pop, vars->mpiInfo);
    //add influx of new particles on boundary
    pPurgeGhost(vars->pop, vars->rho);
    methods->extractEmigrants(vars->pop, vars->mpiInfo);
    puMigrate(vars->pop, vars->mpiInfo, vars->rho);
    pFillGhost(ini, vars->rho, vars->pop, vars->rng);
}

// TODO: Should be adjustable for different modes/based on ini
static void openFiles(dictionary *ini, PINCvariables *vars){
    Units *units = vars->units;
    pOpenH5(ini, vars->pop, units, "pop");
    gOpenH5(ini, vars->rho, vars->mpiInfo, units, units->chargeDensity, "rho");
    gOpenH5(ini, vars->rho_e, vars->mpiInfo, units, units->chargeDensity, "rho_e");
    gOpenH5(ini, vars->rho_i, vars->mpiInfo, units, units->chargeDensity, "rho_i");
    gOpenH5(ini, vars->phi, vars->mpiInfo, units, units->potential, "phi");
    gOpenH5(ini, vars->E, vars->mpiInfo, units, units->eField, "E");
    vars->history = xyOpenH5(ini, "history");
    pCreateEnergyDatasets(vars->history, vars->pop);
    xyCreateDataset(vars->history, "/current/electrons/dataset");
    xyCreateDataset(vars->history, "/current/ions/dataset");
    xyCreateDataset(vars->history, "/potential/dataset");
}

static void closeFiles(PINCvariables *vars){
    pCloseH5(vars->pop);
    gCloseH5(vars->rho);
    gCloseH5(vars->rho_e);
    gCloseH5(vars->rho_i);
    gCloseH5(vars->phi);
    gCloseH5(vars->E);
    xyCloseH5(vars->history);
}

static void setInitialConditions(dictionary *ini, PINCvariables *vars, 
    PINCflags *flags, Methods *methods){
    if (flags->isObject){
        oComputeCapacitanceMatrix(vars->obj, ini, vars->mpiInfo);
    }
    if (flags->isCollisional){
        gZero(vars->rhoNeutral);
        gAdd(vars->rhoNeutral, vars->mccVars->nt);
    }
    if (flags->isObject){
        gZero(vars->rhoObj);
        oCollectObjectCharge(vars->pop, vars->rhoObj, vars->obj, vars->mpiInfo);
        gZero(vars->rhoObj);
    }
    methods->distr(vars->pop, vars->rho, vars->rho_e, vars->rho_i);
    gHaloOp(addSlice, vars->rho, vars->mpiInfo, FROMHALO);
    gHaloOp(addSlice, vars->rho_e, vars->mpiInfo, FROMHALO);
    gHaloOp(addSlice, vars->rho_i, vars->mpiInfo, FROMHALO);
    methods->solve(vars->solver, vars->rho, vars->phi, vars->mpiInfo);
    gFinDiff1st(vars->phi, vars->E);
    gHaloOp(setSlice, vars->E, vars->mpiInfo, TOHALO);
    gMul(vars->E, -1.);
    gBnd(vars->E, vars->mpiInfo);
    if (flags->isBoris){
        puAddEext(ini, vars->pop, vars->E);
        gMul(vars->E, 0.5);
        puGet3DRotationParameters(ini, vars->T, vars->S, 0.5);
        methods->acc(vars->pop, vars->E, vars->T, vars->S);
        gMul(vars->E, 2.0);
        puGet3DRotationParameters(ini, vars->T, vars->S, 1.0);
    }
    if (flags->isCollisional){
        gZero(vars->rhoNeutral);
        gAdd(vars->rhoNeutral, vars->mccVars->nt);
    }
}

/*************************************************
 *		Cycle
 ************************************************/

void cycle(dictionary *ini){
	// TODO: Implement the cycle.

    Methods methods;
    // TODO: Should be passed as an argument
    PINCflags flags;
    defaultPINCFlags(&flags);
    flags.isCollisional = isCollisional(ini);
    // Not present for mgMode, mgModeErrorScaling and sMode
    selectMethods(ini, &methods, &flags);

    // PINC variables
    PINCvariables vars;
    selectPINCvariables(ini, &vars, &flags, &methods);
    // Set grid and ghost layers
    setGridAndGhostLayers(ini, &vars, &methods);
    // Prepare files for writing
    // Open files
    openFiles(ini, &vars);
    // Write initial data
    setInitialConditions(ini, &vars, &flags, &methods);
    // Run the cycle
    // Close files
    closeFiles(&vars);
    // Free memory
    freePINCvariables(&vars, &flags, &methods);
}