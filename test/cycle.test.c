/**
 * @file		cycle.test.c
 * @brief		Unit tests for cycle.c
 * @author		Michal Jan Odorczuk <michaljo@uio.no>,
 */

#include "core.h"
#include "test.h"
#include "iniparser.h"

/*
 * Include cycle.c directly so the test translation unit can call its internal
 * static helpers without changing production linkage.
 */
#define cycle cycleTestTranslationUnitCycle
#include "../src/cycle.c"
#undef cycle

static void accSentinel(){}
static void distrSentinel(){}
static void extractEmigrantsSentinel(){}
static void solverInterfaceSentinel(){}
static void collideSentinel(){}

static int dummySolverStorage = 0;
static int dummySolverAllocCalls = 0;
static int dummySolverFreeCalls = 0;
static const dictionary *dummySolverIniArg = NULL;
static Grid *dummySolverRhoArg = NULL;
static Grid *dummySolverPhiArg = NULL;
static MpiInfo *dummySolverMpiInfoArg = NULL;
static void *dummySolverFreedArg = NULL;

static void resetDummySolverHooks(){
	dummySolverAllocCalls = 0;
	dummySolverFreeCalls = 0;
	dummySolverIniArg = NULL;
	dummySolverRhoArg = NULL;
	dummySolverPhiArg = NULL;
	dummySolverMpiInfoArg = NULL;
	dummySolverFreedArg = NULL;
}

static void *dummySolverAlloc(const dictionary *ini, Grid *rho, Grid *phi,
	const MpiInfo *mpiInfo){
	dummySolverAllocCalls++;
	dummySolverIniArg = ini;
	dummySolverRhoArg = rho;
	dummySolverPhiArg = phi;
	dummySolverMpiInfoArg = (MpiInfo *)mpiInfo;
	return &dummySolverStorage;
}

static void dummySolverFree(void *solver){
	dummySolverFreeCalls++;
	dummySolverFreedArg = solver;
}

static dictionary *makeSelectPINCvariablesIni(){
	dictionary *ini = iniGetDummy();

	iniparser_set(ini, "methods:normalization", "SI");
	iniparser_set(ini, "time:timeStep", "1.0");
	iniparser_set(ini, "population:nSpecies", "2");
	iniparser_set(ini, "population:nParticles", "8,8");
	iniparser_set(ini, "population:nAlloc", "16,16");
	iniparser_set(ini, "population:charge", "-1,1");
	iniparser_set(ini, "population:mass", "1,1836");
	iniparser_set(ini, "population:density", "1,1");
	iniparser_set(ini, "population:thermalVelocity", "1,1");
	iniparser_set(ini, "population:drift", "0,0,0,0,0,0");
	iniparser_set(ini, "fields:BExt", "0,0,0");
	iniparser_set(ini, "fields:EExt", "0,0,0");

	return ini;
}

static void addCollisionalIni(dictionary *ini){
	iniparser_set(ini, "collisions:artificialLoss", "1.0");
	iniparser_set(ini, "collisions:nSpeciesNeutral", "1");
	iniparser_set(ini, "collisions:neutralDrift", "0,0,0");
	iniparser_set(ini, "collisions:numberDensityNeutrals", "2.5");
	iniparser_set(ini, "collisions:thermalVelocityNeutrals", "0.5");
	iniparser_set(ini, "collisions:realElectronMass", "1.0");
	iniparser_set(ini, "collisions:collFrqCex", "0.2");
	iniparser_set(ini, "collisions:collFrqIonElastic", "0.3");
	iniparser_set(ini, "collisions:collFrqElectronElastic", "0.4");
	iniparser_set(ini, "collisions:CEX_a", "0.1");
	iniparser_set(ini, "collisions:CEX_b", "0.2");
	iniparser_set(ini, "collisions:ion_elastic_a", "0.3");
	iniparser_set(ini, "collisions:ion_elastic_b", "0.4");
	iniparser_set(ini, "collisions:electron_a", "0.5");
	iniparser_set(ini, "collisions:electron_b", "0.6");
	iniparser_set(ini, "collisions:electronEnergyMethod", "test");
}

static int gridHasUniformValue(const Grid *grid, double expected, double tol){
	long int nValues = grid->sizeProd[grid->rank];
	for(long int i=0; i<nValues; i++){
		if(fabs(grid->val[i] - expected) > tol){
			return 0;
		}
	}
	return 1;
}

static int testSelectMethodsChoosesConfiguredMethods(){

	dictionary *ini = iniGetDummy();
	Methods methods = {0};
	PINCflags flags = {0};

	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1,1,1");
	iniparser_set(ini, "methods:acc", "puAcc3D1");
	iniparser_set(ini, "methods:distr", "puDistr3D1split");
	iniparser_set(ini, "methods:migrate", "puExtractEmigrants3D");
	iniparser_set(ini, "methods:poisson", "mgSolver");

	selectMethods(ini, &methods, &flags);

	utAssert(methods.acc == (void (*)())puAcc3D1,
		"expected methods:acc to select puAcc3D1");
	utAssert(methods.distr == (void (*)())puDistr3D1split,
		"expected methods:distr to select puDistr3D1split");
	utAssert(methods.extractEmigrants == (void (*)())puExtractEmigrants3D,
		"expected methods:migrate to select puExtractEmigrants3D");
	utAssert(methods.solverInterface == (void (*)())mgSolver_set(ini),
		"expected methods:poisson to select the multigrid solver interface");
	utAssert(methods.solve == (void (*)())mgSolve,
		"expected mgSolver interface to set mgSolve");
	utAssert(methods.solverAlloc == (void *(*)())mgAllocSolver,
		"expected mgSolver interface to set mgAllocSolver");
	utAssert(methods.solverFree == (void (*)())mgFreeSolver,
		"expected mgSolver interface to set mgFreeSolver");

	iniparser_freedict(ini);
	return 0;
}

static int testSelectMethodsChoosesCollisionMethodWhenEnabled(){

	dictionary *ini = iniGetDummy();
	Methods methods = {0};
	PINCflags flags = {0};
	flags.isCollisional = true;

	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1,1,1");
	iniparser_set(ini, "population:nSpecies", "2");
	iniparser_set(ini, "methods:acc", "puAcc3D1");
	iniparser_set(ini, "methods:distr", "puDistr3D1split");
	iniparser_set(ini, "methods:migrate", "puExtractEmigrants3D");
	iniparser_set(ini, "methods:poisson", "mgSolver");
	iniparser_set(ini, "methods:mcc", "mccConstFreq");

	selectMethods(ini, &methods, &flags);

	utAssert(methods.collide == (void (*)())mccConstFreq_set(ini),
		"expected methods:mcc to select mccConstFreq");

	iniparser_freedict(ini);
	return 0;
}

static int testSelectMethodsPreservesCustomSelections(){

	dictionary *ini = iniGetDummy();
	Methods methods = {
		.acc = accSentinel,
		.distr = distrSentinel,
		.solverInterface = solverInterfaceSentinel,
		.extractEmigrants = extractEmigrantsSentinel,
		.collide = collideSentinel,
		.solve = accSentinel,
		.solverAlloc = (void *(*)())solverInterfaceSentinel,
		.solverFree = solverInterfaceSentinel
	};
	PINCflags flags = {0};
	flags.customAcc = true;
	flags.customDistr = true;
	flags.customExtractEmigrants = true;
	flags.customSolver = true;

	selectMethods(ini, &methods, &flags);

	utAssert(methods.acc == accSentinel,
		"custom acc should not be overwritten");
	utAssert(methods.distr == distrSentinel,
		"custom distr should not be overwritten");
	utAssert(methods.extractEmigrants == extractEmigrantsSentinel,
		"custom extractEmigrants should not be overwritten");
	utAssert(methods.solverInterface == solverInterfaceSentinel,
		"custom solver interface should not be overwritten");
	utAssert(methods.collide == collideSentinel,
		"non-collisional run should not change collide");
	utAssert(methods.solve == NULL,
		"custom solver path currently clears solve");
	utAssert(methods.solverAlloc == NULL,
		"custom solver path currently clears solverAlloc");
	utAssert(methods.solverFree == NULL,
		"custom solver path currently clears solverFree");

	iniparser_freedict(ini);
	return 0;
}

static int testSelectPINCvariablesAllocatesCoreState(){

	dictionary *ini = makeSelectPINCvariablesIni();
	Methods methods = {
		.solverAlloc = (void *(*)())dummySolverAlloc,
		.solverFree = dummySolverFree
	};
	PINCflags flags = {.noObject = true};
	PINCvariables vars = {0};
	gsl_rng *expectedRng = NULL;

	resetDummySolverHooks();
	selectPINCvariables(ini, &vars, &flags, &methods);

	utAssert(vars.units != NULL, "expected units to be allocated");
	utAssert(vars.mpiInfo != NULL, "expected mpi info to be allocated");
	utAssert(vars.pop != NULL, "expected population to be allocated");
	utAssert(vars.E != NULL, "expected electric field grid to be allocated");
	utAssert(vars.rho != NULL, "expected rho grid to be allocated");
	utAssert(vars.rho_e != NULL, "expected rho_e grid to be allocated");
	utAssert(vars.rho_i != NULL, "expected rho_i grid to be allocated");
	utAssert(vars.phi != NULL, "expected phi grid to be allocated");
	utAssert(vars.solver == &dummySolverStorage,
		"expected solver pointer returned by solverAlloc to be stored");
	utAssert(vars.rhoObj == NULL, "expected rhoObj to remain NULL when noObject is true");
	utAssert(vars.obj == NULL, "expected object to remain NULL when noObject is true");
	utAssert(vars.mccVars == NULL,
		"expected mccVars to remain NULL when collisions are disabled");
	utAssert(vars.rhoNeutral == NULL,
		"expected rhoNeutral to remain NULL when collisions are disabled");
	utAssert(vars.S == NULL, "expected S to remain NULL when Boris is disabled");
	utAssert(vars.T == NULL, "expected T to remain NULL when Boris is disabled");
	utAssert(vars.rngSync != NULL, "expected synchronized RNG to be allocated");
	utAssert(vars.rng != NULL, "expected RNG to be allocated");
	utAssert(vars.pop->nSpecies == 2, "expected population to use overridden species count");

	utAssert(dummySolverAllocCalls == 1, "expected solverAlloc to be called exactly once");
	utAssert(dummySolverIniArg == ini, "expected solverAlloc to receive the same ini");
	utAssert(dummySolverRhoArg == vars.rho,
		"expected solverAlloc to receive the allocated rho grid");
	utAssert(dummySolverPhiArg == vars.phi,
		"expected solverAlloc to receive the allocated phi grid");
	utAssert(dummySolverMpiInfoArg == vars.mpiInfo,
		"expected solverAlloc to receive the allocated mpi info");

	expectedRng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(expectedRng, vars.mpiInfo->mpiRank + 1);
	utAssert(gsl_rng_get(vars.rng) == gsl_rng_get(expectedRng),
		"expected rng to be seeded from mpiRank + 1");
	gsl_rng_free(expectedRng);

	freePINCvariables(&vars, &flags, &methods);

	utAssert(dummySolverFreeCalls == 1, "expected solverFree to be called exactly once");
	utAssert(dummySolverFreedArg == &dummySolverStorage,
		"expected solverFree to receive the stored solver pointer");

	iniparser_freedict(ini);
	return 0;
}

static int testSelectPINCvariablesAllocatesBorisBuffers(){

	dictionary *ini = makeSelectPINCvariablesIni();
	Methods methods = {
		.solverAlloc = (void *(*)())dummySolverAlloc,
		.solverFree = dummySolverFree
	};
	PINCflags flags = {
		.noObject = true,
		.isBoris = true
	};
	PINCvariables vars = {0};

	resetDummySolverHooks();
	selectPINCvariables(ini, &vars, &flags, &methods);

	utAssert(vars.S != NULL, "expected S to be allocated for Boris runs");
	utAssert(vars.T != NULL, "expected T to be allocated for Boris runs");
	for(int i=0; i<3*vars.pop->nSpecies; i++){
		vars.S[i] = (double)i;
		vars.T[i] = (double)(i + 1);
	}

	freePINCvariables(&vars, &flags, &methods);

	utAssert(dummySolverFreeCalls == 1, "expected Boris path to still free solver state");

	iniparser_freedict(ini);
	return 0;
}

static int testSelectPINCvariablesAllocatesCollisionalState(){

	dictionary *ini = makeSelectPINCvariablesIni();
	Methods methods = {
		.solverAlloc = (void *(*)())dummySolverAlloc,
		.solverFree = dummySolverFree
	};
	PINCflags flags = {
		.noObject = true,
		.isCollisional = true
	};
	PINCvariables vars = {0};

	addCollisionalIni(ini);
	resetDummySolverHooks();
	selectPINCvariables(ini, &vars, &flags, &methods);

	utAssert(vars.mccVars != NULL, "expected collisional runs to allocate mccVars");
	utAssert(vars.rhoNeutral != NULL,
		"expected collisional runs to allocate rhoNeutral");
	utAssert(gridHasUniformValue(vars.rhoNeutral, vars.mccVars->nt, 1e-12),
		"expected rhoNeutral to be initialized from mccVars->nt");

	freePINCvariables(&vars, &flags, &methods);

	utAssert(dummySolverFreeCalls == 1,
		"expected collisional path to still free solver state");

	iniparser_freedict(ini);
	return 0;
}

// All tests for cycle.c is contained in this function
void testCycle(){
	utRun(&testSelectMethodsChoosesConfiguredMethods);
	utRun(&testSelectMethodsChoosesCollisionMethodWhenEnabled);
	utRun(&testSelectMethodsPreservesCustomSelections);
	utRun(&testSelectPINCvariablesAllocatesCoreState);
	utRun(&testSelectPINCvariablesAllocatesBorisBuffers);
	utRun(&testSelectPINCvariablesAllocatesCollisionalState);
}
 