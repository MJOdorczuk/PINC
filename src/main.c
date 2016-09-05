/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 */

#include "core.h"
#include "pusher.h"
#include "multigrid.h"

void regular(dictionary *ini);
funPtr regular_set(dictionary *ini){ return regular; }

void mgRun(dictionary *ini);
funPtr mgRun_set(dictionary *ini){ return mgRun; }

void mgErrorScaling(dictionary *ini);
funPtr mgErrorScaling_set(dictionary *ini){ return mgErrorScaling; }


int main(int argc, char *argv[]){

	/*
	 * INITIALIZE PINC
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv); // No printing before this
	msg(STATUS, "PINC %s started.", VERSION);    // Needs MPI
	MPI_Barrier(MPI_COMM_WORLD);
	parseIndirectInput(ini);

	/*
	 * CHOOSE PINC RUN MODE
	 */
	void (*run)() = select(ini,"methods:mode",regular_set,mgRun_set, mgErrorScaling_set);
	run(ini);

	/*
	 * FINALIZE PINC
	 */
	iniClose(ini);
	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}

void regular(dictionary *ini){

	/*
	 * SELECT METHODS
	 */
	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set);
	void (*distr)() = select(ini,"methods:distr",	puDistr3D1_set,
													puDistrND1_set,
													puDistrND0_set);
	void (*solve)() = select(ini,"methods:poisson", mgSolve_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
																puExtractEmigrantsND_set);

	// char *str;
	//
	// str = iniGetStr("methods:acc");
	// void (*acc)() = NULL;
	// if(!strcmp(str,"puAcc3D1")) acc = puAcc3D1_set();
	// if(!strcmp(str,"puAcc3D1KE")) acc = puAcc3D1KE_set();
	// if(acc==NULL) msg(ERROR,"methods:acc=%s is an invalid option")

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	//Set mgSolve
	MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);

	/*
	 * PREPARE FILES FOR WRITING
	 */
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	pOpenH5(ini, pop, "pop");
	gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);

	/*
	 * INITIAL CONDITIONS
	 */

	// Initalize particles
	// pPosUniform(ini, pop, mpiInfo, rngSync);
	pPosLattice(ini, pop, mpiInfo);
	pVelZero(pop);

	// Perturb particles
	pPosPerturb(ini, pop, mpiInfo);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	/*
	 * INITIALIZATION (E.g. half-step)
	 */

	// Get initial charge density
	distr(pop, rho);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

	// Get initial E-field
	solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);

	// Advance velocities half a step
	gMul(E, 0.5);
	acc(pop, E);
	gMul(E, 2.0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(rank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i",n);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,1.0);

		tStart(t);

		// Move particles
		puMove(pop);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		distr(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);

		gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);

		// Accelerate particle and compute kinetic energy for step n
		acc(pop, E);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		//Write h5 files
		gWriteH5(E, mpiInfo, (double) n);
		gWriteH5(rho, mpiInfo, (double) n);
		gWriteH5(phi, mpiInfo, (double) n);
		// pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		pWriteEnergy(history,pop,(double)n);

	}

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */
	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(phi);
	gCloseH5(E);
	xyCloseH5(history);

	// Free memory
	mgFree(mgRho);
	mgFree(mgPhi);
	mgFree(mgRes);
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);

	gsl_rng_free(rngSync);

}
