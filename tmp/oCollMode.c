static void oCollMode(dictionary *ini){

	/*
	 * SELECT METHODS
	 */
	void (*acc)()   			= select(ini,	"methods:acc",
												puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set,
                        puBoris3D1KETEST_set);

	void (*distr)() 			= select(ini,	"methods:distr",
												puDistr3D1split_set,
												puDistr3D1_set,
												puDistrND1_set,
												puDistrND0_set);


	void (*collide)() = select(ini,	"methods:mcc",
									mccCollissionsOff_set,
									mccConstCrossect_set,
									mccConstFreq_set,
									mccFunctionalCrossect_set);


	void (*extractEmigrants)()	= select(ini,	"methods:migrate",
												puExtractEmigrants3D_set,
												puExtractEmigrantsND_set,
                        						puExtractEmigrants3DOpen_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set
												//sSolver_set
											);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);
	mccNormalize(ini,units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
    	Grid *rhoObj = gAlloc(ini, SCALAR,mpiInfo); // for capMatrix - objects
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	void *solver = solverAlloc(ini, rho, phi, mpiInfo);
	MccVars *mccVars=mccAlloc(ini,units);

    	PincObject *obj = objoAlloc(ini,mpiInfo,units); // for capMatrix - objects
	//TODO: look into multigrid E,rho,rhoObj

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(ini, phi, mpiInfo);
	gSetBndSlicesE(ini, E, mpiInfo);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */

	pOpenH5(ini, pop, units, "pop");
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, rho_e, mpiInfo, units, units->chargeDensity, "rho_e");
	gOpenH5(ini, rho_i, mpiInfo, units, units->chargeDensity, "rho_i");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);
	xyCreateDataset(history,"/current/electrons/dataset");
	xyCreateDataset(history,"/current/ions/dataset");
	xyCreateDataset(history,"/potential/dataset");

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	 * INITIAL CONDITIONS
	 */

	//Compute capacitance matrix
	oComputeCapacitanceMatrix(obj, ini, mpiInfo);

	// Initalize particles
	pPosUniformCell(ini,rho,pop,rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// Perturb particles
	//pPosPerturb(ini, pop, mpiInfo);

	//add influx of new particles on boundary
	pPurgeGhost(pop, rho);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	pFillGhost(ini,rho,pop,rng);

	/*
	 * INITIALIZATION (E.g. half-step)
	 */

    // Clean objects from any charge first.
    gZero(rhoObj);                                   // for capMatrix - objects
    oCollectObjectCharge(pop, rhoObj, obj, mpiInfo); // for capMatrix - objects
    gZero(rhoObj);                                   // for capMatrix - objects

	// Get initial charge density
	distr(pop, rho,rho_e,rho_i);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

	// Get initial E-field

	solve(solver, rho, phi, mpiInfo);

	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);
	gBnd(E, mpiInfo);

  	//Boris parameters
  	int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

  	// add External E
  	puAddEext(ini, pop, E); // adds same value to whole grid

    gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	//-----------------------------------
	//- NEUTRALS - initialization
	//-----------------------------------

	Grid *rhoNeutral = gAlloc(ini, SCALAR,mpiInfo);
	gZero(rhoNeutral);
	gAdd(rhoNeutral,mccVars->nt);

	//-----------------------------------
	//- NEUTRALS - initialization - end
	//-----------------------------------

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i",n);
		tStart(t);

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		// Move particles
		// oRayTrace(pop, obj, deltaRho); <- do we need this still???
		puMove(pop); //puMove(pop, obj); Do not change functions such that PINC does
    	// not work in other run modes!
	    //neMove(neutralPop); // SPH neutrals

		/*
		*   Collisions
		*   Changes velocity component of some particles, not position.
		*/
		//collide(ini,rhoNeutral, pop, mccVars, rng,mpiInfo);

		//add influx of new particles on boundary
		pPurgeGhost(pop, rho);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		pFillGhost(ini,rho,pop,rng);

		// Check that no particle resides out-of-bounds (just for debugging)
		//pPosAssertInLocalFrame(pop, rho); //gives error with open boundary

        // Collect the charges on the objects.
        oCollectObjectCharge(pop, rhoObj, obj, mpiInfo);    // for capMatrix - objects

		// Original
		/*
		*   Collisions
		*   Changes velocity component of some particles, not position.
		*/
		collide(ini,rhoNeutral, pop, mccVars, rng,mpiInfo);

		// Compute charge density
		distr(pop, rho,rho_e,rho_i);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

        // Keep writing Rho here.

        // Add object charge to rho.
        gAddTo(rho, rhoObj);

        //gBnd(phi, mpiInfo);
        solve(solver, rho, phi, mpiInfo); // for capMatrix - objects

        // Second run with solver to account for charges
		oSweepBiasSin( obj, n );
        oApplyCapacitanceMatrix(rho, phi, obj, mpiInfo, units); // for capMatrix - objects

		solve(solver, rho, phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);
		gBnd(E, mpiInfo);

		// Apply external E
        puAddEext(ini, pop, E); // adds same value to whole grid

		// Accelerate particle and compute kinetic energy for step n
        acc(pop, E, T, S);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		if(n >= nTimeSteps-1000 && n%100 == 0 ){
			gWriteH5(rho, mpiInfo, (double) n);
			gWriteH5(rho_e, mpiInfo, (double) n);
			gWriteH5(rho_i, mpiInfo, (double) n);

			gWriteH5(phi, mpiInfo, (double) n);
		}

		pWriteEnergy(history,pop,(double)n,units);
		xyWrite(history,"/current/electrons/dataset",(double)n,units->current*obj->objectCurrent[0],MPI_SUM);
		xyWrite(history,"/current/ions/dataset",(double)n,units->current*obj->objectCurrent[1],MPI_SUM);
		xyWrite(history,"/potential/dataset",(double)n,units->potential*(*obj->bias),MPI_MAX);
	}

    tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */

	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(rho_e);
	gCloseH5(rho_i);

	gCloseH5(phi);
	gCloseH5(E);
    oCloseH5(obj); // for capMatrix - objects

	xyCloseH5(history);

  	// Free memory
	solverFree(solver);
	mccFreeVars(mccVars);
	gFree(rho);
	gFree(rho_e);
	gFree(rho_i);
	gFree(phi);
	free(S);
	free(T);

	gFree(E);

	pFree(pop);
	uFree(units);
    gFree(rhoObj); // for capMatrix - objects
    oFree(obj); // for capMatrix - objects

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);
}