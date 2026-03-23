static void mccMode(dictionary *ini){

	msg(STATUS, "start mcc Test Mode");

	/*
	* SELECT METHODS
	*/
	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
	puAcc3D1KE_set,
	puAccND1_set,
	puAccND1KE_set,
	puAccND0_set,
	puAccND0KE_set,
	puBoris3D1_set,
	puBoris3D1KE_set,
	puBoris3D1KETEST_set);

	void (*distr)() = select(ini,"methods:distr",	puDistr3D1split_set,
													puDistr3D1_set,
													puDistrND1_set,
													puDistrND0_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
	puExtractEmigrantsND_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set
												//sSolver_set
												);

	void (*collide)() = select(ini,	"methods:mcc",
									mccCollissionsOff_set,
									mccConstCrossect_set,
									mccConstFreq_set,
									mccFunctionalCrossect_set);

	//
	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	* INITIALIZE PINC VARIABLES
	*
	*done in "main" for "regular mode"
	*/
	Units *units=uAlloc(ini);
	uNormalize(ini, units);
	// normalize mcc input variables
	//must be done after uNormalize, and before defining mcc vars
	mccNormalize(ini,units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR, mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
	void *solver = solverAlloc(ini, rho, phi);

	/*
	* mcc specific variables
	*/

	MccVars *mccVars=mccAlloc(ini,units);

	// using Boris algo
	int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

	Grid *rhoNeutral = gAlloc(ini, SCALAR,mpiInfo);
	gZero(rhoNeutral);
	gAdd(rhoNeutral,mccVars->nt);

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(ini, phi, mpiInfo);

	//Set mgSolve
	//MgAlgo mgAlgo = getMgAlgo(ini);

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
	hid_t temperature = xyOpenH5(ini,"temperature");
	pCreateTemperatureDatasets(temperature,pop);
	hid_t probe = arrOpenH5(ini,"probe");
	xyzProbeCreateDatasets(probe,phi,mpiInfo);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	* INITIAL CONDITIONS
	*/

	// Initalize particles
	pPosUniformCell(ini,rho,pop,rng);
	pVelMaxwell(ini, pop, rng);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	/*
	* compute initial half-step
	*/

	// Get initial charge density
	distr(pop, rho, rho_e, rho_i); //two species
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

	// Get initial E-field
	solve(solver, rho, phi, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);

	// add External E
	puAddEext(ini, pop, E);

	// Advance velocities half a step
	gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	//Write initial h5 files
	gWriteH5(E, mpiInfo, 0.0);
	gWriteH5(rho, mpiInfo, 0.0);
	gWriteH5(phi, mpiInfo, 0.0);
	pWriteTemperature(temperature,pop,0.0,units,ini);
	pWriteEnergy(history,pop,0.0,units);
	xyzWriteProbe(probe, phi,mpiInfo);

	/*
	* TIME LOOP
	*/

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i of %i",n,nTimeSteps);

		// Check that no particle moves beyond a cell (mostly for debugging)
		//pVelAssertMax(pop,maxVel);

		tStart(t);
		// Move particles
		puMove(pop);
		// oRayTrace(pop, obj);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		/*
		*   Collisions
		*/
		collide(ini,rhoNeutral, pop, mccVars, rng,mpiInfo);

		// Check that no particle resides out-of-bounds (just for debugging)
		//pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		distr(pop, rho); //two species

		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);

		//gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		//solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		solve(solver, rho, phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);

		//gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);
		puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n
		acc(pop, E, T, S);
		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop); // kin_energy per species to pop
		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		if( n%500 == 0){
			// Example of writing another dataset to history.xy.h5
			// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
			//Write h5 files
			gWriteH5(E, mpiInfo, (double) n);
			gWriteH5(rho, mpiInfo, (double) n);
			//gWriteH5(rho_e, mpiInfo, (double) n);
			//gWriteH5(rho_i, mpiInfo, (double) n);
			gWriteH5(phi, mpiInfo, (double) n);
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		}

		if(n == 80000){
			pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5); //0.0001
		}

		pWriteTemperature(temperature,pop,(double)n,units,ini);
		//fillGridIndexes(phi);
		xyzWriteProbe(probe, phi,mpiInfo);
		pWriteEnergy(history,pop,(double)n,units);
	}

	//if(mpiInfo->mpiRank==0) 
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
	xyCloseH5(history);
	xyCloseH5(temperature);
	arrCloseH5(probe);

	// Free memory
	mccFreeVars(mccVars);
	solverFree(solver);
	gFree(rho);
	gFree(rho_e);
	gFree(rho_i);
	gFree(phi);
	gFree(E);
	pFree(pop);
	free(S);
	free(T);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);
}