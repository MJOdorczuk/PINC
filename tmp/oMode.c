static void oMode(dictionary *ini){

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

	void (*extractEmigrants)()	= select(ini,	"methods:migrate",
												puExtractEmigrants3D_set,
												puExtractEmigrantsND_set,
                        						puExtractEmigrants3DOpen_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set,
												sSolver_set);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
    Grid *rhoObj = gAlloc(ini, SCALAR,mpiInfo);     // for capMatrix - objects
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);

	void *solver = solverAlloc(ini, rho, phi, mpiInfo);

    PincObject *obj = objoAlloc(ini,mpiInfo,units); // for capMatrix - objects
    //TODO: look into multigrid E,rho,rhoObj

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

  	// Setting Boundary slices
  	gSetBndSlices(ini, phi, mpiInfo);
	//gSetBndSlices(ini, solver->res, mpiInfo);
	//gSetBndSlices(ini, rho, mpiInfo);
	gSetBndSlicesE(ini, E, mpiInfo);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */

	pOpenH5(ini, pop, units, "pop");
	//double denorm = units->potential;
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, rho_e, mpiInfo, units, units->chargeDensity, "rho_e");
	gOpenH5(ini, rho_i, mpiInfo, units, units->chargeDensity, "rho_i");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");
  // oOpenH5(ini, obj, mpiInfo, units, 1, "test");
  // oReadH5(obj, mpiInfo);


    gOpenH5(ini, rhoObj, mpiInfo, units, units->chargeDensity, "rhoObj"); // for capMatrix - objects
    //oOpenH5(ini, obj, mpiInfo, units, units->chargeDensity, "object"); // for capMatrix - objects
    //oReadH5(obj->domain, mpiInfo, "Object");



    //Count the number of objects and fill the lookup tables.
    // This is done in oAlloc now....

    //oFillLookupTables(obj,mpiInfo);
    // Find all the object nodes which are part of the object surface.
    //oFindObjectSurfaceNodes(obj, mpiInfo);


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
	//pPosUniform(ini, pop, mpiInfo, rngSync);
	//pPosLattice(ini, pop, mpiInfo);
	pPosUniformCell(ini,rho,pop,rng);
	//pVelZero(pop);
	//pVelMaxwell(ini, pop, rng);
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
    //gWriteH5(rho, mpiInfo, (double) 0);

	// Get initial E-field

    //gBnd(phi, mpiInfo);
	solve(solver, rho, phi, mpiInfo);
	//gNeutralizeGrid(phi, mpiInfo);
	//gBnd(phi, mpiInfo);
    //gWriteH5(phi, mpiInfo, (double) 0);
    //pWriteH5(pop, mpiInfo, (double) 0, (double)0+0.5);

	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);
	gBnd(E, mpiInfo);

    //Boris parameters
    int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

    // add External E
	//gZero(E); // for testing Boris
	//gAddTo(Ext); //needs grid definition of Eext
  	puAddEext(ini, pop, E); // adds same value to whole grid

  	gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		long int totPs0 = (pop->iStop[0]- pop->iStart[0]); //debug
		long int totPs1 = (pop->iStop[1]- pop->iStart[1]);
		MPI_Allreduce(MPI_IN_PLACE, &totPs0, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &totPs1, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		msg(STATUS,"Computing time-step %i",n);
        msg(STATUS, "Nr. of particles s=0 %i: ",totPs0);
		msg(STATUS, "Nr. of particles s=1 %i: ",totPs1);

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);

		// Move particles
		// oRayTrace(pop, obj, deltaRho); <- do we need this still???

		puMove(pop); //puMove(pop, obj); Do not change functions such that PINC does
        // not work in other run modes!

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
		//gNeutralizeGrid(phi, mpiInfo);
		//gBnd(phi, mpiInfo);
        // Second run with solver to account for charges
		oSweepBiasSin( obj, n );
		oApplyCapacitanceMatrix(rho, phi, obj, mpiInfo, units); // for capMatrix - objects

		//gBnd(phi, mpiInfo);
		solve(solver, rho, phi, mpiInfo);
		//gNeutralizeGrid(phi, mpiInfo);
		//gBnd(phi, mpiInfo);
		//gHaloOp(setSlice, phi, mpiInfo, TOHALO); // Needed by sSolve but not mgSolve

		double rhoSum = gSumTruegrid(rho);
		double rhoObjSum = gSumTruegrid(rhoObj);
		MPI_Allreduce(MPI_IN_PLACE, &rhoObjSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &rhoSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		msg(STATUS,"total charge = %f PINC values",(rhoSum));
		msg(STATUS,"Object charge = %.6e C",(units->charge*rhoObjSum));
		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);
		gBnd(E, mpiInfo);
		//gBndE(E, mpiInfo); // always neumann cond

		//gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		//gZero(E);
		//gAddTo(Ext); //needs grid definition of Eext
		puAddEext(ini, pop, E); // adds same value to whole grid

		// Accelerate particle and compute kinetic energy for step n
		//acc(pop, E);
		acc(pop, E, T, S);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		if(n%10 == 0 || (n>19500 && n%10==0)){//n>122700){//50614
		//Write h5 files
			//gWriteH5(E, mpiInfo, (double) n);
			gWriteH5(rho, mpiInfo, (double) n);
			gWriteH5(rho_e, mpiInfo, (double) n);
			gWriteH5(rho_i, mpiInfo, (double) n);

			gWriteH5(phi, mpiInfo, (double) n);
			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
			//gWriteH5(rhoObj, mpiInfo, (double) n);
		}
		// if(n%1 == 0){
		// 	pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		// }

		pWriteEnergy(history,pop,(double)n,units);
        xyWrite(history,"/current/electrons/dataset",(double)n,units->current*obj->objectCurrent[0],MPI_SUM);
        xyWrite(history,"/current/ions/dataset",(double)n,units->current*obj->objectCurrent[1],MPI_SUM);
        xyWrite(history,"/potential/dataset",(double)n,units->potential*(*obj->bias),MPI_MAX);
	}

	//if(mpiInfo->mpiRank==0) {
    tMsg(t->total, "Time spent: ");
    //}

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
    gCloseH5(rhoObj);       // for capMatrix - objects
    oCloseH5(obj);          // for capMatrix - objects
    // 11.10.19 segfault seems to link to oClose(), as calling this
    // alters the segfault.

	xyCloseH5(history);

    // Free memory
    // sFree(solver);
    // mgFreeSolver(solver);
    solverFree(solver);
    gFree(rho);
    gFree(rho_e);
    gFree(rho_i);
    gFree(phi);
    free(S);
    free(T);

    gFree(E);
    pFree(pop);
    uFree(units);
    gFree(rhoObj);          // for capMatrix - objects
    oFree(obj);             // for capMatrix - objects

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);


}