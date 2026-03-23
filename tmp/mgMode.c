void mgMode(dictionary *ini){

	Units *units=uAlloc(ini);
	uNormalize(ini, units);
	//Mpi
	MpiInfo *mpiInfo = gAllocMpi(ini);

	//Rand Seed
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	//Grids
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);
	Grid *res= gAlloc(ini, SCALAR,mpiInfo);
	Grid *sol = gAlloc(ini, SCALAR,mpiInfo);
	Grid *E = gAlloc(ini, VECTOR,mpiInfo);
	Grid *error =gAlloc(ini, SCALAR,mpiInfo);

	//Multilevel grids
	Multigrid *mgPhi = mgAlloc(ini, phi);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);

	funPtr mgAlgo = getMgAlgo(ini);

	msg(STATUS, "\nMultigrid settings: \n nLevels = %d \n nPreSmooth = %d \n nCoarseSolve = %d \n nPostSmooth = %d",
	mgRho->nLevels, mgRho->nPreSmooth, mgRho->nCoarseSolve, mgRho->nPostSmooth);
	msg(STATUS, "mgLevels = %d", mgRho->nLevels);

	int rank = rho->rank;
	Timer *t = tAlloc(rank);

	// double tol = 77000.;
	// double err = tol+1.;

	//Compute stuff
	// gFillHeavi(rho, 1, mpiInfo);
	// gFillHeaviSol(sol, 1, mpiInfo);
	// gFillPoint(rho, mpiInfo);
	// gFillPolynomial(rho, mpiInfo);
	// gFillPointSol(sol, mpiInfo);
	// gFillExp(sol, mpiInfo);
	gFillSin(rho, 3, mpiInfo, 0);
	gFillSinSol(sol, 3,  mpiInfo);
	// gFillCst(rho, mpiInfo);
	// gFillRng(rho, mpiInfo, rng);

	// gHaloOp(setSlice, sol, mpiInfo);
	// gFinDiff2nd3D(rho, sol);

	gNeutralizeGrid(rho, mpiInfo);
	//
	double tol = 0.01;
	double avgError = 1;
	double errSquared;
	double resSquared;
	int run = 1;
	//
	while(avgError>tol){
		// Run solver
		tStart(t);
		mgSolveRaw(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		tStop(t);

		//Compute error
		mgCompError(phi, sol, error);
		// avgError = mgAvgError(error, mpiInfo);
		errSquared = mgSumTrueSquared(error);
		avgError = errSquared/gTotTruesize(error, mpiInfo);

		if(!(run%10))	msg(STATUS, "Avg e^2 = %.2e", errSquared);
		run++;
	}

	resSquared = mgSumTrueSquared(res);
	msg(STATUS, "Avg e^2 = %f", avgError);
	msg(STATUS, "Residual squared (res^2) = %f", resSquared);
	msg(STATUS, "Number of Cycles: %d", run);
	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*********************************************************************
	*			STORE GRIDS
	********************************************************************/
	int runNumber = iniGetInt(ini, "multigrid:runNumber");

	if(runNumber == 0){

		//(Re)Compute E, error and residual
		gHaloOp(setSlice, phi, mpiInfo, 0);
		gNeutralizeGrid(phi, mpiInfo);
		gBnd(phi, mpiInfo);
		gFinDiff1st(phi, E);
		mgCompError(phi,sol,error);
		mgResidual(res,rho, phi, mpiInfo);

		gOpenH5(ini, E, mpiInfo, units, 1.0, "E_0");
		gWriteH5(E, mpiInfo, 0.);
		gCloseH5(E);

		gOpenH5(ini, sol, mpiInfo, units, 1.0, "sol_0");
		gWriteH5(sol, mpiInfo, 0.);
		gCloseH5(sol);

		gOpenH5(ini, error, mpiInfo, units, 1.0, "error_0");
		gWriteH5(error, mpiInfo, 0.);
		gCloseH5(error);

		//Saving lvl of grids
		char fName[64];
		for(int lvl = 0; lvl <mgRho->nLevels; lvl ++){

			rho = mgRho->grids[lvl];
			phi = mgPhi->grids[lvl];
			res = mgRes->grids[lvl];

			sprintf(fName, "rho_%d", lvl);
			gOpenH5(ini, rho, mpiInfo, units, 1.0, fName);
			sprintf(fName, "phi_%d", lvl);
			gOpenH5(ini,  phi, mpiInfo, units, 1.0, fName);
			sprintf(fName, "res_%d", lvl);
			gOpenH5(ini, res, mpiInfo, units, 1.0, fName);


			gWriteH5(rho,mpiInfo,0.);
			gWriteH5(phi,mpiInfo,0.);
			gWriteH5(res,mpiInfo,0.);

			gCloseH5(phi);
			gCloseH5(rho);
			gCloseH5(res);
		}
	}

	/**************************************************************
	*		Store time spent
	*************************************************************/

	hid_t timer = xyOpenH5(ini,"timer");
	if(runNumber == 0)	xyCreateDataset(timer,"time");
	if(runNumber == 0)	xyCreateDataset(timer, "cycles");
	//  xyCreateDataset(timer,"time");
	xyWrite(timer,"time",(double) runNumber,(double) t->total,MPI_MAX);
	xyWrite(timer,"cycles",(double) runNumber,(double) run,MPI_MAX);
	xyCloseH5(timer);

	tFree(t);
	gFreeMpi(mpiInfo);

	gsl_rng_free(rng);

	uFree(units);
}