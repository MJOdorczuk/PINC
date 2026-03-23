void mgModeErrorScaling(dictionary *ini){

	Units *units=uAlloc(ini);
	uNormalize(ini, units);

	//Mpi
	MpiInfo *mpiInfo = gAllocMpi(ini);

	//Grids
	Grid *phi 	= gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho 	= gAlloc(ini, SCALAR,mpiInfo);
	Grid *res 	= gAlloc(ini, SCALAR,mpiInfo);
	Grid *E		= gAlloc(ini, VECTOR,mpiInfo);

	//Multigrids
	Multigrid *mgPhi = mgAlloc(ini, phi);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);

	//Error and solutions
	Grid *error = gAlloc(ini, SCALAR,mpiInfo);
	Grid *errorE= gAlloc(ini, VECTOR,mpiInfo);
	Grid *sol 	= gAlloc(ini, SCALAR,mpiInfo);
	Grid *solE	= gAlloc(ini, VECTOR,mpiInfo);

	//Compute stuff
	// gFillHeavi(rho, 1, mpiInfo);

	// gFillHeaviSol(sol, 1, mpiInfo);
	gFillSin(rho, 1, mpiInfo, 0);
	gFillSinSol(sol, 1, mpiInfo);
	gFillSinESol(solE, 1, mpiInfo);

	if(mpiInfo->mpiRank==0)	aiPrint(&rho->trueSize[1], rho->rank-1);

	funPtr mgAlgo = getMgAlgo(ini);

	//Solve
	mgSolveRaw(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

	//(Re)Compute E, error and residual
	gHaloOp(setSlice, phi, mpiInfo, TOHALO);
	gNeutralizeGrid(phi, mpiInfo);
	gBnd(phi, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	mgCompError(phi,sol,error);
	mgCompError(E, solE, errorE);
	mgResidual(res,rho, phi, mpiInfo);

	/*********************************************************************
	*			STORE GRIDS
	********************************************************************/
	int runNumber 	= iniGetInt(ini, "multigrid:runNumber");

	char *fName 	= malloc(8*sizeof(fName));

	sprintf(fName, "rho_%d", runNumber);
	gOpenH5(ini, rho, mpiInfo, units, 1.0, fName);
	gWriteH5(rho, mpiInfo, 0.);
	gCloseH5(rho);

	sprintf(fName, "phi_%d", runNumber);
	gOpenH5(ini, phi, mpiInfo, units, 1.0, fName);
	gWriteH5(phi, mpiInfo, 0.);
	gCloseH5(phi);

	sprintf(fName, "res_%d", runNumber);
	gOpenH5(ini, res, mpiInfo, units, 1.0, fName);
	gWriteH5(res, mpiInfo, 0.);
	gCloseH5(res);

	sprintf(fName, "E_%d", runNumber);
	gOpenH5(ini, E, mpiInfo, units, 1.0, fName);
	gWriteH5(E, mpiInfo, 0.);
	gCloseH5(E);

	sprintf(fName, "sol_%d", runNumber);
	gOpenH5(ini, sol, mpiInfo, units, 1.0, fName);
	gWriteH5(sol, mpiInfo, 0.);
	gCloseH5(sol);

	sprintf(fName, "error_%d", runNumber);
	gOpenH5(ini, error, mpiInfo, units, 1.0, fName);
	gWriteH5(error, mpiInfo, 0.);
	gCloseH5(error);

	sprintf(fName, "solE_%d", runNumber);
	gOpenH5(ini, solE, mpiInfo, units, 1.0, fName);
	gWriteH5(solE, mpiInfo, 0.);
	gCloseH5(solE);

	sprintf(fName, "errorE_%d", runNumber);
	gOpenH5(ini, errorE, mpiInfo, units, 1.0, fName);
	gWriteH5(errorE, mpiInfo, 0.);
	gCloseH5(errorE);

	//Freedom
	gFreeMpi(mpiInfo);
	free(fName);

	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	gFree(error);
	gFree(sol);

	uFree(units);
}