void sMode(dictionary *ini){


	MpiInfo *mpiInfo = gAllocMpi(ini);
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);

	SpectralSolver *solver = sAlloc(ini, rho, phi);

	int *trueSize = rho->trueSize;
	int *nGhostLayers = rho->nGhostLayers;
	double *rhoValStart = &rho->val[nGhostLayers[1]];
	double *phiValStart = &phi->val[nGhostLayers[1]];
	for(long int j=0; j<trueSize[1]; j++){
		rhoValStart[j] = sin(2*M_PI*j/trueSize[1]);
	}

	adPrint(rhoValStart, trueSize[1]);
	sSolve(solver, rho, phi);
	adPrint(phiValStart, trueSize[1]);

	sFree(solver);
	gFree(rho);
	gFree(phi);
	gFreeMpi(mpiInfo);

}