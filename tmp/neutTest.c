static void neutTest(dictionary *ini){

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set
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

	MpiInfo *mpiInfoNeut = gAllocMpi(ini);

	// For SPH neutral particles
	NeutralPopulation *neutralPop = pNeutralAlloc(ini,mpiInfoNeut);
	Grid *V   = gAlloc(ini, VECTOR,mpiInfoNeut);
	Grid *P   = gAlloc(ini, SCALAR,mpiInfoNeut);
	Grid *dKE   = gAlloc(ini, SCALAR,mpiInfoNeut);
	Grid *IE   = gAlloc(ini, SCALAR,mpiInfoNeut);
	Grid *Vtilde   = gAlloc(ini, VECTOR,mpiInfoNeut);
	Grid *Itilde   = gAlloc(ini, SCALAR,mpiInfoNeut);
	Grid *rhoNeutral = gAlloc(ini, SCALAR,mpiInfoNeut);
	Grid *rhoObj = gAlloc(ini, SCALAR,mpiInfoNeut); // for capMatrix - objects

	gZero(rhoNeutral);
   	gZero(P);
	gZero(dKE);
	gZero(IE);
   	gZero(V);
	gZero(Itilde);
   	gZero(Vtilde);

	PincObject *obj =objoAlloc(ini,mpiInfoNeut,units); // for capMatrix - objects

    // for SPH neutrals
	gCreateNeighborhood(ini, mpiInfoNeut, rhoNeutral);
    // We assume same form on neutral density grid and charged density grid

    // need for SPH neutrals a function

	neSetBndSlices( IE, mpiInfoNeut);
	neSetBndSlicesVel(ini, V, mpiInfoNeut);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfoNeut->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */

    // SPH neutrals
    gOpenH5(ini, rhoNeutral, mpiInfoNeut, units, 1, "rhoNeutral");
    gOpenH5(ini, P,   mpiInfoNeut, units, 1, "P");
	gOpenH5(ini, IE,   mpiInfoNeut, units, 1, "IE");
	gOpenH5(ini, V,   mpiInfoNeut, units, units->velocity, "V");

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	 * INITIAL CONDITIONS
	 */

	// SPH neutrals
	nePosUniform(ini, neutralPop, mpiInfoNeut, rngSync);
	neVelDrift(ini, neutralPop);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// SPH neutrals
	nePurgeGhost(neutralPop, rhoNeutral);
	neFillGhost(ini,neutralPop,rngSync,mpiInfoNeut);

	/*
	 * INITIALIZATION (E.g. half-step)
	 */

	nuObjectpurge(neutralPop,rhoObj,obj);

    NeutralDistr3D1(neutralPop, rhoNeutral);
	gHaloOp(addSlice, rhoNeutral, mpiInfoNeut, FROMHALO);
	gHaloOp(setSlice, rhoNeutral, mpiInfoNeut, TOHALO);
	NeutralDistr3D1Vector(neutralPop,V,rhoNeutral);
	gHaloOp(addSlice, V, mpiInfoNeut, FROMHALO);
	nuGBndVel(V,mpiInfoNeut);
	gHaloOp(setSlice, V, mpiInfoNeut, TOHALO);

	neSetI(IE,V,rhoNeutral,ini);
	neSetBndSlicesEnerg(ini,IE,rhoNeutral,mpiInfoNeut);
	gHaloOp(setSlice, IE, mpiInfoNeut, TOHALO);
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",3);
	double multiplyIEBy = 4.;
	int sliceDim = 0;
	neMultiplySlice(IE,(int)(trueSize[0]/2)-1,sliceDim,multiplyIEBy, neutralPop);
	neMultiplySlice(IE,(int)(trueSize[0]/2),sliceDim,multiplyIEBy, neutralPop);
	neMultiplySlice(IE,(int)(trueSize[0]/2)+1,sliceDim,multiplyIEBy, neutralPop);

	// SPH neutrals
	neExtractEmigrants3DOpen(neutralPop, mpiInfoNeut);
	neMigrate(neutralPop, mpiInfoNeut, rhoNeutral);

	gWriteH5(rhoNeutral, mpiInfoNeut, (double) 0);
	gWriteH5(IE, mpiInfoNeut, (double) 0);
	gWriteH5(P, mpiInfoNeut, (double) 0);
	gWriteH5(V, mpiInfoNeut, (double) 0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfoNeut->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS," Computing time-step %i",n);
		msg(STATUS, "Nr. of particles %i: ",(neutralPop->iStop[0]- neutralPop->iStart[0]));
		double gridEnerg = gSumTruegrid(IE);
		double Vsum = gSumTruegrid(V);
		double rhosum = gSumTruegrid(rhoNeutral);
		msg(STATUS,"grid energy = %f",gridEnerg);
		msg(STATUS,"Vsum = %f",Vsum);
		msg(STATUS,"rhosum = %f \n",rhosum);

        neVelAssertMax(neutralPop,maxVel);

		tStart(t);

		nePressureSolve3D(P,IE,rhoNeutral,neutralPop);
		nuObjectSetVal(P,0.,obj);
		neApplyObjI(obj, P );
		gHaloOp(setSlice, P, mpiInfoNeut, TOHALO);

		neAdvectV(V,Vtilde,P,rhoNeutral,neutralPop);
		gHaloOp(setSlice, Vtilde, mpiInfoNeut, TOHALO);

		neAdvectI(IE,Itilde,P,V,rhoNeutral,neutralPop);
		gHaloOp(setSlice, Itilde, mpiInfoNeut, TOHALO);
		neMove(neutralPop,V);

		neExtractEmigrants3DOpen(neutralPop, mpiInfoNeut);
		neMigrate(neutralPop, mpiInfoNeut, rhoNeutral);


		neConvectKE(dKE,Vtilde,rhoNeutral, neutralPop);
		gHaloOp(setSlice, dKE, mpiInfoNeut, TOHALO);
		neConvectV(V,Vtilde,rhoNeutral,neutralPop );
		gHaloOp(setSlice, V, mpiInfoNeut, TOHALO);

		nuGBndVel(V,mpiInfoNeut);

		nePurgeGhost(neutralPop, rhoNeutral);
		neFillGhost(ini,neutralPop,rngSync,mpiInfoNeut);
		NeutralDistr3D1(neutralPop, rhoNeutral);
		gHaloOp(addSlice, rhoNeutral, mpiInfoNeut, FROMHALO);
		neConvectI(IE,Itilde,dKE,rhoNeutral,neutralPop );

		neSetBndSlicesEnerg(ini,IE,rhoNeutral,mpiInfoNeut);
		nuGBnd(IE,mpiInfoNeut);

		gHaloOp(setSlice, IE, mpiInfoNeut, TOHALO);
		nuObjectSetVal(IE,0.,obj);
		neApplyObjI(obj, IE );

		neApplyObjVel(obj,V);

		if(n%10 == 0 || n>4900){
			gWriteH5(V, mpiInfoNeut, (double) n);
			gWriteH5(rhoNeutral, mpiInfoNeut, (double) n);
			gWriteH5(P, mpiInfoNeut, (double) n);
			gWriteH5(IE, mpiInfoNeut, (double) n);
		}
	}

    tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */

	gFreeMpi(mpiInfoNeut);

	// Close h5 files

    oCloseH5(obj);          // for capMatrix - objects

    // SPH neutrals
	gCloseH5(rhoNeutral);
    gCloseH5(P);
	gCloseH5(IE);
	gCloseH5(V);

	gFree(rhoNeutral);
	gFree(P);
	gFree(IE);
	gFree(V);
	gFree(Itilde);
	gFree(Vtilde);

	pNeutralFree(neutralPop);

	uFree(units);
    gFree(rhoObj);          // for capMatrix - objects
    oFree(obj);             // for capMatrix - objects

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);
}