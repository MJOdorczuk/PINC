/**
 * @file		neutrals.c
 * @brief		PIC simulation for neutrals.
 * @author		Steffen Brask <steffen.brask@fys.uio.no>
 *
 * TBD
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "core.h"
#include "pusher.h"
#include "neutrals.h"

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/
void neCutParticle(NeutralPopulation *pop, int s, long int p, double *pos, double *vel)
{

	int nDims = pop->nDims;
	long int pLast = (pop->iStop[s] - 1) * nDims;
	for (int d = 0; d < nDims; d++)
	{
		pos[d] = pop->pos[p + d];
		vel[d] = pop->vel[p + d];
		pop->pos[p + d] = pop->pos[pLast + d];
		pop->vel[p + d] = pop->vel[pLast + d];
	}

	pop->iStop[s]--;
}

void nePToLocalFrame(NeutralPopulation *pop, const MpiInfo *mpiInfo)
{

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{
		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (long int i = iStart; i < iStop; i++)
		{
			double *pos = &pop->pos[i * nDims];
			for (int d = 0; d < nDims; d++)
				pos[d] -= offset[d];
		}
	}
}

void nePToGlobalFrame(NeutralPopulation *pop, const MpiInfo *mpiInfo)
{

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{
		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (long int i = iStart; i < iStop; i++)
		{

			double *pos = &pop->pos[i * nDims];
			for (int d = 0; d < nDims; d++)
				pos[d] += offset[d];
		}
	}
}

void neGNeumann(Grid *grid, const int boundary)
{
	// Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice; // two slices in each dim

	// Compute dimensions and slicesize
	int d = boundary % rank;
	int offset = 1 + (boundary > rank) * (size[d] - 3);

	// Number of elements in slice
	long int nSliceMax = 0; // TODO: this can/should be stored.
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}

	getSlice(&bndSlice[boundary * nSliceMax], grid, d, offset); // edge before halo
	setSlice(&bndSlice[boundary * nSliceMax], grid, d, offset - 1 + (boundary > rank) * 2); // halo

	return;
}

void gDirichletVel(Grid *grid, const int boundary)
{
	// Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice;

	// Compute dimensions and size of slice
	int d = boundary % rank;
	int offset = 1 + (boundary > rank) * (size[d] - 3);

	// Number of elements in slice
	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}
	nSliceMax = nSliceMax * (rank - 1);
	setSlice(&bndSlice[boundary * nSliceMax], grid, d, offset); // edge before halo
	setSlice(&bndSlice[boundary * nSliceMax], grid, d, offset - 1 + (boundary > rank) * 2); // halo

	return;
}

void gDirichletEnerg(Grid *grid, const int boundary)
{
	// Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice;

	// Compute dimensions and size of slice
	int d = boundary % rank;
	int offset = 1 + (boundary > rank) * (size[d] - 3);

	// Number of elements in slice
	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}

	setSlice(&bndSlice[boundary * nSliceMax], grid, d, offset); // edge before halo
	setSlice(&bndSlice[boundary * nSliceMax], grid, d, offset - 1 + (boundary > rank) * 2); // halo

	return;
}

void nuGBnd(Grid *grid, const MpiInfo *mpiInfo)
{
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	// If periodic neutralize phi

	bool periodic = true;
	for (int d = 1; d < rank; d++)
	{
		if (bnd[d] != PERIODIC)
		{
			periodic = false;
		}
	}
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (bnd[d] != PERIODIC)
		{
			periodic = false;
		}
	}

	// Lower edge
	for (int d = 1; d < rank; d++)
	{
		if (subdomain[d - 1] == 0)
		{
			if (bnd[d] == DIRICHLET)
			{
				gDirichletEnerg(grid, d);
			}
			else if (bnd[d] == NEUMANN)
			{
				neGNeumann(grid, d);
			}
		}
	}

	// Higher edge
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (subdomain[d - rank - 1] == nSubdomains[d - rank - 1] - 1)
		{
			if (bnd[d] == DIRICHLET)
				gDirichletEnerg(grid, d);
			if (bnd[d] == NEUMANN)
				neGNeumann(grid, d);
		}
	}

	return;
}

void nuGBndVel(Grid *grid, const MpiInfo *mpiInfo)
{

	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	// If periodic neutralize phi

	bool periodic = true;
	for (int d = 1; d < rank; d++)
	{
		if (bnd[d] != PERIODIC)
		{
			periodic = false;
		}
	}
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (bnd[d] != PERIODIC)
		{
			periodic = false;
		}
	}
	if (periodic == true)

	// Lower edge
	for (int d = 1; d < rank; d++)
	{
		if (subdomain[d - 1] == 0)
		{
			if (bnd[d] == DIRICHLET)
			{
				gDirichletVel(grid, d);
			}
			else if (bnd[d] == NEUMANN)
			{
				neGNeumann(grid, d);
			}
		}
	}

	// Higher edge
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (subdomain[d - rank - 1] == nSubdomains[d - rank - 1] - 1)
		{
			if (bnd[d] == DIRICHLET)
			{
				gDirichletVel(grid, d);
			}
			if (bnd[d] == NEUMANN)
				neGNeumann(grid, d);
		}
	}

	return;
}

void neNewparticle(NeutralPopulation *pop, int s, const double *pos, const double *vel)
{

	int nDims = pop->nDims;
	long int *iStart = pop->iStart;
	long int *iStop = pop->iStop; // New particle added here

	if (iStop[s] >= iStart[s + 1])
		msg(WARNING, "Not enough allocated memory to add new particle to specie"
					 "%i. New particle ignored.",
			s);
	else
	{

		long int p = iStop[s] * nDims;
		for (int d = 0; d < nDims; d++)
		{
			pop->pos[p + d] = pos[d];
			pop->vel[p + d] = vel[d];
		}

		iStop[s]++;
	}
}

void neInjectParticles(int slicePos, int dim, int multiplyDens, const dictionary *ini,
	NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims); // TODO: is it used at all?
	double *velDrift = iniGetDoubleArr(ini, "collisions:neutralDrift", nDims * nSpecies);
	double *velThermal = iniGetDoubleArr(ini, "collisions:thermalVelocityNeutrals", nSpecies);
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	int *nGhostLayers = iniGetIntArr(ini, "grid:nGhostLayers", 2 * nDims);
	int *L = gGetGlobalSize(ini);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	double pos[nDims];
	double vel[nDims];

	nePToGlobalFrame(pop, mpiInfo);

	for (int s = 0; s < nSpecies; s++)
	{
		double velTh = velThermal[s];

		// TODO: why inside loop?
		int d = dim;
		long int globalSizeProd = (L[0] * L[1] * L[2]); // TODO: make nDimensional
		long int sliceSize = (globalSizeProd / (L[d]));

		// slice * particles per cell
		long int newParticles = multiplyDens * ((sliceSize * nParticles[s]) / globalSizeProd);

		for (long int i = 0; i < newParticles; i++)
		{
			// generate velocity for particle
			for (int dd = 0; dd < nDims; dd++)
			{
				vel[dd] = velDrift[(s * nDims) + dd] + gsl_ran_gaussian_ziggurat(rng, velTh);
			}

			// Generate position for particle
			for (int dd = 0; dd < nDims; dd++)
			{
				pos[dd] = (L[dd]) * gsl_rng_uniform_pos(rng) - nGhostLayers[d];
			}
			pos[d] = slicePos + (gsl_rng_uniform_pos(rng)); // in lower ghost
			int correctRange = 0;
			for (int dd = 0; dd < nDims; dd++)
			{
				correctRange += (subdomain[dd] == (int)(posToSubdomain[dd] * pos[dd]));
			}
			// Add only if particle resides in this sub-domain.
			if (correctRange == nDims)
			{
				neNewparticle(pop, s, pos, vel);
			}
		}
	}
	nePToLocalFrame(pop, mpiInfo);

	free(velDrift);
	free(velThermal);
	free(trueSize);
	free(nParticles);
	free(nGhostLayers);
	free(L);

	return;
}

void neMultiplySlice(Grid *target, int slicePos, int dim, double multiplyBy, NeutralPopulation *pop)
{
	int nDims = pop->nDims;

	int rank = nDims + 1;
	int *size = target->size;

	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}

	// TODO: why is it 'addSlice'?
	double *addSlice = malloc(nSliceMax * sizeof(addSlice));

	getSlice(addSlice, target, dim, slicePos);

	// TODO: do we not have helper functions for this?
	for (long int i = 0; i < nSliceMax; i++)
	{
		addSlice[i] = addSlice[i] * multiplyBy;
	}
	setSlice(addSlice, target, dim, slicePos);

	free(addSlice);

	return;
}

void neScatterParticle(NeutralPopulation *pop, double *pos, double *vel)
{

	// neutral object coll
	// TODO: Do these have to be arrays?
	int nDims = pop->nDims;
	double oldPos[nDims];
	double oldVel[nDims];
	int traversed[nDims];
	double delta[nDims];

	for (int d = 0; d < nDims; d++)
	{
		oldVel[d] = vel[d]; // using oldVel temporarily.
		oldPos[d] = pos[d] - oldVel[d];
		if ((int)(pos[d]) - (int)(oldPos[d]) != 0)
		{
			traversed[d] = (int)(pos[d]) - (int)(oldPos[d]);
			delta[d] = pos[d] - (int)pos[d];
			pos[d] = ((int)pos[d] - delta[d]);
			vel[d] = -oldVel[d];
			if (traversed[d] == -1)
			{ // hotfix
				pos[d] += 2;
			}
		}
	}
}

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

NeutralPopulation *pNeutralAlloc(const dictionary *ini, const MpiInfo *mpiInfo)
{

	// Get MPI info
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int mpiRank = mpiInfo->mpiRank;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

	// Load data
	int nSpecies = iniGetInt(ini, "collisions:nSpeciesNeutral");
	int nDims = iniGetInt(ini, "grid:nDims");
	char **boundaries = iniGetStrArr(ini, "grid:boundaries", 2 * nDims);

	// Number of particles to allocate for (for all computing nodes)
	long int *nAllocTotal = iniGetLongIntArr(ini, "population:nAlloc", nSpecies);

	// Determine memory to allocate for this node
	long int *nAlloc = malloc(nSpecies * sizeof(long int));
	for (int s = 0; s < nSpecies; s++)
	{
		nAlloc[s] = ceil((double)nAllocTotal[s] / size);
		if (nAlloc[s] * size != nAllocTotal[s])
			msg(WARNING, "increased number of allocated particles from %i to %i"
						 "to get integer per computing node",
				nAllocTotal[s], nAlloc[s] * size);
	}

	long int *iStart = malloc((nSpecies + 1) * sizeof(long int));
	long int *iStop = malloc(nSpecies * sizeof(long int));

	iStart[0] = 0;
	for (int s = 1; s < nSpecies + 1; s++)
		iStart[s] = iStart[s - 1] + nAlloc[s - 1];
	for (int s = 0; s < nSpecies; s++)
		iStop[s] = iStart[s]; // No particles yet

	int rank = nDims + 1;
	bndType *bnd = malloc(2 * rank * sizeof(*bnd));

	bnd[0] = NONE;	  // initialize
	bnd[rank] = NONE; // initialize

	int b = 0;
	// TODO: Why not just dd from zero?
	for (int d = 1; d < rank; d++)
	{
		int dd = d - 1;
		int firstElem = mpiRank - subdomain[dd] * nSubdomainsProd[dd];

		int upperSubdomain = firstElem + ((subdomain[dd] + 1) % nSubdomains[dd]) * nSubdomainsProd[dd];
		int lowerSubdomain = firstElem + ((subdomain[dd] - 1 + nSubdomains[dd]) % nSubdomains[dd]) * nSubdomainsProd[dd];

		// lower
		int r = dd + 1;
		if (lowerSubdomain >= mpiRank)
		{
			if (!strcmp(boundaries[b], "PERIODIC"))
				bnd[r] = PERIODIC;
			else if (!strcmp(boundaries[b], "DIRICHLET"))
				bnd[r] = DIRICHLET;
			else if (!strcmp(boundaries[b], "NEUMANN"))
				bnd[r] = NEUMANN;
			else
				msg(ERROR, "%s invalid value for grid:boundaries", boundaries[b]);
		}
		else if (lowerSubdomain < mpiRank)
		{
			bnd[r] = PERIODIC;
		}
		else
		{
			bnd[r] = NONE; // initialize
		}
		// upper
		r = rank + dd + 1;
		if (upperSubdomain <= mpiRank)
		{
			if (!strcmp(boundaries[b + rank - 1], "PERIODIC"))
				bnd[r] = PERIODIC;
			else if (!strcmp(boundaries[b + rank - 1], "DIRICHLET"))
				bnd[r] = DIRICHLET;
			else if (!strcmp(boundaries[b + rank - 1], "NEUMANN"))
				bnd[r] = NEUMANN;
			else
				msg(ERROR, "%s invalid value for grid:boundaries", boundaries[b]);
		}
		else if (upperSubdomain > mpiRank)
		{
			bnd[r] = PERIODIC;
		}
		else
		{
			bnd[r] = NONE; // initialize
		}
		b++;
	}

	double *mass = iniGetDoubleArr(ini, "collisions:neutralMass", nSpecies);
	double *density = iniGetDoubleArr(ini, "collisions:numberDensityNeutrals", nSpecies);
	double rho0 = density[0]; //*mass[0];

	double NvelThermal = iniGetDouble(ini, "collisions:thermalVelocityNeutrals");
	double stiffnessC = (NvelThermal); //((7.)*(NvelThermal*NvelThermal))/(mass[0]*rho0); // Monaghan (1994)
	NeutralPopulation *neutralPop = malloc(sizeof(*neutralPop));
	neutralPop->pos = malloc((long int)nDims * iStart[nSpecies] * sizeof(double));
	neutralPop->vel = malloc((long int)nDims * iStart[nSpecies] * sizeof(double));
	neutralPop->nSpeciesNeutral = nSpecies;
	neutralPop->nDims = nDims;
	neutralPop->iStart = iStart;
	neutralPop->iStop = iStop;
	neutralPop->mass = mass; // iniGetDoubleArr(ini,"collisions:neutralMass",nSpecies);
	neutralPop->bnd = bnd;
	neutralPop->stiffnessConstant = stiffnessC; // 0.00001
	neutralPop->rho0 = rho0;

	free(nAlloc);
	free(nAllocTotal);
	free(density);
	free(boundaries);
	return neutralPop;
}

void pNeutralFree(NeutralPopulation *pop)
{

	free(pop->pos);
	free(pop->vel);
	free(pop->iStart);
	free(pop->iStop);
	free(pop->mass);
	free(pop->bnd);
	free(pop);
}

void nePosLattice(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo)
{

	// Read from ini
	int nDims = pop->nDims;
	int nSpecies = pop->nSpeciesNeutral;
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = gGetGlobalSize(ini);
	long int V = gGetGlobalVolume(ini);

	for (int s = 0; s < nSpecies; s++)
	{

		// Particle-particle distance in lattice
		double l = pow(V / (double)nParticles[s], 1.0 / nDims);

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart;
		double *pos = &pop->pos[iStart * nDims];

		// Iterate through all particles to be generated
		// Generate particles on global frame on all nodes and discard the ones
		// out of range. This is simpler as it resembles pPosUniform()
		for (long int i = 0; i < nParticles[s]; i++)
		{

			double linearPos = l * i;
			for (int d = 0; d < nDims; d++)
			{
				pos[d] = fmod(linearPos, L[d]);
				// TODO: Is this really needed?
				linearPos /= L[d];
			}

			// Count the number of dimensions where the particle resides in
			// the range of this node
			int correctRange = 0;
			for (int d = 0; d < nDims; d++)
				correctRange += (subdomain[d] == (int)(posToSubdomain[d] * pos[d]));

			// Iterate only if particle resides in this sub-domain.
			if (correctRange == nDims)
			{
				pos += nDims;
				iStop++;
			}
		}

		if (iStop > pop->iStart[s + 1])
		{
			int allocated = pop->iStart[s + 1] - iStart;
			int generated = iStop - iStart;
			msg(ERROR, "allocated only %i particles of specie %i per node but"
					   "%i generated",
				allocated, s, generated);
		}

		pop->iStop[s] = iStop;
	}

	nePToLocalFrame(pop, mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);
}

void nePosUniform(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo, const gsl_rng *rng)
{

	// Read from ini
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = gGetGlobalSize(ini);

	for (int s = 0; s < nSpecies; s++)
	{

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart;
		double *pos = &pop->pos[iStart * nDims];

		// Iterate through all particles to be generated. Same seed on all MPI
		// nodes ensure same particles are generated everywhere.
		for (long int i = 0; i < nParticles[s]; i++)
		{

			// Generate position for particle i
			for (int d = 0; d < nDims; d++)
			{
				pos[d] = L[d] * gsl_rng_uniform_pos(rng);
			}

			// Count the number of dimensions where the particle resides in
			// the range of this node
			int correctRange = 0;
			for (int d = 0; d < nDims; d++)
				correctRange += (subdomain[d] == (int)(posToSubdomain[d] * pos[d]));

			// Iterate only if particle resides in this sub-domain.
			if (correctRange == nDims)
			{
				pos += nDims;
				iStop++;
			}
		}

		if (iStop > pop->iStart[s + 1])
		{
			int allocated = pop->iStart[s + 1] - iStart;
			int generated = iStop - iStart;
			msg(ERROR, "allocated only %i particles of specie %i per node but"
					   "%i generated",
				allocated, s, generated);
		}

		pop->iStop[s] = iStop;
	}

	nePToLocalFrame(pop, mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);
}

void neVelMaxwell(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	double *velDrift = iniGetDoubleArr(ini, "collisions:neutralDrift", nDims * nSpecies);
	double *velThermal = iniGetDoubleArr(ini, "collisions:thermalVelocityNeutrals", nSpecies);
	long int index = 0;

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = velThermal[s];
		velTh = velTh;

		for (long int i = iStart; i < iStop; i++)
		{

			double *vel = &pop->vel[i * nDims];
			for (int d = 0; d < nDims; d++)
			{
				index = (s * nDims) + d;
				vel[d] = velDrift[index] + gsl_ran_gaussian_ziggurat(rng, velTh);
			}
		}
	}
	free(velDrift);
	free(velThermal);
}

void neVelDrift(const dictionary *ini, NeutralPopulation *pop)
{

	// test function. takes only two species
	int nDims = pop->nDims;
	int nSpecies = pop->nSpeciesNeutral;
	double *velDrift = iniGetDoubleArr(ini, "collisions:neutralDrift", nDims * nSpecies);

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (long int i = iStart; i < iStop; i++)
		{
			for (int d = 0; d < nDims; d++)
			{
				pop->vel[i * nDims + d] = velDrift[d];
			}
		}
	}
}

void nePurgeGhost(NeutralPopulation *pop, const Grid *grid)
{

	// this will delete all particles residing on Dirichlet boundaries
	int *size = grid->size;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	double *pos = pop->pos;
	double *vel = pop->vel;
	bool cut = false;
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];
		for (int i = iStart; i < iStop; i++)
		{
			for (int d = 0; d < nDims; d++)
			{
				if ((pos[i * nDims + d] > size[d + 1] - nGhostLayers[d + 1] - 1
				 && (bnd[d + rank + 1] == DIRICHLET || bnd[d + rank + 1] == NEUMANN)))
				{
					cut = true;
				}
				if ((pos[i * nDims + d] < nGhostLayers[d + rank + 1]
				 && (bnd[d + 1] == DIRICHLET || bnd[d + 1] == NEUMANN)))
				{
					cut = true;
				}
			}
			if (cut == true)
			{
				neCutParticle(pop, s, i * nDims, pos, vel);
				cut = false;
				iStop--;
				i--;
			}
		}
	}
}

void neFillGhost(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	bndType *bnd = pop->bnd;
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims);
	double *velDrift = iniGetDoubleArr(ini, "collisions:neutralDrift", nDims * nSpecies);
	double *velThermal = iniGetDoubleArr(ini, "collisions:thermalVelocityNeutrals", nSpecies);
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	int *nGhostLayers = iniGetIntArr(ini, "grid:nGhostLayers", 2 * nDims);
	int *L = gGetGlobalSize(ini);
	int rank = nDims + 1;

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	double pos[nDims];
	double vel[nDims];

	nePToGlobalFrame(pop, mpiInfo);

	for (int s = 0; s < nSpecies; s++)
	{
		double velTh = velThermal[s];

		for (int d = 0; d < nDims; d++)
		{
			long int globalSizeProd = (L[0] * L[1] * L[2]); // TODO: make nDimensional

			long int sliceSize = (globalSizeProd / (L[d]));
			long int newParticles = ((sliceSize * nParticles[s]) / globalSizeProd); // slice * particles per cell

			// Lower ghost slice
			if (bnd[d + 1] == DIRICHLET || bnd[d + 1] == NEUMANN)
			{
				for (long int i = 0; i < newParticles; i++)
				{
					// generate velocity for particle
					for (int dd = 0; dd < nDims; dd++)
					{
						vel[dd] = velDrift[(s * nDims) + dd] + gsl_ran_gaussian_ziggurat(rng, velTh);
					}

					// Generate position for particle
					for (int dd = 0; dd < nDims; dd++)
					{
						pos[dd] = (L[dd]) * gsl_rng_uniform_pos(rng) - nGhostLayers[d] * 0.5;
					}
					pos[d] = nGhostLayers[d] * (gsl_rng_uniform_pos(rng)) - nGhostLayers[d]; // in lower ghost
					int correctRange = 0;
					for (int dd = 0; dd < nDims; dd++)
					{
						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd] * pos[dd]));
					}
					// Add only if particle resides in this sub-domain.
					if (correctRange == nDims)
					{
						neNewparticle(pop, s, pos, vel);
					}
				}
			}

			// Upper ghost slice
			if (bnd[d + rank + 1] == DIRICHLET || bnd[d + rank + 1] == NEUMANN)
			{
				for (long int i = 0; i < newParticles; i++)
				{
					// generate velocity for particle

					for (int dd = 0; dd < nDims; dd++)
					{
						vel[dd] = velDrift[(s * nDims) + dd] + gsl_ran_gaussian_ziggurat(rng, velTh);
					}

					// Generate position for particle
					for (int dd = 0; dd < nDims; dd++)
					{
						pos[dd] = (L[dd]) * gsl_rng_uniform_pos(rng) - nGhostLayers[d] * 0.5;
					}
					pos[d] = L[d] + nGhostLayers[d] * (gsl_rng_uniform_pos(rng)) - nGhostLayers[d]; // in lower ghost

					int correctRange = 0;
					for (int dd = 0; dd < nDims; dd++)
						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd] * (pos[dd])));

					// Add only if particle resides in this sub-domain.
					if (correctRange == nDims)
					{
						neNewparticle(pop, s, pos, vel);
					}
				}
			}
		}
	}
	nePToLocalFrame(pop, mpiInfo);

	free(velDrift);
	free(velThermal);
	free(trueSize);
	free(nParticles);
	free(nGhostLayers);
	free(L);

	return;
}

// #########################################
//  Distributer
//  ########################################

// TODO: Should the name not be neSanity?
static void puSanity(dictionary *ini, const char *name, int dim, int order)
{

	int nDims = iniGetInt(ini, "grid:nDims");
	int *nGhostLayers = iniGetIntArr(ini, "grid:nGhostLayers", 2 * nDims);
	double *thresholds = iniGetDoubleArr(ini, "grid:thresholds", 2 * nDims);

	// TBD: can be improved by checking dimensions separately
	int minLayers = aiMin(nGhostLayers, 2 * nDims);
	double minThreshold = adMin(thresholds, 2 * nDims);
	double maxThreshold = adMax(thresholds, 2 * nDims);

	if (nDims != dim && dim != 0)
		msg(ERROR, "%s only supports grid:nDims=%d", name, dim);

	int reqLayers = 0;
	// TODO: Why is the 0 case needed? It is already zero
	if (order == 0)
		reqLayers = 0;
	if (order == 1)
		reqLayers = 1;
	if (order == 2)
		reqLayers = 1;

	// TODO: Should it not be < reqLayers?
	if (minLayers < 1)
		msg(ERROR, "%s requires grid:nGhostLayers >=%d", name, reqLayers);

	double reqMinThreshold = 0;
	if (order == 0)
		reqMinThreshold = -0.5;
	if (order == 1)
		reqMinThreshold = 0;
	if (order == 2)
		reqMinThreshold = 0.5;

	if (minThreshold < reqMinThreshold)
		msg(ERROR, "%s requires grid:thresholds >=%.1f", name, reqMinThreshold);

	double reqMaxThreshold = minLayers - 0.5;

	if (maxThreshold > reqMaxThreshold)
		msg(ERROR, "%s requires grid:thresholds <= grid:nGhostLayers - 0.5", name);

	if (minThreshold == reqMaxThreshold)
		msg(WARNING, "%s is not very well tested for grid:thresholds of exactly equal to grid:nGhostLayers - 0.5", name);

	free(nGhostLayers);
	free(thresholds);
}

funPtr NeutralDistr3D1_set(dictionary *ini)
{
	puSanity(ini, "puDistr3D1", 3, 1); // reuse puSanity for now
	return NeutralDistr3D1;
}
void NeutralDistr3D1(const NeutralPopulation *pop, Grid *rho)
{

	gZero(rho);
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;
	int nSpecies = pop->nSpeciesNeutral;

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (int i = iStart; i < iStop; i++)
		{

			double *pos = &pop->pos[3 * i];

			int j = (int)(pos[0] + 0.5); // for NGP
			int k = (int)(pos[1] + 0.5);
			int l = (int)(pos[2] + 0.5);

			//  Index of neighbouring nodes
			long int p = j + k * sizeProd[2] + l * sizeProd[3];
			val[p] += 1.0; // xcomp*ycomp*zcomp;
		}
	}
}

funPtr NeutralDistr3D1Vector_set(dictionary *ini)
{
	puSanity(ini, "puDistr3D1", 3, 1);
	return NeutralDistr3D1Vector;
}
void NeutralDistr3D1Vector(const NeutralPopulation *pop, Grid *bulkV, Grid *rho)
{

	gZero(bulkV);
	double *val = bulkV->val;
	long int *sizeProd = bulkV->sizeProd;
	long int *scalarSizeProd = rho->sizeProd;
	double *rhoVal = rho->val;
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{
		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (int i = iStart; i < iStop; i++)
		{

			double *pos = &pop->pos[3 * i];
			double *vel = &pop->vel[3 * i];

			// Integer parts of position
			int j = (int)(pos[0] + 0.5); // for NGP
			int k = (int)(pos[1] + 0.5);
			int l = (int)(pos[2] + 0.5);

			long int scalarp = j + k * scalarSizeProd[2] + l * scalarSizeProd[3];
			long int p = j * sizeProd[1] + k * sizeProd[2] + l * sizeProd[3];

			if (p >= sizeProd[4])
			{
				msg(ERROR, "Particle %i at (%f,%f,%f) out-of-bounds, tried to access node %li", i, pos[0], pos[1], pos[2], p);
			}
			for (int d = 0; d < nDims; d++)
			{
				double K = vel[d];
				val[p + d] += K / rhoVal[scalarp];
			}
		}
	}
}

// #############################
//  accelerator
// #############################

static inline void neInterp3D1(double *result, const double *pos,
							   const double *val, const long int *sizeProd)
{

	// Integer parts of position
	int j = (int)pos[0];
	int k = (int)pos[1];
	int l = (int)pos[2];

	// Decimal (cell-referenced) parts of position and their complement
	double x = pos[0] - j;
	double y = pos[1] - k;
	double z = pos[2] - l;
	double xcomp = 1 - x;
	double ycomp = 1 - y;
	double zcomp = 1 - z;

	// Index of neighbouring nodes
	long int p = j * 3 + k * sizeProd[2] + l * sizeProd[3];
	long int pj = p + 3; // sizeProd[1];
	long int pk = p + sizeProd[2];
	long int pjk = pk + 3; // sizeProd[1];
	long int pl = p + sizeProd[3];
	long int pjl = pl + 3; // sizeProd[1];
	long int pkl = pl + sizeProd[2];
	long int pjkl = pkl + 3; // sizeProd[1];

	// Linear interpolation
	for (int v = 0; v < 3; v++)
		result[v] = zcomp * (ycomp * (xcomp * val[p + v] + x * val[pj + v]) + y * (xcomp * val[pk + v] + x * val[pjk + v]))
				  + z * (ycomp * (xcomp * val[pl + v] + x * val[pjl + v]) + y * (xcomp * val[pkl + v] + x * val[pjkl + v]));
}

static inline void neInterp3D1scalar(double *result, const double *pos,
									 const double *val, const long int *sizeProd)
{

	result[0] = 0;
	// Integer parts of position
	int j = (int)pos[0];
	int k = (int)pos[1];
	int l = (int)pos[2];

	// Decimal (cell-referenced) parts of position and their complement
	double x = pos[0] - j;
	double y = pos[1] - k;
	double z = pos[2] - l;
	double xcomp = 1 - x;
	double ycomp = 1 - y;
	double zcomp = 1 - z;

	// Index of neighbouring nodes
	long int p = j * sizeProd[1] + k * sizeProd[2] + l * sizeProd[3];
	long int pj = p + sizeProd[1]; // sizeProd[1];
	long int pk = p + sizeProd[2];
	long int pjk = pk + sizeProd[1]; // sizeProd[1];
	long int pl = p + sizeProd[3];
	long int pjl = pl + sizeProd[1]; // sizeProd[1];
	long int pkl = pl + sizeProd[2];
	long int pjkl = pkl + sizeProd[1]; // sizeProd[1];

	// Linear interpolation
	// for(int v=0;v<3;v++)
	result[0] = zcomp * (ycomp * (xcomp * val[p] + x * val[pj]) + y * (xcomp * val[pk] + x * val[pjk]))
			  + z * (ycomp * (xcomp * val[pl] + x * val[pjl]) + y * (xcomp * val[pkl] + x * val[pjkl]));
}

funPtr neAcc3D1_set(dictionary *ini)
{
	puSanity(ini, "puAcc3D1", 3, 1); // reuse puSanity for now
	return neAcc3D1;
}

void neAcc3D1(NeutralPopulation *pop, Grid *Pgrad, Grid *divBulkV)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = Pgrad->sizeProd;
	double *val = Pgrad->val;
	long int *divSizeProd = divBulkV->sizeProd;
	double *divVal = divBulkV->val;

	for (int s = 0; s < nSpecies; s++)
	{

		long int pStart = pop->iStart[s] * nDims;
		long int pStop = pop->iStop[s] * nDims;

		for (long int p = pStart; p < pStop; p += nDims)
		{
			double dv[3];
			double divergence[3];
			neInterp3D1(dv, &pos[p], val, sizeProd);
			neInterp3D1(divergence, &pos[p], divVal, divSizeProd); // grad not div
			for (int d = 0; d < nDims; d++)
			{
				// TODO: What if equal to 1? Treated like 0, but without warning?
				if (divergence[d] > 1.)
				{
					divergence[d] = 0;
					msg(WARNING, "Huge divergence encountered, correcting for safety, energy is not conserved");
				}
				if (divergence[d] > 1e-61 || divergence[d] < 1.)
				{
					vel[p + d] += (dv[d] + dv[d] * divergence[d]); //+vel[p+d]*divergence[d]);
				}
				else
				{
					vel[p + d] += dv[d];
				}

				if (vel[p + d] > 1.)
				{
					vel[p + d] = 0.001;
				}
			}
		}
	}
}

// ###########################
//  Mover
// ###########################

void neMove(NeutralPopulation *pop, Grid *V)
{

	int nSpecies = pop->nSpeciesNeutral;

	double *val = V->val;
	long int *sizeProd = V->sizeProd;

	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	for (int s = 0; s < nSpecies; s++)
	{

		long int pStart = pop->iStart[s] * nDims;
		long int pStop = pop->iStop[s] * nDims;

		for (long int p = pStart; p < pStop; p += nDims)
		{
			double dv[3];
			neInterp3D1(dv, &pos[p], val, sizeProd);
			for (int d = 0; d < nDims; d++)
			{
				// TODO: dv suggests it is acceleration
				vel[p + d] = dv[d];	 // val[index+d];
				pos[p + d] += dv[d]; // val[index+d];
			}
		}
	}
}

/******************************************************************************
 *	FINITE DIFFERENCE
 *****************************************************************************/

void divFinDiff1st(Grid *result, Grid *field, Grid *rho)
{

	// Performs first order centered finite difference on field and returns a scalar
	gZero(result);
	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	long int *fieldSizeProd = field->sizeProd;

	double *resultVal = result->val;
	double *fieldVal = field->val;
	double *rhoVal = rho->val;

	// Scalar indices
	long int fNext, fPrev;
	long int s;

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	long int scalarStart = alSum(&rho->sizeProd[1], rank - 1);
	long int scalarEnd = rho->sizeProd[rank] - scalarStart;

	// Centered Finite difference
	for (int d = 1; d < rank; d++)
	{
		fNext = start + fieldSizeProd[d] + (d - 1); //*fieldSizeProd[1]
		fPrev = start - fieldSizeProd[d] + (d - 1); //*fieldSizeProd[1]
		s = scalarStart;

		for (int g = start; g < end; g += rank - 1)
		{
			if (rhoVal[s] < 1e-62)
			{
				// TODO: Is it really needed?
				resultVal[g + d - 1] += 0;
				msg(WARNING, "Low density encountered %f, on node %li energy will not be conserved!", rhoVal[s], s);
			}
			// TODO: What about lower values?
			if (rhoVal[s] >= 1e-12)
			{
				resultVal[g + d - 1] += (1. / (2 * rhoVal[s])) * (fieldVal[fNext] - fieldVal[fPrev]);
			}
			fNext += rank - 1;
			fPrev += rank - 1;
			s++; // fNext;

			if (s > scalarEnd)
			{
				msg(ERROR, "index out of bounds in divFinDiff1st");
			}
			if (fNext > fieldSizeProd[rank])
			{
				msg(ERROR, "index out of bounds in divFinDiff1st, index: %li max: %li, scalar index: %li",
					fNext, fieldSizeProd[rank], s);
			}
		}
	}
	return;
}

// #############################
//  Migration
// #############################

funPtr neExtractEmigrants3DOpen_set(const dictionary *ini)
{
	int nDims = iniGetInt(ini, "grid:nDims");
	if (nDims != 3)
		msg(ERROR, "neExtractEmigrants3DOpen requires grid:nDims=3");
	return neExtractEmigrants3DOpen;
}
void neExtractEmigrants3DOpen(NeutralPopulation *pop, MpiInfo *mpiInfo)
{

	// TODO: Needs more testing
	int nSpecies = pop->nSpeciesNeutral;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *thresholds = mpiInfo->thresholds;
	const int neighborhoodCenter = 13;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;
	int *trueSize = mpiInfo->trueSize;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;
	int *nSubdomains = malloc(3 * sizeof(*nSubdomains));
	bndType *bnd = pop->bnd;

	double dummyPos[3];

	int *offset = mpiInfo->offset;

	nSubdomains[0] = nSubdomainsProd[1];
	nSubdomains[1] = nSubdomainsProd[2] / nSubdomainsProd[1];
	nSubdomains[2] = (nSubdomainsProd[3] / nSubdomainsProd[2]);

	double **emigrants = mpiInfo->emigrantsDummy;
	for (int ne = 0; ne < nNeighbors; ne++)
	{
		emigrants[ne] = mpiInfo->emigrants[ne];
	}
	alSetAll(nEmigrants, nSpecies * nNeighbors, 0);

	double lx = thresholds[0];
	double ly = thresholds[1];
	double lz = thresholds[2];
	double ux = thresholds[3];
	double uy = thresholds[4];
	double uz = thresholds[5];


	for (int s = 0; s < nSpecies; s++)
	{
		long int pStart = pop->iStart[s] * 3;
		long int pStop = pop->iStop[s] * 3;
		long int removedupp = 0; // debug
		long int removedlow = 0; // debug
		long int exhanged = 0;	 // debug
		for (long int p = pStart; p < pStop; p += 3)
		{

			for (int d = 0; d < 3; d++)
			{
				dummyPos[d] = pop->pos[p + d] + offset[d];
			}

			if ((dummyPos[0] > trueSize[0] * (nSubdomains[0]) && bnd[5] != PERIODIC) ||
				(dummyPos[1] > trueSize[1] * (nSubdomains[1]) && bnd[6] != PERIODIC) ||
				(dummyPos[2] > trueSize[2] * (nSubdomains[2]) && bnd[7] != PERIODIC))
			{
				removedupp += 1; // debug
				pos[p] = pos[pStop - 3];
				pos[p + 1] = pos[pStop - 2];
				pos[p + 2] = pos[pStop - 1];
				vel[p] = vel[pStop - 3];
				vel[p + 1] = vel[pStop - 2];
				vel[p + 2] = vel[pStop - 1];

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;
			}
			else if ((dummyPos[0] < -1. && bnd[1] != PERIODIC) ||
					 (dummyPos[1] < -1. && bnd[2] != PERIODIC) ||
					 (dummyPos[2] < -1. && bnd[3] != PERIODIC))
			{

				removedlow += 1; // debug
				pos[p] = pos[pStop - 3];
				pos[p + 1] = pos[pStop - 2];
				pos[p + 2] = pos[pStop - 1];
				vel[p] = vel[pStop - 3];
				vel[p + 1] = vel[pStop - 2];
				vel[p + 2] = vel[pStop - 1];

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;
			}
			else if ((pos[p] >= ux - 1 && bnd[5] == PERIODIC) ||
					 (pos[p + 1] >= uy - 1 && bnd[6] == PERIODIC) ||
					 (pos[p + 2] >= uz - 1 && bnd[7] == PERIODIC) ||
					 (pos[p] < lx + 1 && bnd[1] == PERIODIC) ||
					 (pos[p + 1] < ly + 1 && bnd[2] == PERIODIC) ||
					 (pos[p + 2] < lz + 1 && bnd[3] == PERIODIC))
			{
				// Should look for a better implementation ^

				double x = pos[p];
				double y = pos[p + 1];
				double z = pos[p + 2];
				int nx = -(x < lx) + (x >= ux);
				int ny = -(y < ly) + (y >= uy);
				int nz = -(z < lz) + (z >= uz);
				int ne = neighborhoodCenter + nx + 3 * ny + 9 * nz;

				if (ne != neighborhoodCenter)
				{
					exhanged += 1;
					*(emigrants[ne]++) = x;
					*(emigrants[ne]++) = y;
					*(emigrants[ne]++) = z;
					*(emigrants[ne]++) = vel[p];
					*(emigrants[ne]++) = vel[p + 1];
					*(emigrants[ne]++) = vel[p + 2];
					nEmigrants[ne * nSpecies + s]++;

					pos[p] = pos[pStop - 3];
					pos[p + 1] = pos[pStop - 2];
					pos[p + 2] = pos[pStop - 1];
					vel[p] = vel[pStop - 3];
					vel[p + 1] = vel[pStop - 2];
					vel[p + 2] = vel[pStop - 1];

					pStop -= 3;
					p -= 3;
					pop->iStop[s]--;
				}
			}
		}
	}
	free(nSubdomains);
}

static inline void neShiftImmigrants(MpiInfo *mpiInfo, Grid *grid, int ne, int nSpecies)
{

	double *immigrants = mpiInfo->immigrants;
	long int nImmigrantsTotal = alSum(&mpiInfo->nImmigrants[ne * nSpecies], nSpecies);
	int nDims = mpiInfo->nDims;

	for (int d = 0; d < nDims; d++)
	{
		int n = ne % 3 - 1;
		// TODO: why?
		ne /= 3;

		double shift = n * grid->trueSize[d + 1];

		for (int i = 0; i < nImmigrantsTotal; i++)
		{
			immigrants[d + 2 * nDims * i] += shift;

			double pos = immigrants[d + 2 * nDims * i];

			if (pos > grid->trueSize[d + 1] + 2)
			{
				msg(ERROR, "particle %i skipped two domains, pos = %f", i, pos);
			}
		}
	}
}

static inline void neImportParticles(NeutralPopulation *pop, double *particles, long int *nParticles, int nSpecies)
{

	int nDims = pop->nDims;
	long int *iStop = pop->iStop;

	for (int s = 0; s < nSpecies; s++)
	{
		double *pos = &pop->pos[nDims * iStop[s]];
		double *vel = &pop->vel[nDims * iStop[s]];

		for (int i = 0; i < nParticles[s]; i++)
		{
			for (int d = 0; d < nDims; d++)
				*(pos++) = *(particles++);
			for (int d = 0; d < nDims; d++)
			{
				*(vel++) = *(particles++);
			}
		}

		iStop[s] += nParticles[s];
	}
}

static inline void neExchangeMigrants(NeutralPopulation *pop, MpiInfo *mpiInfo, Grid *grid)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nNeighbors = mpiInfo->nNeighbors;
	long int nImmigrantsAlloc = mpiInfo->nImmigrantsAlloc;
	int nDims = mpiInfo->nDims;
	double **emigrants = mpiInfo->emigrants;
	double *immigrants = mpiInfo->immigrants;
	long int *nImmigrants = mpiInfo->nImmigrants;
	MPI_Request *send = mpiInfo->send;

	for (int ne = 0; ne < nNeighbors; ne++)
	{
		if (ne != mpiInfo->neighborhoodCenter)
		{
			int rank = puNeighborToRank(mpiInfo, ne);
			int reciprocal = puNeighborToReciprocal(ne, nDims);
			long int *nEmigrants = &mpiInfo->nEmigrants[nSpecies * ne];
			long int length = alSum(nEmigrants, nSpecies) * 2 * nDims;
			MPI_Isend(emigrants[ne], length, MPI_DOUBLE, rank, reciprocal, MPI_COMM_WORLD, &send[ne]);
		}
	}

	// Since "immigrants" is reused for every receive operation MPI_Irecv cannot
	// be used. However, in order to receive and process whichever comes first
	// MPI_ANY_SOURCE is used.
	for (int a = 0; a < nNeighbors - 1; a++)
	{

		MPI_Status status;
		MPI_Recv(immigrants, nImmigrantsAlloc, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int ne = status.MPI_TAG; // Which neighbor it is from equals the tag

		neShiftImmigrants(mpiInfo, grid, ne, nSpecies);
		neImportParticles(pop, immigrants, &nImmigrants[ne * nSpecies], nSpecies);
	}

	MPI_Waitall(nNeighbors, send, MPI_STATUS_IGNORE);
}

// TODO: Add fault-handling in case of too small Population struct
static inline void neExchangeNMigrants(NeutralPopulation *pop, MpiInfo *mpiInfo)
{

	int nSpecies = pop->nSpeciesNeutral;
	int nNeighbors = mpiInfo->nNeighbors;
	MPI_Request *send = mpiInfo->send;
	MPI_Request *recv = mpiInfo->recv;

	// Send number of emigrants and receive number of immigrants.
	// Order of reception is not necessarily as determined by the for-loop since
	// we're using non-blocking send/receive.
	for (int ne = 0; ne < nNeighbors; ne++)
	{
		if (ne != mpiInfo->neighborhoodCenter)
		{
			int rank = puNeighborToRank(mpiInfo, ne);
			int reciprocal = puNeighborToReciprocal(ne, mpiInfo->nDims);
			long int *nEmigrants = &mpiInfo->nEmigrants[nSpecies * ne];
			long int *nImmigrants = &mpiInfo->nImmigrants[nSpecies * ne];
			// We get segfault on large simulation systems, as a safety measure we can
			// try to use a different tag here than in exchangeMigrants()
			// hence the "2*nNeighbors"
			MPI_Isend(nEmigrants, nSpecies, MPI_LONG, rank, reciprocal + 2 * nNeighbors, MPI_COMM_WORLD, &send[ne]);
			MPI_Irecv(nImmigrants, nSpecies, MPI_LONG, rank, ne + 2 * nNeighbors, MPI_COMM_WORLD, &recv[ne]);
		}
	}

	MPI_Waitall(nNeighbors, send, MPI_STATUS_IGNORE);
	MPI_Waitall(nNeighbors, recv, MPI_STATUS_IGNORE);
}

void neMigrate(NeutralPopulation *pop, MpiInfo *mpiInfo, Grid *grid)
{

	neExchangeNMigrants(pop, mpiInfo);
	neExchangeMigrants(pop, mpiInfo, grid);
}

void nePNew(NeutralPopulation *pop, int s, const double *pos, const double *vel)
{

	int nDims = pop->nDims;
	long int *iStart = pop->iStart;
	long int *iStop = pop->iStop; // New particle added here

	if (iStop[s] >= iStart[s + 1])
		msg(WARNING, "Not enough allocated memory to add new particle to specie"
					 "%i. New particle ignored.",
			s);
	else
	{

		long int p = iStop[s] * nDims;
		for (int d = 0; d < nDims; d++)
		{
			pop->pos[p + d] = pos[d];
			pop->vel[p + d] = vel[d];
		}

		iStop[s]++;
	}
}

// #############################
//  Pressure solver
// #############################

void neSetI(Grid *I, Grid *V, Grid *rho, const dictionary *ini)
{

	gZero(I);
	double *IEVal = I->val;
	double *rhoVal = rho->val;

	int rank = I->rank;
	long int *sizeProd = I->sizeProd;
	long int *fieldSizeProd = V->sizeProd;

	int nDims = rank - 1;
	int nSpecies = iniGetInt(ini, "collisions:nSpeciesNeutral");
	double *velTherm = iniGetDoubleArr(ini, "collisions:thermalVelocityNeutrals", nDims * nSpecies);

	long int f;
	long int s;

	// TODO: If 0 is hardcoded, is the subtraction necessary?
	long int start = 0;
	long int end = sizeProd[rank] - start;

	long int fieldstart = 0;

	// Centered Finite difference
	for (int d = 1; d < rank; d++)
	{
		f = fieldstart + (d - 1);

		s = start;

		for (int g = start; g < end; g++)
		{
			// TODO: Why recalculate most of it all the time?
			IEVal[g] = 0.5 * ((rhoVal[g])) * velTherm[0] * velTherm[0];
			if (g > sizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neSetI");
			}
			if (f > fieldSizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neSetI, index: %li max: %li, scalar index: %li", f, fieldSizeProd[rank], s);
			}
			f += rank - 1;

			s++; // fNext;

		}
	}
}

void nePressureSolve3D(Grid *P, Grid *IE, Grid *rho, NeutralPopulation *pop)
{
	gZero(P);

	double *mass = pop->mass;
	int *nGhostLayers = P->nGhostLayers;
	int *trueSize = P->trueSize;
	long int *sizeProd = P->sizeProd;
	double rho0 = pop->rho0;
	// Seperate values
	double *PVal = P->val;
	double *IEVal = IE->val;
	double *rhoVal = rho->val;

	double gamma = 5. / 3.; // Adiabatic index

	int index = 0;

	for (int k = nGhostLayers[1]; k < trueSize[3] + nGhostLayers[7]; k++)
	{
		for (int j = nGhostLayers[2]; j < trueSize[2] + nGhostLayers[6]; j++)
		{
			for (int i = nGhostLayers[3]; i < trueSize[1] + nGhostLayers[5]; i++)
			{
				index = i + j * sizeProd[2] + k * sizeProd[3];
				if (rhoVal[index] == 0)
				{
					PVal[index] = 0;
				}
				// TODO: else? It will be multiplied by zero anyway, is it really needed at all?
				if (rhoVal[index] != 0)
				{
					// should be multiplied by massss for several species
					PVal[index] = (gamma - 1) * (rhoVal[index] / rho0) * mass[0] * IEVal[index];
				}
			}
		}
	}
}

void neAdvectI(Grid *IE, Grid *Itilde, Grid *P, Grid *V, Grid *rho, NeutralPopulation *pop)
{
	gZero(Itilde);
	double *PVal = P->val;
	double *IEVal = IE->val;
	double *ItildeVal = Itilde->val;
	double *bulkVVal = V->val;
	double *rhoVal = rho->val;
	double *mass = pop->mass;
	int nDims = pop->nDims;

	int rank = IE->rank;
	long int *sizeProd = IE->sizeProd;
	long int *fieldSizeProd = V->sizeProd;

	long int fNext, fPrev, sNext, sPrev;
	long int s, f;

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	long int fieldstart = alSum(&fieldSizeProd[1], rank - 1);

	// Centered Finite difference
	for (int d = 1; d < rank; d++)
	{
		fNext = fieldstart + fieldSizeProd[d] + (d - 1);
		fPrev = fieldstart - fieldSizeProd[d] + (d - 1);

		sNext = start + sizeProd[d];
		sPrev = start - sizeProd[d];

		s = start;
		f = fieldstart + (d - 1);

		for (int g = start; g < end; g++)
		{
			if (rhoVal[g] == 0.)
			{
				ItildeVal[g] = 0;
			}
			// TODO: else?
			if (rhoVal[g] != 0.)
			{
				ItildeVal[g] += IEVal[g] / nDims + (1. / (2 * mass[0] * rhoVal[g]))
												 * (PVal[sPrev] * (bulkVVal[fPrev] - bulkVVal[f])
												  - PVal[sNext] * (bulkVVal[fNext] - bulkVVal[f]));
			}
			fNext += rank - 1;
			fPrev += rank - 1;
			sNext++;
			sPrev++;
			s++; // fNext;
			f += rank - 1;

			if (g > sizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve");
			}
			if (fNext > fieldSizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve, index: %li max: %li, scalar index: %li",
					fNext, fieldSizeProd[rank], s);
			}
		}
	}
}

void neAdvectV(Grid *V, Grid *Vtilde, Grid *P, Grid *rho, NeutralPopulation *pop)
{
	gZero(Vtilde);
	double *PVal = P->val;
	double *bulkVVal = V->val;
	double *bulkVtildeVal = Vtilde->val;
	double *rhoVal = rho->val;
	double *mass = pop->mass;

	int rank = P->rank;
	long int *sizeProd = P->sizeProd;
	long int *fieldSizeProd = V->sizeProd;

	long int fNext, fPrev, sNext, sPrev;
	long int s, f;

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	long int fieldstart = alSum(&fieldSizeProd[1], rank - 1);

	// Centered Finite difference
	for (int d = 1; d < rank; d++)
	{
		fNext = fieldstart + fieldSizeProd[d] + (d - 1);
		fPrev = fieldstart - fieldSizeProd[d] + (d - 1);

		sNext = start + sizeProd[d];
		sPrev = start - sizeProd[d];

		s = start;
		f = fieldstart;

		for (int g = start; g < end; g++)
		{
			if (rhoVal[g] == 0.)
			{
				bulkVtildeVal[f + d - 1] = 0;
			}
			// TODO: else?
			if (rhoVal[g] != 0.)
			{
				bulkVtildeVal[f + d - 1] = bulkVVal[f + d - 1] + (1. / (2 * mass[0] * rhoVal[g])) * (PVal[sPrev] - PVal[sNext]);
			}
			fNext += rank - 1;
			fPrev += rank - 1;
			sNext++;
			sPrev++;
			s++; // fNext;
			f += rank - 1;
			if (g > sizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve");
			}
			if (fNext > fieldSizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve, index: %li max: %li, scalar index: %li",
					fNext, fieldSizeProd[rank], s);
			}
		}
	}
}

void neConvectV(Grid *V, Grid *Vtilde, Grid *rhoNeutral, NeutralPopulation *pop)
{

	gCopy(Vtilde, V); // V <--- Vtilde

	int nSpecies = pop->nSpeciesNeutral;

	double *vVal = V->val;
	double *vtildeVal = Vtilde->val;
	double *rhoVal = rhoNeutral->val;
	long int *fieldSizeProd = V->sizeProd;
	long int *sizeProd = rhoNeutral->sizeProd;

	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	for (int s = 0; s < nSpecies; s++)
	{

		long int pStart = pop->iStart[s] * nDims;
		long int pStop = pop->iStop[s] * nDims;

		for (long int p = pStart; p < pStop; p += nDims)
		{

			int jPrev = (int)(pos[p] - vel[p] + 0.5);
			int kPrev = (int)(pos[p + 1] - vel[p + 1] + 0.5);
			int lPrev = (int)(pos[p + 2] - vel[p + 2] + 0.5);

			int j = (int)(pos[p] + 0.5);
			int k = (int)(pos[p + 1] + 0.5);
			int l = (int)(pos[p + 2] + 0.5);

			long int indexPrev = jPrev + kPrev * sizeProd[2] + lPrev * sizeProd[3];
			long int index = j + k * sizeProd[2] + l * sizeProd[3];

			if (indexPrev != index)
			{
				long int fieldIndexPrev = jPrev * fieldSizeProd[1] + kPrev * fieldSizeProd[2] + lPrev * fieldSizeProd[3];
				long int fieldIndex = j * fieldSizeProd[1] + k * fieldSizeProd[2] + l * fieldSizeProd[3];
				for (int d = 0; d < nDims; d++)
				{
					vVal[fieldIndex + d] = (rhoVal[index] * vtildeVal[fieldIndex + d] + vtildeVal[fieldIndexPrev + d])
										 / (rhoVal[index] + 1);
				}
				rhoVal[index] += 1;
				rhoVal[indexPrev] -= 1;
			}
		}
	}
}

void neConvectKE(Grid *dKE, Grid *Vtilde, Grid *rhoNeutral, NeutralPopulation *pop)
{

	gZero(dKE);

	int nSpecies = pop->nSpeciesNeutral;

	double *dKEVal = dKE->val;
	double *vtildeVal = Vtilde->val;
	double *rhoVal = rhoNeutral->val;
	long int *fieldSizeProd = Vtilde->sizeProd;
	long int *sizeProd = dKE->sizeProd;

	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;

	for (int s = 0; s < nSpecies; s++)
	{

		long int pStart = pop->iStart[s] * nDims;
		long int pStop = pop->iStop[s] * nDims;

		for (long int p = pStart; p < pStop; p += nDims)
		{

			int jPrev = (int)(pos[p] - vel[p] + 0.5);
			int kPrev = (int)(pos[p + 1] - vel[p + 1] + 0.5);
			int lPrev = (int)(pos[p + 2] - vel[p + 2] + 0.5);

			int j = (int)(pos[p] + 0.5);
			int k = (int)(pos[p + 1] + 0.5);
			int l = (int)(pos[p + 2] + 0.5);

			long int indexPrev = jPrev + kPrev * sizeProd[2] + lPrev * sizeProd[3];
			long int index = j + k * sizeProd[2] + l * sizeProd[3];
			if (indexPrev != index)
			{
				long int fieldIndexPrev = jPrev * fieldSizeProd[1] + kPrev * fieldSizeProd[2] + lPrev * fieldSizeProd[3];
				long int fieldIndex = j * fieldSizeProd[1] + k * fieldSizeProd[2] + l * fieldSizeProd[3];
				for (int d = 0; d < nDims; d++)
				{
					dKEVal[index] -= 0.5 * mass[0] * (rhoVal[index] / (rhoVal[index] + 1))
								   * pow((vtildeVal[fieldIndexPrev + d] - vtildeVal[fieldIndex + d]), 2);
				}
			}
		}
	}
}

void neConvectI(Grid *IE, Grid *Itilde, Grid *dKE, Grid *rhoNeutral, NeutralPopulation *pop)
{

	gCopy(Itilde, IE); // IE <--- Itilde

	int nSpecies = pop->nSpeciesNeutral;

	double *dKEVal = dKE->val;
	double *IVal = IE->val;
	double *ItildeVal = Itilde->val;
	double *rhoVal = rhoNeutral->val;
	long int *sizeProd = dKE->sizeProd;

	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;

	for (int s = 0; s < nSpecies; s++)
	{

		long int pStart = pop->iStart[s] * nDims;
		long int pStop = pop->iStop[s] * nDims;

		for (long int p = pStart; p < pStop; p += nDims)
		{

			int jPrev = (int)(pos[p] - vel[p] + 0.5);
			int kPrev = (int)(pos[p + 1] - vel[p + 1] + 0.5);
			int lPrev = (int)(pos[p + 2] - vel[p + 2] + 0.5);

			int j = (int)(pos[p] + 0.5);
			int k = (int)(pos[p + 1] + 0.5);
			int l = (int)(pos[p + 2] + 0.5);

			long int indexPrev = jPrev + kPrev * sizeProd[2] + lPrev * sizeProd[3];
			long int index = j + k * sizeProd[2] + l * sizeProd[3];
			if (indexPrev != index)
			{
				IVal[index] = (rhoVal[index] * IVal[index] + ItildeVal[indexPrev]) / (rhoVal[index] + 1) + 0;
				if (rhoVal[indexPrev] == 0.0)
				{
					IVal[indexPrev] = rhoVal[indexPrev] * ItildeVal[indexPrev];
				}
			}
		}
	}
	for (int i = 0; i < sizeProd[nDims + 1]; i++)
	{
		IVal[i] += sqrt(pow(dKEVal[i], 2)) / ((rhoVal[i] + 1) * mass[0]);
	}
}

void neAddPressure(Grid *bulkV, Grid *Pgrad, Grid *rho)
{

	int rank = bulkV->rank;
	long int *fieldSizeProd = bulkV->sizeProd;
	long int *sizeProd = rho->sizeProd;
	double *bulkVVal = bulkV->val;
	double *PgradVal = Pgrad->val;
	double *rhoVal = rho->val;

	// indices
	long int s;
	long int f;
	int fNext = fieldSizeProd[1];

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	// Excluding edge values where Pgrad is not defined
	for (int d = 1; d < rank; d++)
	{
		s = start;
		f = start * fieldSizeProd[1] + (d - 1);

		for (int g = start; g < end; g++)
		{
			bulkVVal[f] += (1. / (2. * rhoVal[s])) * PgradVal[f];
			s++;
			f += fNext;
		}
	}
}

/***********************************
 *	Boundary functions
 ***********************************/

void neSetBndSlices(Grid *grid, const MpiInfo *mpiInfo)
{

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	// Number of elements in slice
	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}

	// TODO: Always 0?
	double constant2 = 0.; // solution to bnd cond = 0.

	// Lower edge
	for (int d = 1; d < rank; d++)
	{
		if (subdomain[d - 1] == 0)
		{
			if (bnd[d] == DIRICHLET)
			{
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * d)] = 0.;
				}
			}
			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					// Solution to equation. constant for now
				}
			}
		}
	}

	// Upper edge
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (subdomain[d - rank - 1] == nSubdomains[d - rank - 1] - 1)
		{
			if (bnd[d] == DIRICHLET)
			{

				for (int s = 0; s < nSliceMax; s++)
				{

					bndSlice[s + (nSliceMax * (d))] = 0;
					bndSlice[s + (nSliceMax * (d))] = 0.;
				}
			}

			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	return;
}

void neSetBndSlicesEnerg(const dictionary *ini, Grid *grid, Grid *rho, const MpiInfo *mpiInfo)
{

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	double *sendSlice = rho->sendSlice;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.
	int nSpecies = iniGetInt(ini, "collisions:nSpeciesNeutral");
	double *velTherm = iniGetDoubleArr(ini, "collisions:thermalVelocityNeutrals", nDims * nSpecies);

	double veld[nDims];

	for (int s = 0; s < nSpecies; s++)
	{
		for (int d = 0; d < nDims; d++)
		{
			veld[d] += (1. / nSpecies) * velTherm[d + s * nDims];
		}
	}

	double Dspeed = 0;
	for (int d = 0; d < nDims; d++)
	{
		Dspeed += sqrt(veld[d] * veld[d]);
	}

	// Number of elements in slice
	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}

	double constant2 = 0.; // solution to bnd cond = 0.

	// Lower edge
	for (int d = 1; d < rank; d++)
	{
		if (subdomain[d - 1] == 0)
		{
			if (bnd[d] == DIRICHLET)
			{
				getSlice(sendSlice, rho, d, 1);
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * d)] = 0.5 * (sendSlice[s]) * velTherm[0] * velTherm[0];
				}
			}
			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					// Solution to equation. constant for now
				}
			}
		}
	}

	// Upper edge
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (subdomain[d - rank - 1] == nSubdomains[d - rank - 1] - 1)
		{
			if (bnd[d] == DIRICHLET)
			{
				getSlice(sendSlice, rho, d - rank, size[d - rank] - 2);
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * (d))] = 0.5 * (sendSlice[s]) * velTherm[0] * velTherm[0];
				}
			}

			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	return;
}

void neSetBndSlicesVel(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo)
{

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	// int *nGhostLayers = grid->nGhostLayers;
	// double *bndSolution = grid->bndSolution;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.
	int nSpecies = iniGetInt(ini, "collisions:nSpeciesNeutral");
	double *velDrift = iniGetDoubleArr(ini, "collisions:neutralDrift", nDims * nSpecies);

	double veld[nDims];

	for (int s = 0; s < nSpecies; s++)
	{
		for (int d = 0; d < nDims; d++)
		{
			// msg(STATUS,"d = %i, d+s*nDims = %i",d,d+s*nDims);
			veld[d] += (1. / nSpecies) * velDrift[d + s * nDims];
		}
	}

	// Number of elements in slice
	long int nSliceMax = 0;
	for (int d = 1; d < rank; d++)
	{
		long int nSlice = 1;
		for (int dd = 0; dd < rank; dd++)
		{
			if (dd != d)
				nSlice *= size[dd];
		}
		if (nSlice > nSliceMax)
			nSliceMax = nSlice;
	}
	nSliceMax = nSliceMax * nDims;

	// double constant1 = 0.; //Dirichlet
	double constant2 = 0.; // solution to bnd cond = 0.

	// Lower edge
	for (int d = 1; d < rank; d++)
	{
		if (subdomain[d - 1] == 0)
		{
			if (bnd[d] == DIRICHLET)
			{
				for (int s = 0; s < nSliceMax; s += nDims)
				{
					for (int dd = 0; dd < rank - 1; dd++)
					{ // dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + dd + (nSliceMax * d)] = veld[dd];
					}
				}
			}
			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					// Solution to equation. constant for now
				}
			}
		}
	}

	// Upper edge
	for (int d = rank + 1; d < 2 * rank; d++)
	{
		if (subdomain[d - rank - 1] == nSubdomains[d - rank - 1] - 1)
		{
			if (bnd[d] == DIRICHLET)
			{

				for (int s = 0; s < nSliceMax; s += nDims)
				{

					bndSlice[s + (nSliceMax * (d))] = 0;
					for (int dd = 0; dd < rank - 1; dd++)
					{
						bndSlice[s + dd + (nSliceMax * d)] = veld[dd];
					}
				}
			}

			if (bnd[d] == NEUMANN)
			{
				for (int s = 0; s < nSliceMax; s++)
				{
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	// msg(STATUS,"nSliceMax = %li",nSliceMax);
	// adPrint(&bndSlice[nSliceMax], nSliceMax*(rank));
	free(velDrift);
	// free(B);
	return;
}

// #########################################
//  Object functions
// #########################################
void neIndexToPos3D(Grid *grid, long int index, long int *pos)
{

	long int *sizeProd = grid->sizeProd;
	long int i, j, k;
	long int p = index; /// sizeProd[1];

	k = (int)(p / sizeProd[3]);
	j = (int)((p - k * sizeProd[3]) / (sizeProd[2]));
	i = (int)(p - j * sizeProd[2] - k * sizeProd[3]) / sizeProd[1];
	// printf("pos = %li,%li,%li \n",i,j,k);
	pos[0] = i;
	pos[1] = j;
	pos[2] = k;
}

void neApplyObjVel(PincObject *obj, Grid *V)
{

	// neInternalEnergySolve(IE,P,bulkV,rhoNeutral,neutralPop);
	double *val = V->val;
	// int nDims = pop->nDims;

	int rank = obj->domain->rank;
	long int *lookupSurfOff = obj->lookupSurfaceOffset;
	long int *lookupIntOff = obj->lookupInteriorOffset;
	long int *sizeProd = obj->domain->sizeProd;
	long int *fieldSizeProd = V->sizeProd;

	// int index = 0;
	// printf("sizeProd[1] = %i, %i, %i, %i \n",sizeProd[1],sizeProd[2],sizeProd[3],sizeProd[4]);
	long int fNext, fPrev, sNext, sPrev;
	long int s, f;
	// int fNext = fieldSizeProd[1];

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	long int fieldstart = alSum(&fieldSizeProd[1], rank - 1);
	// long int fieldend = fieldSizeProd[rank]-start;

	bool prevInside, nextInside;
	for (int d = 1; d < rank; d++)
	{
		// fNext = fieldstart + fieldSizeProd[d]+(d-1);
		// fPrev = fieldstart - fieldSizeProd[d]+(d-1);

		sNext = start + sizeProd[d];
		sPrev = start - sizeProd[d];

		fNext = fieldstart + fieldSizeProd[d] + (d - 1);
		fPrev = fieldstart - fieldSizeProd[d] + (d - 1);

		s = start;
		f = fieldstart + (d - 1);

		// PVal[fPrev] is a scalar not field
		// printf("\n \n");
		for (int g = start; g < end; g++)
		{

			for (long int a = 0; a < obj->nObjects; a++)
			{
				for (long int b = lookupSurfOff[a]; b < lookupSurfOff[a + 1]; b++)
				{
					if ((obj->lookupSurface[b]) == g)
					{
						prevInside = false;
						nextInside = false;
						for (long int a2 = 0; a2 < obj->nObjects; a2++)
						{

							for (long int b2 = lookupIntOff[a2]; b2 < lookupIntOff[a2 + 1]; b2++)
							{
								if ((obj->lookupInterior[b2]) == sNext)
								{
									nextInside = true;
									// printf("setting val \n");
								}
								if ((obj->lookupInterior[b2]) == sPrev)
								{
									// printf("setting val");
									prevInside = true;
								}
							}
						}
						for (long int a2 = 0; a2 < obj->nObjects; a2++)
						{

							for (long int b2 = lookupSurfOff[a2]; b2 < lookupSurfOff[a2 + 1]; b2++)
							{
								if ((obj->lookupSurface[b2]) == sNext)
								{
									nextInside = true;
									// printf("setting val \n");
								}
								if ((obj->lookupSurface[b2]) == sPrev)
								{
									// printf("setting val");
									prevInside = true;
								}
							}
						}
						if (nextInside == true)
						{
							if (prevInside == false)
							{
								val[f] = -1. * sqrt(val[fPrev] * val[fPrev]);
								// long int gridpos[3];
								// long int gridposPrev[3];
								// long int gridposNext[3];
								// neIndexToPos3D(obj->domain,g,gridpos);
								// neIndexToPos3D(obj->domain,sPrev,gridposPrev);
								// neIndexToPos3D(obj->domain,sNext,gridposNext);
								// printf("in dim %i for node %li,%li,%li, using prev=%li,%li,%li, next=%li,%li,%li, \n",d,gridpos[0],gridpos[1],gridpos[2],gridposPrev[0],gridposPrev[1],gridposPrev[2],gridposNext[0],gridposNext[1],gridposNext[2]);
							}
						}
						if (prevInside == true)
						{
							if (nextInside == false)
							{
								val[f] = 1. * sqrt(val[fNext] * val[fNext]);
							}
						}
					}
				}
			}
			// exit(0);
			fNext += rank - 1;
			fPrev += rank - 1;
			sNext++;
			sPrev++;
			s++; // fNext;
			f += rank - 1;
			// printf("IEVal[g] = %f, rhoVal[g] = %f \n",IEVal[g],rhoVal[g] );
			// printf("%f, %f \n",PVal[sPrev], bulkVVal[fPrev]  );
			if (g > sizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve");
			}
			if (fNext > fieldSizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve, index: %li max: %li, scalar index: %li", fNext, fieldSizeProd[rank], s);
			}
			// printf("node: %li, using %li, and %li \n",g,fNext,fPrev);
			// printf("fieldSizeProd = %li, scalarsizeProd = %li \n \n",fieldSizeProd[4],sizeProd[4]);
		}
	}
	// exit(0);
}

void neApplyObjI(PincObject *obj, Grid *IE)
{

	// neInternalEnergySolve(IE,P,bulkV,rhoNeutral,neutralPop);
	double *val = IE->val;
	// int nDims = pop->nDims;

	int rank = obj->domain->rank;
	long int *lookupSurfOff = obj->lookupSurfaceOffset;
	long int *lookupIntOff = obj->lookupInteriorOffset;
	long int *sizeProd = obj->domain->sizeProd;

	// int index = 0;
	// printf("sizeProd[1] = %i, %i, %i, %i \n",sizeProd[1],sizeProd[2],sizeProd[3],sizeProd[4]);
	long int sNext, sPrev;
	long int s;
	// int fNext = fieldSizeProd[1];

	long int start = alSum(&sizeProd[1], rank - 1);
	long int end = sizeProd[rank] - start;

	bool nextInside, prevInside;
	for (int d = 1; d < rank; d++)
	{

		sNext = start + sizeProd[d];
		sPrev = start - sizeProd[d];

		s = start;

		// PVal[fPrev] is a scalar not field
		// printf("\n \n");
		for (int g = start; g < end; g++)
		{

			for (long int a = 0; a < obj->nObjects; a++)
			{
				for (long int b = lookupSurfOff[a]; b < lookupSurfOff[a + 1]; b++)
				{
					if ((obj->lookupSurface[b]) == g)
					{
						// for (long int a2=0; a2<obj->nObjects; a2++) {
						// 	for (long int b2=lookupIntOff[a2]; b2<lookupIntOff[a2+1]; b2++) {
						// 		if ((obj->lookupInterior[b2])==sNext) {
						// 			val[sNext] = 20*val[sPrev];
						// 			//printf("setting val \n");
						// 		}
						// 		if ((obj->lookupInterior[b2])==sPrev) {
						// 			val[sPrev] = 20*val[sNext];
						// 			//printf("setting val");
						// 		}
						// 	}
						// }
						prevInside = false;
						nextInside = false;
						for (long int a2 = 0; a2 < obj->nObjects; a2++)
						{

							for (long int b2 = lookupIntOff[a2]; b2 < lookupIntOff[a2 + 1]; b2++)
							{
								if ((obj->lookupInterior[b2]) == sNext)
								{
									nextInside = true;
									// printf("setting val \n");
								}
								if ((obj->lookupInterior[b2]) == sPrev)
								{
									// printf("setting val");
									prevInside = true;
								}
							}
						}
						for (long int a2 = 0; a2 < obj->nObjects; a2++)
						{

							for (long int b2 = lookupSurfOff[a2]; b2 < lookupSurfOff[a2 + 1]; b2++)
							{
								if ((obj->lookupSurface[b2]) == sNext)
								{
									nextInside = true;
									// printf("setting val \n");
								}
								if ((obj->lookupSurface[b2]) == sPrev)
								{
									// printf("setting val");
									prevInside = true;
								}
							}
						}
						if (nextInside == true)
						{
							if (prevInside == false)
							{
								val[s] = val[sPrev];
								// long int gridpos[3];
								// long int gridposPrev[3];
								// long int gridposNext[3];
								// neIndexToPos3D(obj->domain,g,gridpos);
								// neIndexToPos3D(obj->domain,sPrev,gridposPrev);
								// neIndexToPos3D(obj->domain,sNext,gridposNext);
								// printf("in dim %i for node %li,%li,%li, using prev=%li,%li,%li, next=%li,%li,%li, \n",d,gridpos[0],gridpos[1],gridpos[2],gridposPrev[0],gridposPrev[1],gridposPrev[2],gridposNext[0],gridposNext[1],gridposNext[2]);
							}
						}
						if (prevInside == true)
						{
							if (nextInside == false)
							{
								val[s] = val[sNext];
							}
						}
					}
				}
			}
			// exit(0);

			sNext++;
			sPrev++;
			s++; // fNext;

			// printf("IEVal[g] = %f, rhoVal[g] = %f \n",IEVal[g],rhoVal[g] );
			// printf("%f, %f \n",PVal[sPrev], bulkVVal[fPrev]  );
			if (g > sizeProd[rank])
			{
				msg(ERROR, "index out of bounds in neInternalEnergySolve");
			}
			// printf("node: %li, using %li, and %li \n",g,fNext,fPrev);
			// printf("fieldSizeProd = %li, scalarsizeProd = %li \n \n",fieldSizeProd[4],sizeProd[4]);
		}
	}
	// exit(0);
}

void nuObjectpurge(NeutralPopulation *pop, Grid *rhoObj, PincObject *obj)
{

	// int rank = mpiInfo->mpiRank;
	// int size = mpiInfo->mpiSize;

	// double *val = rhoObj->val;
	long int *sizeProd = rhoObj->sizeProd;
	long int nDims = pop->nDims;

	int nSpecies = pop->nSpeciesNeutral;

	long int *lookupIntOff = obj->lookupInteriorOffset;
	long int *lookupSurfOff = obj->lookupSurfaceOffset;

	int cutNumber = 0;
	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (long int i = iStart; i < iStop; i++)
		{

			double *pos = &pop->pos[i * nDims];
			// double *vel = &pop->vel[i*nDims];

			// Integer parts of position
			int j = (int)pos[0];
			int k = (int)pos[1];
			int l = (int)pos[2];

			long int p = j + k * sizeProd[2] + l * sizeProd[3];
			long int pIndex = i * nDims; // j + k*sizeProd[2] + l*sizeProd[3];
			// msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
			//  Check whether p is one of the object nodes and collect the charge if so.
			for (long int a = 0; a < obj->nObjects; a++)
			{
				for (long int b = lookupIntOff[a]; b < lookupIntOff[a + 1]; b++)
				{
					if ((obj->lookupInterior[b]) == p)
					{
						// msg(STATUS,"p, pIndex: %li,%li, %li",p,(pIndex-iStart*nDims),(iStop-iStart));
						// msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
						// msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);
						neCutParticle(pop, s, pIndex, pop->pos, pop->vel);
						cutNumber += 1;
						// msg(STATUS,"iStop = %li",iStop);
						iStop--;
					}
				}
			}
			for (long int a = 0; a < obj->nObjects; a++)
			{
				for (long int b = lookupSurfOff[a]; b < lookupSurfOff[a + 1]; b++)
				{
					if ((obj->lookupSurface[b]) == p)
					{
						// msg(STATUS,"p, pIndex: %li,%li, %li",p,(pIndex-iStart*nDims),(iStop-iStart));
						// msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
						// msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);
						neCutParticle(pop, s, pIndex, pop->pos, pop->vel);
						cutNumber += 1;
						// msg(STATUS,"iStop = %li",iStop);
						iStop--;
					}
				}
			}
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &cutNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	printf("Neutrals cutNumber = %i \n", cutNumber);
	cutNumber = 0;
}

void nuObjectCollide(NeutralPopulation *pop, Grid *rhoObj, PincObject *obj)
{

	// int rank = mpiInfo->mpiRank;
	// int size = mpiInfo->mpiSize;

	// double *val = rhoObj->val;
	long int *sizeProd = rhoObj->sizeProd;
	long int nDims = pop->nDims;

	int nSpecies = pop->nSpeciesNeutral;

	long int *lookupIntOff = obj->lookupInteriorOffset;
	// long int *lookupSurfOff = obj->lookupSurfaceOffset;

	int cutNumber = 0;
	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for (long int i = iStart; i < iStop; i++)
		{

			double *pos = &pop->pos[i * nDims];
			double *vel = &pop->vel[i * nDims];

			// Integer parts of position
			int j = (int)pos[0];
			int k = (int)pos[1];
			int l = (int)pos[2];

			long int p = j + k * sizeProd[2] + l * sizeProd[3];
			// long int pIndex = i*nDims; //j + k*sizeProd[2] + l*sizeProd[3];
			// msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
			//  Check whether p is one of the object nodes and collect the charge if so.
			//  for (long int a=0; a<obj->nObjects; a++) {
			//  	for (long int b=lookupSurfOff[a]; b<lookupSurfOff[a+1]; b++) {
			//  		if ((obj->lookupSurface[b])==p) {

			for (long int a = 0; a < obj->nObjects; a++)
			{
				for (long int b = lookupIntOff[a]; b < lookupIntOff[a + 1]; b++)
				{
					if ((obj->lookupInterior[b]) == p)
					{
						// msg(STATUS,"j,p, pIndex: %i, %li,%li, %li",j,p,(pIndex-iStart*nDims),(iStop-iStart));
						// msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
						// msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);

						// printf("Before: pos = %f,%f,%f \n",pos[0],pos[1],pos[2]);
						// printf("Before: vel = %f,%f,%f \n",vel[0],vel[1],vel[2]);

						neScatterParticle(pop, pos, vel);

						// printf("After: pos = %f,%f,%f \n",pos[0],pos[1],pos[2]);
						// printf("After: vel = %f,%f,%f \n",vel[0],vel[1],vel[2]);
						// printf("\n");

						cutNumber += 1;
						// msg(STATUS,"iStop = %li",iStop);
						iStop--;

						// Unit testing

						// Integer parts of position
						int j = (int)pos[0];
						int k = (int)pos[1];
						int l = (int)pos[2];

						long int p = j + k * sizeProd[2] + l * sizeProd[3];
						for (long int a = 0; a < obj->nObjects; a++)
						{
							for (long int b = lookupIntOff[a]; b < lookupIntOff[a + 1]; b++)
							{
								if ((obj->lookupInterior[b]) == p)
								{
									msg(WARNING, "Neutral particle-object did not collide correct");
									printf("After: pos = %f,%f,%f \n", pos[0], pos[1], pos[2]);
									printf("After: vel = %f,%f,%f \n", vel[0], vel[1], vel[2]);
									printf("\n");
								}
							}
						} // Unit testing
					}
				}
			}
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &cutNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	msg(STATUS, "collided Number = %i \n", cutNumber);
	cutNumber = 0;
}

void nuObjectSetVal(Grid *rho, double constant, PincObject *obj)
{

	// int rank = mpiInfo->mpiRank;
	// int size = mpiInfo->mpiSize;

	// double *val = rhoObj->val;
	double *rhoVal = rho->val;

	// long int *sizeProd = rhoObj->sizeProd;

	long int *lookupIntOff = obj->lookupInteriorOffset;
	// long int *lookupSurfOff = obj->lookupSurfaceOffset;

	// adPrint(rhoVal,rho->sizeProd[4]);
	// msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
	//  Check whether p is one of the object nodes and collect the charge if so.
	for (long int a = 0; a < obj->nObjects; a++)
	{
		for (long int b = lookupIntOff[a]; b < lookupIntOff[a + 1]; b++)
		{
			rhoVal[obj->lookupInterior[b]] = constant;
			// printf("rhoVal[obj->lookupInterior[b]] = %f \n",rhoVal[obj->lookupInterior[b]] );
		}
	}
	// for (long int a=0; a<obj->nObjects; a++) {
	// 	for (long int b=lookupSurfOff[a]; b<lookupSurfOff[a+1]; b++) {
	// 		rhoVal[obj->lookupSurface[b]] = constant;
	// 		//printf("obj->lookupSurfOff[b]= %li \n",obj->lookupSurface[b]);
	// 	}
	// }
	// adPrint(rhoVal,rho->sizeProd[4]);
	// exit(0);
}

void neVelAssertMax(const NeutralPopulation *pop, double max)
{

	double *vel = pop->vel;

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];
		for (int i = iStart; i < iStop; i++)
		{

			for (int d = 0; d < nDims; d++)
			{

				if (vel[i * nDims + d] > max)
				{
					msg(ERROR, "Particle i=%li (of specie %i) travels too"
							   "fast in dimension %i: %f>%f",
						i, s, d, vel[i * nDims + d], max);
				}
			}
		}
	}
}

void nePosAssertInLocalFrame(const NeutralPopulation *pop, const Grid *grid)
{

	int *size = grid->size;
	double *pos = pop->pos;

	int nSpecies = pop->nSpeciesNeutral;

	int nDims = pop->nDims;

	for (int s = 0; s < nSpecies; s++)
	{

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];
		for (int i = iStart; i < iStop; i++)
		{

			for (int d = 0; d < nDims; d++)
			{

				if (pos[i * nDims + d] > size[d + 1] + 1 || pos[i * nDims + d] < -1)
				{
					printf("Particle i=%i (of neutral specie %i) is out of bounds"
						   "in dimension %i: %f>%i",
						   i, s, d, pos[i * nDims + d], size[d + 1] - 1);
					exit(1);
				}
			}
		}
	}
}
