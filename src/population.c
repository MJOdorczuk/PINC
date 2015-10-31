/**
 * @file		particles.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle handling.
 * @date		26.10.15
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#include "pinc.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

Population *allocPopulation(const dictionary *ini){

	// Sanity check
	iniAssertEqualNElements(ini,4,"population:nParticles","population:nAlloc","population:q","population:m");

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data
	int nSpecies;
	long int *nAllocTotal = iniGetLongIntArr(ini,"population:nAlloc",&nSpecies);	// This is for all computing nodes
	int nDims = iniGetNElements(ini,"grid:nTGPoints");
	if(nDims==0)	msg(ERROR,"grid:nTGPoints not found");

	// Determine memory to allocate for this node
	long int *nAlloc = malloc(nSpecies*sizeof(long int));
	for(int s=0;s<nSpecies;s++){
		nAlloc[s] = ceil(nAllocTotal[s]/size);
		if(nAlloc[s]*size!=nAllocTotal[s])
			msg(WARNING,"increased number of allocated particles from %i to %i to get integer per computing node",
				nAllocTotal[s],nAlloc[s]*size);
	}

	// Determine iStart and iStop
	long int *iStart = malloc((nSpecies+1)*sizeof(long int));
	long int *iStop  = malloc(nSpecies*sizeof(long int));

	iStart[0] = 0;
	for(int s=1;s<nSpecies+1;s++) iStart[s]=iStart[s-1]+nAlloc[s-1];
	for(int s=0;s<nSpecies;s++) iStop[s]=iStart[s]-1; // No particles yet, set stop before iStart.

	free(nAlloc);

	// Store in struct
	Population *pop = malloc(sizeof(Population));
	pop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->q = iniGetDoubleArr(ini,"population:q",&nSpecies);
	pop->m = iniGetDoubleArr(ini,"population:m",&nSpecies);
	pop->nSpecies = nSpecies;
	pop->nDims = nDims;
	pop->iStart = iStart;
	pop->iStop = iStop;
	pop->energy = NULL;	// Assume unused for now

	return pop;

}

void freePopulation(Population *pop){

	free(pop->pos);
	free(pop->vel);
	free(pop->q);
	free(pop->m);
	free(pop->energy);
	free(pop);

}

void populateUniformly(const dictionary *ini, Population *pop, const Grid *grid, const gsl_rng *rng){

	// Read from ini
	int nSpecies, nDims;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",&nSpecies);
	int *nTGPoints = iniGetIntArr(ini,"grid:nTGPoints",&nDims);
	

	// Read from grid
	int *node = grid->node;
	int *nNodes = grid->nNodes;
	double *posToNode = grid->posToNode;

	// Compute normalized length of global reference frame
	int *L = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++){
		L[d] = nNodes[d]*nTGPoints[d]-1;
		msg(STATUS,"%i=L=nNodes*(nTGPoints-1)=%i*(%i-1)",L[d],nNodes[d],nTGPoints[d]);
	}
	msg(STATUS,"Domain size: %i,%i,%i",L[0],L[1],L[2]);
	msg(STATUS,"nDims=%i",nDims);

	for(int s=0;s<nSpecies;s++){

		msg(STATUS,"specie %i",s);

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart-1;
		double *pos = &pop->pos[iStart*nDims];

		// Iterate through all particles to be generated
		// Same seed on all MPI nodes ensure same particles are generated everywhere.
		for(long int i=0;i<nParticles[s];i++){	// Generate nParticles of this specie

			// Generate position for particle i
			for(int d=0;d<nDims;d++) pos[d] = L[d]*gsl_rng_uniform_pos(rng);

			// Count the number of dimensions where the particle resides in the range of this node
			int correctRange = 0;
			for(int d=0;d<nDims;d++) correctRange += (node[d] == (int)(posToNode[d]*pos[d]));
//			msg(STATUS,"pos of particle %i: %f,%f,%f",i,pos[0],pos[1],pos[2]);

			// Iterate only if particle resides in this sub-domain.
			if(correctRange==nDims){
				pos += nDims;
				iStop++;
//				msg(STATUS,"Keeping this");
			}

		}

		if(iStop>=pop->iStart[s+1]){
			int allocated = pop->iStart[s+1]-iStart;
			int generated = iStop-iStart+1;
			msg(ERROR,"allocated only %i particles of specie %i per node but %i generated",allocated,s,generated);
		}

		pop->iStop[s]=iStop;



	}

//	free(L);
//	free(nParticles);

}
