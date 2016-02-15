#ifndef MULTIGRID_H
#define MULTIGRID_H

/**
 * @file		multigrid.h
 * @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Poisson Solver, multigrid.
 * @date		26.10.15
 *
 *
 * Functions dealing with the initialisation and destruction of multigrid structures and
 * a multigrid solver containing restriction, prolongation operatorors and smoothers
 *
 */


/**
  * @brief Contains the grids needed in a multigrid solver, as well as other specifications
  * @param ini 			Input file, contains specification for the run
  * @param gridQuantity  Grid with quantities and grid specifications as a memmber
  *
  * The finest grid, grid 0, links to the grid used in the rest of the program, the other grids
  * in the grid array is coarser grids used in the
  * The preSmooth, postSmooth and coarseSolv is set in the input.ini file, for now the only
  * options are Gauss-Seidel Red-Black (gaussSeidel). More TBD.
  */
 typedef struct {
    Grid **grids;   ///< Array of Grid structs of decreasing coarseness
    int nLevels;         			///< #Grid levels
    int nMGCycles;         			///< Multigrid cycles we want to run
	int nPreSmooth;					///<
	int nPostSmooth;
	int nCoarseSolve;

    void (*coarseSolv)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Coarse Grid Solver function
    void (*postSmooth)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Post Smooth function
    void (*preSmooth)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Pre Smooth function
	void (*restrictor)(const Grid *fine, Grid *coarse);	///< Function pointer to restrictor
	void (*prolongator)(Grid *fine, const Grid *coarse, const MpiInfo *mpiInfo);	///< Function pointer to prolongator
} Multigrid;


/**
 * @brief Allocates multigrid struct
 * @param grid 		Finest grid
 * @param ini		Dictionary
 *
 *	Allocates memory for a multgrid struct. A multigrid struct consist of an array of
 *	grids, where each subsequent grid has half the true grid points of the previous.
 *
 *	In addition it stores a few additional parameters that are useful during the
 *	handling of the multigrid algorithm.
 * 	nLevels:	Number of grids
 *	nMGCycles:	Number of MG cycles to run
 *	nPreSmooth:	Number of cycles the presmoother to run
 *	nPostSmooth:Number of cycles the postsmoother to run
 *	nCoarseSolve:Number of cycles for the coarse solver to run
 *
 *	The algorithms for the solver, restrictors and prolongators are set in the allocation according
 *	to a input file, then it is handled by a function pointer stored in the MG struct.
 *'	So if we want to use the prolongator, the function pointer can just be used as a regular function
 * 	and the chosen prolongator algorithm will be used.
 * 	\code
		multigrid->prolongator(fine, coarse, mpiInfo);
 *	\endcode
 *
 *
 *	NB!The number of true grid points used in the finest grid needs to be a multiple
 *	nLevels*2, to make it possible to half the grid points down to the coarsest grid.
 */

Multigrid *mgAlloc(const dictionary *ini, Grid *grid);

 /**
  * @brief Free multigrid struct, top gridQuantity needs to be freed seperately
  * @param 	multigrid
  *
  * Since the finest grid is allocated seperately and used on it's own without
  * the multigrid struct, it is not freed in this destructor.
  * Variables freed: gridQuantity [1->end]
  *
  */
void mgFree(Multigrid *multigrid);

/**
 * @brief Solves Poissons equation for electric potential, with multigrid V cycles
 * @param	mrRho	Source term
 * @param	mgPhi	Solution term
 * @param	mgRes	Residual
 * @param	mpiInfo	Subdomain information
 * @return	mgPhi
 *
 *	This is an implementation of a Multigrid V Cycle solver. See "DOC" for more information.
 */

void linearMGSolv(Multigrid *mgRho, Multigrid *mgPhi, Multigrid *mgRes, const MpiInfo *mpiInfo);

/**
 * @brief Gauss-Seidel Red and Black 3D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	3D dimensional implementation of Gauss-Seidel RB, which does several sweeps through the
 *	grid trying to simplify the iteration through the grid.
 *
 */
void gaussSeidel3DNew(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief Gauss-Seidel Red and Black 3D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	3D dimensional implementation of Gauss-Seidel RB, which does one sweep through the grid for
 *	each color, but has slightly more complicated behaviour on the edges, due to needing to skip
 *	the ghostlayers.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void gaussSeidel3D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief Gauss-Seidel Red and Black 2D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	2D dimensional implementation of Gauss-Seidel RB, which does several sweeps through the
 *	grid trying to simplify the iteration through the grid.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void gaussSeidel2D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief Jacobian method
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 * Non-optimized implementation of a jacobian algorithm to solve poissons equation.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void jacobian(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);


/**
 * @brief Half weight restriction, 2D
 * @param	fine	Source term
 * @param	coarse	Solution term
 * @param	mpiInfo	Subdomain information
 * @return	coarse
 *
 *	Implementation of a half weight restriction scheme, copying the fine
 *	down to the coarser grid.
 *
 */

void halfWeightRestrict2D(const Grid *fine, Grid *coarse);


/**
 * @brief Half weight restriction, 3D
 * @param	fine	Source term
 * @param	coarse	Solution term
 * @return	coarse
 *
 *	Implementation of a half weight restriction scheme, copying the fine
 *	down to the coarser grid.
 *
 */
void halfWeightRestrict3D(const Grid *fine, Grid *coarse);

/**
 * @brief Bilinear interpolation, 2D
 * @param	fine	Fine grid
 * @param	coarse	Coarse grid
 * @param	mpiInfo	Subdomain information
 * @return	fine
 *
 *	Implementation of a bilnear interpolation scheme, interpolating the coarse
 *	grid onto the fine grid.
 *
 */

void bilinearProlong2D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);

/**
 * @brief Bilinear interpolation, 3D
 * @param	fine	Fine grid
 * @param	coarse	Coarse grid
 * @param	mpiInfo	Subdomain information
 * @return	fine
 *
 *	Implementation of a bilnear interpolation scheme, interpolating the coarse
 *	grid onto the fine grid.
 *
 */
void bilinearProlong3D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);

/**
 * @brief Computes residual
 * @param	res		Residual grid
 * @param	phi		Phi	grid
 * @param	rho		Rho grid
 * @param	mpiInfo	Subdomain information
 * @return	res
 *
 *	Computes the residual on a grid level.
 *	\f[
 *		d_l = \nabla^2_l\phi_l - \rho_l
 *	\f]
 */
void mgResidual(Grid *res, const Grid *rho, const Grid *phi,const MpiInfo *mpiInfo);

/**
 * @brief Returns mass of a grid
 * @param	grid		Grid struct
 * @param	ini			Dictionary
 *
 *	Computes a mass of a grid:
 *	\f[
 *		M = \frac{1}{N} \sum_N d^2_N
 *	\f]
 *
 * NB! Not written for efficiency
 */
 double mgResMass3D(Grid *grid);


 #endif // MULTIGRID_H
