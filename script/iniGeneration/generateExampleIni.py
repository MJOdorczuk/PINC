def generateFilesPart(
        # location of object.grid.h5 and the simulation output
        output="data/"
        ):
    return f'''[files]
output                  = {output}
'''

# copied from previous examples - todo: check the purpose and describe
def generateMsgfilesPart(
        parsedumpFile="parsedump.txt",
        collisionFile="CollisionDump.txt"):
    return f'''[msgfiles]
parsedump               = {parsedumpFile}
collision               = {collisionFile}
'''

def generateTimePart(
        # number of timesteps of the simulation
        timesteps=3000,
        # lenfth of the timestep in seconds(todo: check the unit)
        timestep=1e-8):
    return f'''[time]
nTimeSteps              = {timesteps}
timeStep                = {timestep}
'''

def generateGridPart(
        # number of subdomains in each direction, this scirpts assumes 3D case
        nSubdomains = [2, 2, 2],
        # number of particles to allocate for (corner, edge, face) in particles per cell(todo: check the unit)
        nEmigrantsAlloc = [2, 4, 16],
        # size of the domain in each direction
        domainSize = [32, 32, 32],
        # step size in Debye lengths of specie 0(todo: check the unit)
        # Is it really Debye lengths? It seems to be the size of the cell in meters
        stepSize = 0.0142,
        # number of ghost layers (assume same amount in all directions)
        nGhostLayers = 1,
        # threshold for particle migration(todo: check the specific purpose)
        threshold = 0.1,
        # boundary conditions for each direction [-x, +x, -y, +y, -z, +z](todo: check correctness)
        boundaries = "DIRICHLET,DIRICHLET,DIRICHLET,DIRICHLET,DIRICHLET,DIRICHLET"
        ):
    return f'''[grid]
nDims                   = 3
nSubdomains             = {nSubdomains[0]},{nSubdomains[1]},{nSubdomains[2]}
nEmigrantsAlloc         = {nEmigrantsAlloc[0]} pc,{nEmigrantsAlloc[1]} pc,{nEmigrantsAlloc[2]} pc
trueSize                = {",".join(f"{int(domainSize[i]/nSubdomains[i])}" for i in range(3))}
stepSize                = {stepSize}
nGhostLayers            = {nGhostLayers}
thresholds              = {threshold}
boundaries              = {boundaries}
'''

def generateFieldsPart(
        # external magnetic field in Teslas(todo: check the unit)
        Bext = [0, 0, 0],
        # external electric field in V/m(todo: check the unit)
        Eext = [0, 0, 0]
        ):
    return f'''[fields]
BExt                    = {Bext[0]},{Bext[1]},{Bext[2]}
EExt                    = {Eext[0]},{Eext[1]},{Eext[2]}
'''

def generatePopulationPart(
        # number of particles per cell(todo: check the unit)
        nParticles = 16,
        # number of particles per cell to allocate memory for(todo: check the unit)
        nAlloc = 32,
        # charge of the particles in Coulombs(todo: check the unit)
        charge = [-1.60217662e-19, 1.60217662e-19],
        # mass of the particles in kg(todo: check the unit)
        mass = [9.10938356e-31, 1.6726219e-27],
        # density of the particles in m^-3(todo: check the unit and specific purpose)
        density = [5.8977e9, 5.8977e9],
        # drift of the particles in m/s(todo: check the unit)
        drift = [0, 0, 0, 0, 0, 0],
        # amplitude of the perturbation of the particles in m(todo: check the unit and specific purpose)
        perturbAmplitude = [1e-5, 0, 0, 0, 0, 0],
        # mode of the perturbation of the particles(todo: check the specific purpose)
        perturbMode = [1, 0, 0, 0, 0, 0],
        # thermal velocity of the particles in m/s(todo: check the unit)
        thermalVelocity = [123111, 2873],
        # maximum velocity of the particles in m/s(todo: check the unit, it cannot be true)
        maxVel = 1
        ):
    return f'''[population]
nSpecies                = 2
nParticles              = {nParticles} pc
nAlloc                  = {nAlloc} pc
charge                  = {charge[0]},{charge[1]}
mass                    = {mass[0]},{mass[1]}
density                 = {density[0]},{density[1]}
drift                   = {",".join(map(str, drift))}
perturbAmplitude        = {",".join(map(str, perturbAmplitude))}
perturbMode             = {",".join(map(str, perturbMode))}
thermalVelocity         = {",".join(map(str, thermalVelocity))}
maxVel                  = {maxVel}
'''

def generateMethodsPart(
        # collider method
        collisionMode = "oCollMode",
        # Poisson solver method(todo: check the specific purpose)
        poisson = "mgSolver",
        # accelerator method
        acc = "puBoris3D1KETEST",
        # distribution method
        distr = "puDistr3D1split",
        # migration method
        migrate = "puExtractEmigrants3DOpen",
        # Monte Carlo Collision method(todo: verify the options)
        mcc = "mccConstCrossect"
        ):
    return f'''[methods]
mode                    = {collisionMode}
normalization           = SI
poisson                 = {poisson}
acc                     = {acc}
distr                   = {distr}
migrate                 = {migrate}
mcc                     = {mcc}
'''

def generateMultigridPart(
        # Choice of multigrid cycle type
        mgCycle = "mgVRecursive",
        # Choice of presmoother method
        preSmooth = "gaussSeidelRBND",
        # Choice of postsmoother method
        postSmooth = "gaussSeidelRBND",
        # Choice of coarse grid solver
        coarseSolver = "gaussSeidelRBND",
        # Number of multigrid levels
        mgLevels = 4,
        # Number of multigrid cycles
        mgCycles = 5,
        # Number of iterations for the presmoother
        nPreSmooth = 10,
        # Number of iterations for the postsmoother
        nPostSmooth = 10,
        # Number of iterations for the coarse grid solver
        nCoarseSolve = 10,
        # Choice of prolongation stencil(todo: read about it)
        prolongator = "bilinear",
        # Choice of restrictor stencil(todo: read about it)
        restrictor = "halfWeight",
        # Run number(todo: check the specific purpose)
        runNumber = 0.0,
        # Tolerance for the multigrid solver
        tol = 1e-6,
        # Objective tolerance for the multigrid solver
        objTol = 1e-9
        ):
    return f'''[multigrid]
cycle                   = {mgCycle}
preSmooth               = {preSmooth}
postSmooth              = {postSmooth}
coarseSolver            = {coarseSolver}
mgLevels                = {mgLevels}
mgCycles                = {mgCycles}
nPreSmooth              = {nPreSmooth}
nPostSmooth             = {nPostSmooth}
nCoarseSolve            = {nCoarseSolve}
prolongator             = {prolongator}
restrictor              = {restrictor}
runNumber		        = {runNumber}
tol 		            = {tol}
objTol		            = {objTol}
'''

def generateObjectPart(
        # Bias included
        biasOn = False,
        # Bias value in Volts
        bias = 2.0,
        # Sweep included
        sweepOn = False,
        # Sweep time in seconds(todo: check the unit, comment stated hertz but it does not make sense)
        sweepTime = 25e-8,
        # Sweep range in Volts (only max)
        sweepRange = 6,
        # Sweep offset in Volts
        sweepOffset = 2,
        # Sweep start time in seconds
        sweepStart = 1.2e-8,
        # Sweep steps
        sweepSteps = 50
        ):
    return f'''[object]
biasOn                  = {1 if biasOn else 0}
bias                    = {bias}
sweepOn                 = {1 if sweepOn else 0}
sweepTime               = {sweepTime}
sweepRange              = {sweepRange}
sweepOffset             = {sweepOffset}
sweepStart              = {sweepStart}
sweepSteps              = {sweepSteps}
'''

def generateCollisionsPart(
        # Choice of energy transfer in electron collisions(todo: verify)
        electronEnergyMethod = "conservative",
        # Number of neutral species
        nSpeciesNeutral = 1,
        # Mass of the neutral species in kg(todo: check the unit)
        neutralMass = 1.6724828e-27,
        # Drift of the neutral species in m/s(todo: check the unit)
        neutralDrift = [0, 0, 0],
        # Number density of the neutral species in m^-3(todo: check the unit)
        numberDensityNeutrals = 5.8977e9,
        # Thermal velocity of the neutral species in m/s(todo: check the unit)
        thermalVelocityNeutrals = 2873,
        # Artificial loss factor for electrons(todo: check the specific purpose)
        artificialLoss = 1.0,
        # Real electron mass in kg(todo: check the unit and specific purpose)
        realElectronMass = 9.10938356e-31,
        # Collision frequency for charge exchange collisions
        collFrqCex = 1.44e8,
        # Collision frequency for ion elastic collisions
        collFreIonElastic = 2.89e8,
        # Collision frequency for electron elastic collisions
        collFrqElectronElastic = 3.27e9,
        # If using functional form of crossections e.g mccGetPmax...() freq propto v, this is experimental and overwrites above collfreqs, if using method "mccFunctionalCrossect".
        CEX_a = 0.00241,
        CEX_b = 57.06,
        ion_elastic_a = 0.00081,
        ion_elastic_b = 150.0,
        electron_a = 0.001205,
        electron_b = 1.2758
        ):
    return f'''[collisions]
electronEnergyMethod    = {electronEnergyMethod}
nSpeciesNeutral         = {nSpeciesNeutral}
neutralMass             = {neutralMass}
neutralDrift            = {",".join(map(str, neutralDrift))}
numberDensityNeutrals   = {numberDensityNeutrals}
thermalVelocityNeutrals = {thermalVelocityNeutrals}
artificialLoss          = {artificialLoss}
realElectronMass        = {realElectronMass}
collFrqCex              = {collFrqCex}
collFrqIonElastic       = {collFreIonElastic}
collFrqElectronElastic  = {collFrqElectronElastic}
CEX_a                   = {CEX_a}
CEX_b                   = {CEX_b}
ion_elastic_a           = {ion_elastic_a}
ion_elastic_b           = {ion_elastic_b}
electron_a              = {electron_a}
electron_b              = {electron_b}
'''


if __name__ == "__main__":
    import h5py
    import sys

    drift = [11492, 0, 0]

    if len(sys.argv) > 4:
        with h5py.File(sys.argv[1], "r") as objectFile:
            z, y, x = objectFile["Object"][()].shape
            nx, ny, nz = sys.argv[2:5]

    iniContent = f'''
; @file			input.ini
; @brief		PINC input file template.
; @author		Michał Jan Odorczuk <michaljo@uio.no>
;
{generateFilesPart()}
{generateMsgfilesPart()}
{generateTimePart()}
{generateGridPart(nSubdomains=[int(nx), int(ny), int(nz)], domainSize=[x, y, z])}
{generateFieldsPart()}
{generatePopulationPart(drift=drift+drift)}
{generateMethodsPart()}
{generateMultigridPart()}
{generateObjectPart()}
{generateCollisionsPart(neutralDrift=drift)}
'''
    with open("input.ini", "w") as f:
        f.write(iniContent)