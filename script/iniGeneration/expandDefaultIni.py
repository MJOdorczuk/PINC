def overwrite(params, label, default):
    return params[label] if label in params else default

def generateFilesPart(files_params):
    # location of object.grid.h5 and the simulation output
    output = overwrite(files_params, 'output', 'data/')
    return f'''[files]
output                  = {output}
'''

# copied from previous examples - todo: check the purpose and describe
def generateMsgfilesPart(msgfiles_params):
    parsedumpFile = overwrite(msgfiles_params, "parsedump", "parsedump.txt")
    collisionFile = overwrite(msgfiles_params, "collision", "CollisionDump.txt")
    return f'''[msgfiles]
parsedump               = {parsedumpFile}
collision               = {collisionFile}
'''

def generateTimePart(time_params):
    # number of timesteps of the simulation
    timesteps=overwrite(time_params, "nTimeSteps", 3000)
    # length of the timestep in seconds
    timestep=overwrite(time_params, "timeStep", 2.5e-9)
    return f'''[time]
nTimeSteps              = {timesteps}
timeStep                = {timestep}
'''

def generateGridPart(grid_params):
    # number of dimensions of the grid
    nDims = overwrite(grid_params, "nDims", 3)
    # number of subdomains in each direction, this scirpts assumes 3D case
    nSubdomains = overwrite(grid_params, "nSubdomains", [1, 1, 1])
    # number of particles to allocate for (corner, edge, face) in particles per cell
    nEmigrantsAlloc = overwrite(grid_params, "nEmigrantsAlloc", [2, 4, 16])
    # size of the domain in each direction
    domainSize = overwrite(grid_params, "domainSize", [128, 96, 32])
    # step size in metres
    # roughly half the expected Debye length
    stepSize = overwrite(grid_params, "stepSize", 0.005)
    # number of ghost layers (assume same amount in all directions)
    nGhostLayers = overwrite(grid_params, "nGhostLayers", 1)
    # threshold for particle migration(todo: check the specific purpose)
    threshold = overwrite(grid_params, "threshold", 0.1)
    # boundary conditions for each direction [-x, +x, -y, +y, -z, +z]
    # NEUMANN means zero gradient boundaries - results in divergence (why?)
    # DIRICHLET means zero value boundaries... incorrect for downstream
    boundaries = overwrite(grid_params, "boundaries", ["DIRICHLET"] * 6)
    return f'''[grid]
nDims                   = {nDims}
nSubdomains             = {nSubdomains[0]},{nSubdomains[1]},{nSubdomains[2]}
nEmigrantsAlloc         = {nEmigrantsAlloc[0]} pc,{nEmigrantsAlloc[1]} pc,{nEmigrantsAlloc[2]} pc
trueSize                = {",".join(f"{int(domainSize[i]/nSubdomains[i])}" for i in range(3))}
stepSize                = {stepSize}
nGhostLayers            = {nGhostLayers}
thresholds              = {threshold}
boundaries              = {",".join(boundaries)}
'''

def generateFieldsPart(fields_params):
    # external magnetic field in Teslas(todo: check the unit)
    Bext = overwrite(fields_params, "BExt", [0, 0, 0])
    # external electric field in V/m(todo: check the unit)
    Eext = overwrite(fields_params, "EExt", [0, 0, 0])
    return f'''[fields]
BExt                    = {Bext[0]},{Bext[1]},{Bext[2]}
EExt                    = {Eext[0]},{Eext[1]},{Eext[2]}
'''

def generatePopulationPart(population_params):
    # number of species
    nSpecies = overwrite(population_params, "nSpecies", 2)
    # number of particles per cell(todo: check the unit)
    nParticles = overwrite(population_params, "nParticles", 16)
    # number of particles per cell to allocate memory for(todo: check the unit)
    nAlloc = overwrite(population_params, "nAlloc", 32)
    # charge of the particles in Coulombs
    # electrons and monoatomic monopositive oxygen
    charge = overwrite(population_params, "charge", [-1.60217663e-19, 1.60217663e-19])
    # mass of the particles in kilograms
    # electrons and monoatomic monopositive oxygen
    mass = overwrite(population_params, "mass", [9.1093837e-31, 2.6566962e-26])
    # density of the particles in 1/m^3
    density = overwrite(population_params, "density", [1e11, 1e11])
    # drift of the particles in electron thermal velocities
    drift = overwrite(population_params, "drift", [0] * 6)
    # Deprecated
    # # amplitude of the perturbation of the particles in m(todo: check the unit and specific purpose)
    # perturbAmplitude = overwrite(population_params, "perturbAmplitude", [1e-5, 0, 0, 0, 0, 0])
    # # mode of the perturbation of the particles(todo: check the specific purpose)
    # perturbMode = overwrite(population_params, "perturbMode", [1, 0, 0, 0, 0, 0])
    # thermal velocity of the particles in m/s
    # 188 km/s for electrons, 750 m/s for O^+
    thermalVelocity = overwrite(population_params, "thermalVelocity", [1.88e5, 750])
    # maximum velocity of the particles in Courant numbers
    # (todo: should I not make it slightly lower than 1?)
    maxVel = overwrite(population_params, "maxVel", 1)
    return f'''[population]
nSpecies                = {nSpecies}
nParticles              = {nParticles} pc
nAlloc                  = {nAlloc} pc
charge                  = {charge[0]},{charge[1]}
mass                    = {mass[0]},{mass[1]}
density                 = {density[0]},{density[1]}
drift                   = {",".join(map(str, drift))}
thermalVelocity         = {",".join(map(str, thermalVelocity))}
maxVel                  = {maxVel}
'''
# perturbAmplitude        = {",".join(map(str, perturbAmplitude))}
# perturbMode             = {",".join(map(str, perturbMode))}

def generateMethodsPart(methods_params):
    # collider method
    collisionMode = overwrite(methods_params, "mode", "oCollMode")
    # normalisation method
    normalization = overwrite(methods_params, "normalization", "SI")
    # Poisson solver method(todo: check the specific purpose)
    poisson = overwrite(methods_params, "poisson", "mgSolver")
    # accelerator method
    acc = overwrite(methods_params, "acc", "puAcc3D1KE")
    # distribution method
    distr = overwrite(methods_params, "distr", "puDistr3D1split")
    # migration method
    migrate = overwrite(methods_params, "migrate", "puExtractEmigrants3DOpen")
    # Monte Carlo Collision method(todo: verify the options)
    mcc = overwrite(methods_params, "mcc", "mccConstCrossect")
    return f'''[methods]
mode                    = {collisionMode}
normalization           = {normalization}
poisson                 = {poisson}
acc                     = {acc}
distr                   = {distr}
migrate                 = {migrate}
mcc                     = {mcc}
'''

def generateMultigridPart(multigrid_params):
    # Choice of multigrid cycle type
    mgCycle = overwrite(multigrid_params, "mgCycle", "mgVRecursive")
    # Choice of presmoother method
    preSmooth = overwrite(multigrid_params, "preSmooth", "gaussSeidelRBND")
    # Choice of postsmoother method
    postSmooth = overwrite(multigrid_params, "postSmooth", "gaussSeidelRBND")
    # Choice of coarse grid solver
    coarseSolver = overwrite(multigrid_params, "coarseSolver", "gaussSeidelRBND")
    # Number of multigrid levels
    mgLevels = overwrite(multigrid_params, "mgLevels", 4)
    # Number of multigrid cycles
    mgCycles = overwrite(multigrid_params, "mgCycles", 5)
    # Number of iterations for the presmoother
    nPreSmooth = overwrite(multigrid_params, "nPreSmooth", 10)
    # Number of iterations for the postsmoother
    nPostSmooth = overwrite(multigrid_params, "nPostSmooth", 10)
    # Number of iterations for the coarse grid solver
    nCoarseSolve = overwrite(multigrid_params, "nCoarseSolve", 10)
    # Choice of prolongation stencil(todo: read about it)
    prolongator = overwrite(multigrid_params, "prolongator", "bilinear")
    # Choice of restrictor stencil(todo: read about it)
    restrictor = overwrite(multigrid_params, "restrictor", "halfWeight")
    # Run number(todo: check the specific purpose)
    runNumber = overwrite(multigrid_params, "runNumber", 0.0)
    # Tolerance for the multigrid solver
    tol = overwrite(multigrid_params, "tol", 1e-6)
    # Objective tolerance for the multigrid solver
    objTol = overwrite(multigrid_params, "objTol", 1e-9)
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

def generateObjectPart(object_params):
    # Bias included
    biasOn = overwrite(object_params, "biasOn", 0)
    # Bias value in Volts
    bias = overwrite(object_params, "bias", 0)
    # Sweep included
    sweepOn = overwrite(object_params, "sweepOn", 0)
    # Sweep time in seconds(todo: check the unit, comment stated hertz but it does not make sense)
    sweepTime = overwrite(object_params, "sweepTime", 25e-8)
    # Sweep range in Volts (only max)
    sweepRange = overwrite(object_params, "sweepRange", 6)
    # Sweep offset in Volts
    sweepOffset = overwrite(object_params, "sweepOffset", 2)
    # Sweep start time in seconds
    sweepStart = overwrite(object_params, "sweepStart", 1.2e-8)
    # Sweep steps
    sweepSteps = overwrite(object_params, "sweepSteps", 50)
    return f'''[object]
biasOn                  = {biasOn}
bias                    = {bias}
sweepOn                 = {sweepOn}
sweepTime               = {sweepTime}
sweepRange              = {sweepRange}
sweepOffset             = {sweepOffset}
sweepStart              = {sweepStart}
sweepSteps              = {sweepSteps}
'''

def generateCollisionsPart(collisions_params):
    # Choice of energy transfer in electron collisions(todo: verify)
    electronEnergyMethod = overwrite(collisions_params, "electronEnergyMethod", "conservative")
    # Number of neutral species
    nSpeciesNeutral = overwrite(collisions_params, "nSpeciesNeutral", 1)
    # Mass of the neutral species in kg (oxygen)
    neutralMass = overwrite(collisions_params, "neutralMass", 2.6566962e-26)
    # Drift of the neutral species in electron thermal velocities
    neutralDrift = overwrite(collisions_params, "neutralDrift", [0, 0, 0])
    # Number density of the neutral species in 1/m^3
    numberDensityNeutrals = overwrite(collisions_params, "numberDensityNeutrals", 1e15)
    # Thermal velocity of the neutral species in electron thermal velocities
    thermalVelocityNeutrals = overwrite(collisions_params, "thermalVelocityNeutrals", 750)
    # Artificial loss factor for electrons(todo: check the specific purpose)
    artificialLoss = overwrite(collisions_params, "artificialLoss", 1.0)
    # Collision frequency for charge exchange collisions
    collFrqCex = overwrite(collisions_params, "collFrqCex", 1.45)
    # Collision frequency for ion elastic collisions
    collFreIonElastic = overwrite(collisions_params, "collFrqIonElastic", 3.674)
    # Collision frequency for electron elastic collisions
    collFrqElectronElastic = overwrite(collisions_params, "collFrqElectronElastic", 12.48)
    # Charge exchange collisional cross section
    sigmaCEX = overwrite(collisions_params, "sigmaCEX", 5e-19)
    # Ion elastic collisional cross section
    # Normalised by multiplying by Debye length and ion number density
    sigmaIonElastic = overwrite(collisions_params, "sigmaIonElastic", 5.5e-19)
    # Electron elastic collisional cross section
    # Normalised by multiplying by Debye length and ion number density
    sigmaElectronElastic = overwrite(collisions_params, "sigmaElectronElastic", 4e-20)
    # If using functional form of crossections e.g mccGetPmax...() freq propto v,
    # this is experimental and overwrites above collfreqs, if using method "mccFunctionalCrossect".
    # sigma_adj = a * exp(- b * v^2)
    # a parameter given in m^2?
    # b parameter given in s^2/m^2?
    CEX_a = overwrite(collisions_params, "CEX_a", 0.00241)
    CEX_b = overwrite(collisions_params, "CEX_b", 57.06)
    ion_elastic_a = overwrite(collisions_params, "ion_elastic_a", 0.00081)
    ion_elastic_b = overwrite(collisions_params, "ion_elastic_b", 150.0)
    electron_a = overwrite(collisions_params, "electron_a", 0.001205)
    electron_b = overwrite(collisions_params, "electron_b", 1.2758)
    return f'''[collisions]
electronEnergyMethod    = {electronEnergyMethod}
nSpeciesNeutral         = {nSpeciesNeutral}
neutralMass             = {neutralMass}
neutralDrift            = {",".join(map(str, neutralDrift))}
numberDensityNeutrals   = {numberDensityNeutrals}
thermalVelocityNeutrals = {thermalVelocityNeutrals}
artificialLoss          = {artificialLoss}
collFrqCex              = {collFrqCex}
collFrqIonElastic       = {collFreIonElastic}
collFrqElectronElastic  = {collFrqElectronElastic}
sigmaCEX                = {sigmaCEX}
sigmaIonElastic         = {sigmaIonElastic}
sigmaElectronElastic    = {sigmaElectronElastic}
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
    import json

    drift = [7800, 0, 0]
    params = "{}"

    if len(sys.argv) > 4:
        with h5py.File(sys.argv[1], "r") as objectFile:
            z, y, x = objectFile["Object"][()].shape
            nx, ny, nz = sys.argv[2:5]
    if len(sys.argv) > 5:
        params = sys.argv[5]
    params = json.loads(params)
    if "grid" not in params:
        params["grid"] = {}
    params["grid"]["domainSize"] = [int(x), int(y), int(z)]
    params["grid"]["nSubdomains"] = [int(nx), int(ny), int(nz)]

    iniContent = f'''
; @file			input.ini
; @brief		PINC input file template.
; @author		Michał Jan Odorczuk <michaljo@uio.no>
;
{generateFilesPart(overwrite(params, "files", {}))}
{generateMsgfilesPart(overwrite(params, "msgfiles", {}))}
{generateTimePart(overwrite(params, "time", {}))}
{generateGridPart(params["grid"])}
{generateFieldsPart(overwrite(params, "fields", {}))}
{generatePopulationPart(overwrite(params, "population", {}))}
{generateMethodsPart(overwrite(params, "methods", {}))}
{generateMultigridPart(overwrite(params, "multigrid", {}))}
{generateObjectPart(overwrite(params, "object", {}))}
{generateCollisionsPart(overwrite(params, "collisions", {}))}
'''
    with open("input.ini", "w") as f:
        f.write(iniContent)