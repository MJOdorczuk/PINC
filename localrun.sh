make veryclean
make
rm data/* 2> /dev/null
# cp test_object_64_32_32.grid.h5 data/object.grid.h5
# Dimensions in x, y and z, length in x-axis, front location in x-axis
# 0 for just the main body, 1 for the body and the probes, 2 for a sphere
# 3 for the main body spanning almost all the width
python3 script/ObjectCreation/NorSat.py 128 96 32 10 20 0 data/object.grid.h5
python3 script/iniGeneration/expandDefaultIni.py data/object.grid.h5 2 3 1 \
'{  "time" : { "nTimeSteps" : 100000},
    "grid" : { "nEmigrantsAlloc" : [8, 16, 64], "nGhostLayers" : 1 },
    "population" : { "nParticles" : 16, "nAlloc" : 32, "drift" : [7725, 0, 0, 7725, 0, 0]},
    "methods" : { "mode" : "oCollCustomRhoMode", "mcc" : "mccConstCrossect" },
    "multigrid" : { "mgLevels" : 4 },
    "collisions" : { "neutralDrift" : [7725, 0, 0], "numberDensityNeutrals" : 1e19 }}'
# TODO: put in a field initiatilastion script
cp NorSatflow_rhoNeutral.grid.h5 data/rhoNeutral.grid.h5
cp NorSatFlow_v0Neutral.grid.h5 data/v0Neutral.grid.h5
cp NorSatFlow_v1Neutral.grid.h5 data/v1Neutral.grid.h5
cp NorSatFlow_v2Neutral.grid.h5 data/v2Neutral.grid.h5
cp NorSatFlow_vthNeutral.grid.h5 data/vthNeutral.grid.h5
# Uniform field case
# python3 script/fieldsInitialisation/uniformFields.py 128 96 32 data \
# '{"rhoNeutral": 1, "v0Neutral": 1, "v1Neutral": 0, "v2Neutral": 0, "vthNeutral": 1}'
echo Running PINC
mpirun -np 6 pinc input.ini > out.log 2> err.log
echo simulation done
