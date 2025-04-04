make veryclean
make
rm data/* 2> /dev/null
# cp test_object_64_32_32.grid.h5 data/object.grid.h5
# Dimensions in x, y and z, length in x-axis, front location in x-axis
# 0 for just the main body, 1 for the body and the probes
python3 script/ObjectCreation/NorSat.py 128 96 32 10 20 0 data/object.grid.h5
python3 script/iniGeneration/expandDefaultIni.py data/object.grid.h5 2 3 1 \
'{  "time" : { "nTimeSteps" : 40000},
    "grid" : { "nEmigrantsAlloc" : [8, 16, 64], "nGhostLayers" : 1 },
    "population" : { "nParticles" : 16, "nAlloc" : 32, "drift" : [7725, 0, 0, 7725, 0, 0]},
    "methods" : { "mode" : "oCollCustomRhoMode", "mcc" : "mccConstFreq" },
    "collisions" : { "neutralDrift" : [7725, 0, 0], "numberDensityNeutrals" : 1e16, "collFrqCex" : 14.5, "collFrqIonElastic" : 36.74, "collFrqElectronElastic" : 124.8 }}'
python3 script/fieldsInitialisation/uniformFields.py 128 96 32 data \
'{"rhoNeutral": 1, "v0Neutral": 1, "v1Neutral": 0, "v2Neutral": 0, "vthNeutral": 1}'
# The values should be given via .ini file
# TODO: change the uniform fields to be of value 1, I guess
# '{"rhoNeutral": 5.8977e9, "v0Neutral": 7800, "v1Neutral": 0, "v2Neutral": 0, "vthNeutral": 2873}'
echo Running PINC
mpirun -np 6 pinc input.ini > out.log 2> err.log
echo simulation done