make veryclean
make
rm data/* 2> /dev/null
# cp test_object_64_32_32.grid.h5 data/object.grid.h5
python3 script/ObjectCreation/NorSat.py 32 32 32 10 5 0 data/object.grid.h5
python3 script/iniGeneration/generateExampleIni.py data/object.grid.h5 2 2 2
python3 script/fieldsInitialisation/uniformFields.py 32 32 32 data \
'{"rhoNeutral": 1, "v0Neutral": 1, "v1Neutral": 0, "v2Neutral": 0, "vthNeutral": 1}'
# The values should be given via .ini file
# TODO: change the uniform fields to be of value 1, I guess
# '{"rhoNeutral": 5.8977e9, "v0Neutral": 7800, "v1Neutral": 0, "v2Neutral": 0, "vthNeutral": 2873}'
echo Running PINC
mpirun -np 8 pinc input.ini > out.log 2> err.log
echo simulation done