rm data/*
cp test_object_64_32_32.grid.h5 data/object.grid.h5
python3 script/iniGeneration/generateExampleIni.py data/object.grid.h5 2 2 2
echo Running PINC
mpirun -np 8 pinc input.ini