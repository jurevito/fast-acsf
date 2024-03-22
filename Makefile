compile_windows:
	gcc -shared -O2 -ffast-math -fopenmp -static -Wall -o ./acsf/featurizer.dll -fPIC ./acsf/featurizer.c ./acsf/util.c

compile_unix:
	gcc -shared -O2 -ffast-math -fopenmp -Wall -o ./acsf/featurizer.so -fPIC ./acsf/featurizer.c ./acsf/util.c
