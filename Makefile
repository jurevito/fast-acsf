compile:
	gcc -shared -O2 -fopenmp -static -Wall -o .\acsf\featurizer.dll -fPIC .\acsf\featurizer.c .\acsf\util.c

