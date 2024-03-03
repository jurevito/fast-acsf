compile:
	gcc -shared -O2 -ffast-math -Wall -o .\acsf\featurizer.dll -fPIC .\acsf\featurizer.c .\acsf\util.c
