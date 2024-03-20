compile_windows:
	gcc -shared -O2 -ffast-math -Wall -o ./acsf/featurizer.dll -fPIC ./acsf/featurizer.c ./acsf/util.c

compile_unix:
	gcc -shared -O2 -ffast-math -Wall -o ./acsf/featurizer.so -fPIC ./acsf/featurizer.c ./acsf/util.c
