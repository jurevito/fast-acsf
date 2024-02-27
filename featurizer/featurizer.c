#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
    double x;
    double y;
    double z;
    int atom_index;
} Coord;

double eucl_dist(Coord p1, Coord p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

int square(int x) {
	return x*x;
}

void featurize(Coord* mol_atoms, int num_mol_atom, Coord* protein_atoms, int num_protein_atom) {
    int** dist = (int**)malloc(num_mol_atom*sizeof(int*));
    for (int i = 0; i < num_mol_atom; i++) {
        dist[i] = (int*)malloc(num_protein_atom*sizeof(int));
    }

    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
			dist[i][j] = eucl_dist(mol_atoms[i], protein_atoms[j]);
		}
    }

	printf("Matrix has %d element\n", num_mol_atom*num_protein_atom);
}