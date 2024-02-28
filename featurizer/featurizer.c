#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double RADIAL_CUTOFF = 12;
const double ANGULAR_CUTOFF = 6;

const double RADIAL_STEP = 0.5;
const double ANGULAR_STEP = 2.0;
const int N_THETA = 8;

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

    // Setting radial steps.
    int rs_radial_length = ceil((RADIAL_CUTOFF - 0.5) / RADIAL_STEP);
    double* rs_radial = (double*)malloc(rs_radial_length*sizeof(double));

    double rad = 0.5;
    for (int i = 0; i<rs_radial_length; i++) {
        rs_radial[i] = rad;
        rad += RADIAL_STEP;
    }

    // Setting angular steps.
    int rs_angular_length = ceil((ANGULAR_CUTOFF - 0.5) / ANGULAR_STEP);
    double* rs_angular = (double*)malloc(rs_angular_length*sizeof(double));

    double ang = 0.5;
    for (int i = 0; i<rs_angular_length; i++) {
        rs_angular[i] = ang;
        ang += ANGULAR_STEP;
    }

    // Setting theta list.
    double* theta_list = (double*)malloc(N_THETA*sizeof(double));
    double interval = 2*M_PI / (N_THETA);
    for (int i = 0; i < N_THETA; i++) {
        theta_list[i] = interval*i;
    }


    double** dist = (double**)malloc(num_mol_atom*sizeof(double*));
    for (int i = 0; i < num_mol_atom; i++) {
        dist[i] = (double*)malloc(num_protein_atom*sizeof(double));
    }

    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
			dist[i][j] = eucl_dist(mol_atoms[i], protein_atoms[j]);
		}
        break;
    }

	printf("Matrix has %d element\n", num_mol_atom*num_protein_atom);
}