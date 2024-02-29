#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double RADIAL_CUTOFF = 12;
const double ANGULAR_CUTOFF = 6;

const double RADIAL_STEP = 0.5;
const double ANGULAR_STEP = 2.0;
const int N_THETA = 8;
const int N_ELEMENTS = 9;

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

double radial_sym_func(double Rij, double Rs) {
    const double Rc = RADIAL_CUTOFF;
    const double eta = 4.0;

    double fc = 0.5 * cos(M_PI * Rij / Rc) + 0.5;
    return exp(-eta * pow(Rij - Rs, 2)) * fc;
}

double angular_sym_func(double Rij, double Rik, double theta_ijk, double theta, double Rs) {
    double Rc = ANGULAR_CUTOFF;
    double zeta = 8;
    double eta = 4;

    double fc_Rij = 0.5 * cos(M_PI * Rij / Rc) + 0.5;
    double fc_Rik = 0.5 * cos(M_PI * Rik / Rc) + 0.5;
    return pow(2, 1-zeta) * pow(1 + cos(theta_ijk - theta), zeta) * exp(-eta * pow(((Rij + Rik) / 2 - Rs), 2)) * fc_Rij * fc_Rik;
}

void featurize(Coord* mol_atoms, int num_mol_atom, Coord* protein_atoms, int num_protein_atom) {

    // Setting radial steps.
    int rs_radial_length = ceil((RADIAL_CUTOFF - 0.5) / RADIAL_STEP);
    double* rs_radial = malloc(rs_radial_length*sizeof(double));

    double rad = 0.5;
    for (int i = 0; i<rs_radial_length; i++) {
        rs_radial[i] = rad;
        rad += RADIAL_STEP;
    }

    // Setting angular steps.
    int rs_angular_length = ceil((ANGULAR_CUTOFF - 0.5) / ANGULAR_STEP);
    double* rs_angular = malloc(rs_angular_length*sizeof(double));

    double ang = 0.5;
    for (int i = 0; i<rs_angular_length; i++) {
        rs_angular[i] = ang;
        ang += ANGULAR_STEP;
    }

    // Setting theta list.
    double* theta_list = malloc(N_THETA*sizeof(double));
    double interval = 2*M_PI / (N_THETA);
    for (int i = 0; i < N_THETA; i++) {
        theta_list[i] = interval*i;
    }

    //Setting PL distance matrix.
    double** dist = malloc(num_mol_atom*sizeof(double*));
    for (int i = 0; i < num_mol_atom; i++) {
        dist[i] = malloc(num_protein_atom*sizeof(double));
    }

    // Setting result vector.
    int radial_length = N_ELEMENTS*N_ELEMENTS * rs_radial_length;
    double* result = calloc(radial_length, sizeof(double));

    int smaller = 0;
    double sum = 0;

    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
			dist[i][j] = eucl_dist(mol_atoms[i], protein_atoms[j]);

            // Calculate features only for pairs within radius.
            if (dist[i][j] < RADIAL_CUTOFF) {
                
                // Calculate radial features.
                for (int k = 0; k<rs_radial_length; k++) {
                    int index = (rs_radial_length*N_ELEMENTS*mol_atoms[i].atom_index) + (rs_radial_length*protein_atoms[j].atom_index) + k;
                    result[index] += radial_sym_func(dist[i][j], rs_radial[k]);
                }
            }
		}
    }

    
    // printf("radial_length: %d\n", radial_length);
    //for (int i = 0 ; i<20 ; i++) {
    //    printf("result: %f\n", result[i]);
    //}
    //printf("Matrix has %d element (%d)\n", num_mol_atom*num_protein_atom, smaller);

    // Free all memory that was used.
    free(rs_radial);
    free(rs_angular);
    free(theta_list);
    for (int i = 0; i < num_mol_atom; i++) {
        free(dist[i]);
    }
    free(dist);
    free(result);
}