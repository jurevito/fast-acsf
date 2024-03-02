#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"

const double RADIAL_CUTOFF = 12;
const double ANGULAR_CUTOFF = 6;

const double RADIAL_STEP = 0.5;
const double ANGULAR_STEP = 2.0;
const int N_THETA = 8;
const int N_ELEMENTS = 9;

double radial_sym_func(double Rij, double Rs) {
    const double Rc = RADIAL_CUTOFF;
    const double eta = 4.0;

    double fc = 0.5 * cos(M_PI * Rij / Rc) + 0.5;
    return exp(-eta * pow(Rij - Rs, 2)) * fc;
}

double angular_sym_func(double Rij, double Rik, double theta_ijk, double theta, double Rs) {
    const double Rc = ANGULAR_CUTOFF;
    const double zeta = 8;
    const double eta = 4;

    double fc_Rij = 0.5 * cos(M_PI * Rij / Rc) + 0.5;
    double fc_Rik = 0.5 * cos(M_PI * Rik / Rc) + 0.5;
    return pow(2, 1-zeta) * pow(1 + cos(theta_ijk - theta), zeta) * exp(-eta * pow(((Rij + Rik) / 2 - Rs), 2)) * fc_Rij * fc_Rik;
}

int find_angular_index(int lig_atom_idx, int p1_atom_idx, int p2_atom_idx, int th, int rs, int rs_angular_length) {
    if (p1_atom_idx > p2_atom_idx) {
        int tmp = p1_atom_idx;
        p1_atom_idx = p2_atom_idx;
        p2_atom_idx = tmp;
    }
    
    int part1 = (N_ELEMENTS*(N_ELEMENTS+1) / 2 * lig_atom_idx);
    int part2 = (p1_atom_idx / 2.0f * (2 * N_ELEMENTS - p1_atom_idx + 1));
    int part3 = (p2_atom_idx-p1_atom_idx);
    return N_THETA*rs_angular_length*(part1 + part2 + part3) + th*rs_angular_length + rs;
}

double* featurize(Atom* mol_atoms, int num_mol_atom, Atom* protein_atoms, int num_protein_atom) {

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
    double** dist_lp = malloc(num_mol_atom*sizeof(double*));
    for (int i = 0; i < num_mol_atom; i++) {
        dist_lp[i] = malloc(num_protein_atom*sizeof(double));
    }

    // Setting result vector.
    int radial_length = N_ELEMENTS*N_ELEMENTS * rs_radial_length;
    int angular_length = N_ELEMENTS * binom(N_ELEMENTS+1, 2) * N_THETA * rs_angular_length;
    double* result = calloc(radial_length + angular_length, sizeof(double));

    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
			dist_lp[i][j] = eucl_dist(mol_atoms[i], protein_atoms[j]);

            // Calculate features only for pairs within radius.
            if (dist_lp[i][j] < RADIAL_CUTOFF) {
                
                // Calculate radial features.
                for (int k = 0; k<rs_radial_length; k++) {
                    int index = (rs_radial_length*N_ELEMENTS*mol_atoms[i].atom_index) + (rs_radial_length*protein_atoms[j].atom_index) + k;
                    result[index] += radial_sym_func(dist_lp[i][j], rs_radial[k]);
                }
            }
		}
    }

    // This is how to get indeces of those I need to get angles.
    // Much better. Instead of 1.017.000 it is only 42.025 which is 4%.
    int num_close_triplets = 0;
    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
            if (dist_lp[i][j] < ANGULAR_CUTOFF) {
                for (int k = j+1; k<num_protein_atom; k++) { 
                    if (dist_lp[i][k] < ANGULAR_CUTOFF) {
                        num_close_triplets++;
                    }
                }
            }
		}
    }

    int index = 0;
    Triplet* triplet_idxs = malloc(num_close_triplets * sizeof(Triplet));
    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
            if (dist_lp[i][j] < ANGULAR_CUTOFF) {
                for (int k = j+1; k<num_protein_atom; k++) { 
                    if (dist_lp[i][k] < ANGULAR_CUTOFF) {
                        triplet_idxs[index].lig_idx = i;
                        triplet_idxs[index].p1_idx = j;
                        triplet_idxs[index].p2_idx = k;
                        index++;
                    }
                }
            }
		}
    }

    double feat_sum = 0;
    for (int i = 0; i < num_close_triplets; i++) {
        int lig_idx = triplet_idxs[i].lig_idx;
        int p1_idx = triplet_idxs[i].p1_idx;
        int p2_idx = triplet_idxs[i].p2_idx;

        int lig_atom_idx = mol_atoms[lig_idx].atom_index;
        int p1_atom_idx = protein_atoms[p1_idx].atom_index;
        int p2_atom_idx = protein_atoms[p2_idx].atom_index;

        double angle = calc_angle(mol_atoms[lig_idx], protein_atoms[p1_idx], protein_atoms[p2_idx]);
        for (int l = 0 ; l<rs_angular_length ; l++) {
            for (int m = 0 ; m<N_THETA; m++) {
                int index = find_angular_index(lig_atom_idx, p1_atom_idx, p2_atom_idx, m, l, rs_angular_length);
                double feat = angular_sym_func(dist_lp[lig_idx][p1_idx], dist_lp[lig_idx][p2_idx], angle, theta_list[m], rs_angular[l]);
                result[index + radial_length] += feat;
                if (p1_atom_idx == 0 &&
                    lig_atom_idx == 0 &&
                    p2_atom_idx == 0 &&
                    m == 1 &&
                    l == 1) {
                        feat_sum += feat;
                }

                if (index == 4 && (p1_atom_idx != 0 ||
                    lig_atom_idx != 0 ||
                    p2_atom_idx != 0 ||
                    m != 1 ||
                    l != 1)) {
                        printf("FUCK: (%d, %d, %d) (th = %d, rs = %d)\n", lig_atom_idx, p1_atom_idx, p2_atom_idx, m, l);
                }
                //if (index == 200) {
                //    printf("INFO: (%d, %d, %d) (th = %d, rs = %d)\n", lig_atom_idx, p1_atom_idx, p2_atom_idx, m, l);
                //}
            }
        }
    }

    int ang_index = find_angular_index(0, 0, 8, 2, 2, rs_angular_length);
    printf("Specific: %.8f, Feat Sum: %f, Index: %d, radial_len: %d\n", result[1867], feat_sum, ang_index, radial_length);
    printf("Result has %d element\n", radial_length + angular_length);

    // Free all memory that was used.
    free(rs_radial);
    free(rs_angular);
    free(theta_list);
    for (int i = 0; i < num_mol_atom; i++) {
        free(dist_lp[i]);
    }
    free(dist_lp);
    free(triplet_idxs);

    return result;
}