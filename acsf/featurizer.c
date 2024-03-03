#include "util.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double radial_sym_func(double Rij, double Rs, double radial_cutoff) {
    const double Rc = radial_cutoff;
    const double eta = 4.0;

    double fc = 0.5 * cos(M_PI * Rij / Rc) + 0.5;
    return exp(-eta * pow(Rij - Rs, 2)) * fc;
}

double angular_sym_func(double Rij, double Rik, double theta_ijk, double theta, double Rs, double fc_Rij, double fc_Rik, double angular_cutoff) {
    const double zeta = 8;
    const double eta = 4;

    return pow(2, 1 - zeta) * pow(1 + cos(theta_ijk - theta), zeta) * exp(-eta * pow(((Rij + Rik) / 2 - Rs), 2)) * fc_Rij * fc_Rik;
}

int angular_index(int p1_atom_idx, int lig_atom_idx, int p2_atom_idx, int th, int rs, int rs_angular_length, int num_elem, int num_theta) {
    if (p1_atom_idx > p2_atom_idx) {
        int tmp = p1_atom_idx;
        p1_atom_idx = p2_atom_idx;
        p2_atom_idx = tmp;
    }

    int index = (num_elem * (num_elem + 1) / 2 * lig_atom_idx);
    index += (p1_atom_idx / 2.0f * (2 * num_elem - p1_atom_idx + 1));
    index += (p2_atom_idx - p1_atom_idx);
    index *= num_theta * rs_angular_length;
    return index + th * rs_angular_length + rs;
}

void free_features(double *features) { free(features); }

Result featurize(Atom *mol_atoms, int num_mol_atom, Atom *protein_atoms, int num_protein_atom, Config config) {

    // Initializing radial steps.
    int rs_radial_length = ceil((config.radial_cutoff - 0.5) / config.radial_step);
    double *rs_radial = malloc(rs_radial_length * sizeof(double));

    double rad = 0.5;
    for (int i = 0; i < rs_radial_length; i++) {
        rs_radial[i] = rad;
        rad += config.radial_step;
    }

    // Initializing angular steps.
    int rs_angular_length = ceil((config.angular_cutoff - 0.5) / config.angular_step);
    double *rs_angular = malloc(rs_angular_length * sizeof(double));

    double ang = 0.5;
    for (int i = 0; i < rs_angular_length; i++) {
        rs_angular[i] = ang;
        ang += config.angular_step;
    }

    // Initializing PL distance matrix.
    double **dist = malloc(num_mol_atom * sizeof(double *));
    for (int i = 0; i < num_mol_atom; i++) {
        dist[i] = malloc(num_protein_atom * sizeof(double));
    }

    // Initializing features vector.
    int radial_length = config.num_elems * config.num_elems * rs_radial_length;
    int angular_length = config.num_elems * binom(config.num_elems + 1, 2) * config.num_theta * rs_angular_length;
    double *features = calloc(radial_length + angular_length, sizeof(double));

    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
            dist[i][j] = eucl_dist(mol_atoms[i], protein_atoms[j]);

            // Calculate features only for pairs within radius.
            if (dist[i][j] < config.radial_cutoff) {

                // Calculate radial features.
                for (int k = 0; k < rs_radial_length; k++) {
                    int index = (rs_radial_length * config.num_elems * mol_atoms[i].atom_index) + (rs_radial_length * protein_atoms[j].atom_index) + k;
                    features[index] += radial_sym_func(dist[i][j], rs_radial[k], config.radial_cutoff);
                }
            }
        }
    }

    // Initializing theta list.
    double *theta_list = malloc(config.num_theta * sizeof(double));
    double interval = 2 * M_PI / (config.num_theta);
    for (int i = 0; i < config.num_theta; i++) {
        theta_list[i] = interval * i;
    }

    // Find triplets of atoms where distance between ligand atom
    // and two protein atoms is small enough.
    for (int i = 0; i < num_mol_atom; i++) {
        for (int j = 0; j < num_protein_atom; j++) {
            if (dist[i][j] < config.angular_cutoff) {
                for (int k = j + 1; k < num_protein_atom; k++) {
                    if (dist[i][k] < config.angular_cutoff) {

                        int lig_atom_idx = mol_atoms[i].atom_index;
                        int p1_atom_idx = protein_atoms[j].atom_index;
                        int p2_atom_idx = protein_atoms[k].atom_index;

                        // Calculate angular features.
                        double angle = calc_angle(mol_atoms[i], protein_atoms[j], protein_atoms[k]);
                        double fc_Rij = 0.5 * cos(M_PI * dist[i][j] / config.angular_cutoff) + 0.5;
                        double fc_Rik = 0.5 * cos(M_PI * dist[i][k] / config.angular_cutoff) + 0.5;

                        for (int l = 0; l < rs_angular_length; l++) {
                            for (int m = 0; m < config.num_theta; m++) {
                                int index = angular_index(p1_atom_idx, lig_atom_idx, p2_atom_idx, m, l, rs_angular_length, config.num_elems, config.num_theta);
                                features[index + radial_length] +=
                                    angular_sym_func(dist[i][j], dist[i][k], angle, theta_list[m], rs_angular[l], fc_Rij, fc_Rik, config.angular_cutoff);
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < num_mol_atom; i++) {
        free(dist[i]);
    }
    free(dist);
    free(rs_radial);
    free(rs_angular);
    free(theta_list);

    Result result = {features, radial_length + angular_length};
    return result;
}