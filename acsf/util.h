#ifndef UTIL_H
#define UTIL_H

typedef struct {
    double x;
    double y;
    double z;
    int atom_index;
} Atom;

typedef struct {
    double radial_cutoff;
    double angular_cutoff;
    double radial_step;
    double angular_step;
    int num_theta;
    int num_elems;
    int num_mols;
} Config;

typedef struct {
    double* features;
    int num_rows;
    int num_cols;
} Result;

double eucl_dist(Atom p1, Atom p2);
double dot_product(Atom a, Atom b);
double magnitude(Atom a);
double calc_angle(Atom a, Atom b, Atom c);
int fac(int n);
int binom(int n, int k);
#endif
