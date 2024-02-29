#ifndef UTIL_H
#define UTIL_H

typedef struct {
    double x;
    double y;
    double z;
    int atom_index;
} Atom;

typedef struct {
    int lig_idx;
    int p1_idx;
    int p2_idx;
} Triplet;

double eucl_dist(Atom p1, Atom p2);
double dot_product(Atom a, Atom b);
double magnitude(Atom a);
double calc_angle(Atom a, Atom b, Atom c);
#endif
