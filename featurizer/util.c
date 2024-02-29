#include <math.h>
#include "util.h"

double eucl_dist(Atom p1, Atom p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;

    return sqrt(dx*dx + dy*dy + dz*dz);
}

double dot_product(Atom a, Atom b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Atom cross_product(Atom a, Atom b) {
    Atom result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

double magnitude(Atom a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double calc_angle(Atom a, Atom b, Atom c) {
    Atom ab, ac;
    double dot_ab_ac, mag_ab, mag_ac, angle_rad;

    ab.x = b.x - a.x;
    ab.y = b.y - a.y;
    ab.z = b.z - a.z;
    ac.x = c.x - a.x;
    ac.y = c.y - a.y;
    ac.z = c.z - a.z;

    dot_ab_ac = dot_product(ab, ac);

    mag_ab = magnitude(ab);
    mag_ac = magnitude(ac);

    angle_rad = acos(dot_ab_ac / (mag_ab * mag_ac));

    return angle_rad;
}
