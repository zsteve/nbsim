/*  Particle class
    Models a point particle of given mass and momentum
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#ifndef __PARTICLE__
#define __PARTICLE__

#include "r3vec.h"

class particle{
public:
    r3vec r;    // position
    r3vec p;    // momentum
    double m;   // mass

    r3vec f;    // force; more of an auxiliary variable 

    uint32_t c; // color, for fun!

    particle();
    particle(const r3vec& r, const r3vec& p, double m = 1);
    ~particle();

    r3vec v();  // fetch velocity

    /* get displacement of b from a */
    static r3vec get_disp(particle& a, particle& b);
    static double get_dist(particle& a, particle& b);
};

#endif