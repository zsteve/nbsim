/*  Universe of point particles
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#ifndef __UNIVERSE__
#define __UNIVERSE__

#include <cmath>
#include <vector>
#include "particle.h"

using namespace std;

class universe{
private:
    vector<particle> system;
    vector<particle> system2;

    int size;
    
    vector<particle> *current, *next;
public:
    /* universe class - contains all particles */
    universe(int num_particles);
    ~universe();

    int get_size();

    particle& operator[](int index);
    particle& add_particle(particle& p);

    vector<particle>& get_next();
    vector<particle>& get_current();

    vector<particle>::iterator begin();
    vector<particle>::iterator end();

    vector<particle>::iterator begin_next();
    vector<particle>::iterator end_next();

    void swap();
};

r3vec get_grav_atr(particle& a, particle& b);

#endif