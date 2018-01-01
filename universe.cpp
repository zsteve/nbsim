/*  Universe of point particles
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#include "universe.h"

double grav_const = 6.67408e-2 /*6.67408e-11*/;

universe::universe(int num_particles){
    /* we ask for the num_particles a priori so we can reserve */
    system.reserve(num_particles);
    system2.reserve(num_particles);
    size = 0;

    current = &system;
    next = &system2;
}

universe::~universe(){
    return;
}

int universe::get_size(){
    return size;
}

particle& universe::operator[](int index){
    // assume that we are passed a valid index into the particle system
    return (*current)[index];
}

particle& universe::add_particle(particle& p){
    // add a particle to the system and returns a reference to inserted particle
    current->push_back(p);
    next->push_back(p);
    ++size;
    return current->back();
}

vector<particle>& universe::get_next(){ return *next; }
vector<particle>& universe::get_current(){ return *current; }

vector<particle>::iterator universe::begin(){
    return current->begin();
}

vector<particle>::iterator universe::end(){
    return current->end();
}

vector<particle>::iterator universe::begin_next(){
    return next->begin();
}

vector<particle>::iterator universe::end_next(){
    return next->end();
}

r3vec get_grav_atr(particle& a, particle& b){
    /*  calculate gravitational attraction between a pair of particles
        For each _pair_ of particles, we only have to compute this once.

        Returns the force on particle A, due to particle B, i.e. F(A, B)
        Force on particle B due to particle A is F(B, A) = -F(A, B) by symmetry
    */
    const double epsilon = 0.001;	// so we have a 'softened' gravitational attraction
					// correcting for the singularity at 0
    double r = particle::get_dist(a, b);
    return (grav_const*(a.m)*(b.m)/(r*r*r+epsilon))*particle::get_disp(a, b);
}

void universe::swap(){
    vector<particle>* temp = current;
    current = next;
    next = temp;
    return;
}
