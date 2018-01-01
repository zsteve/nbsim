/*  Particle class
    Models a point particle of given mass and momentum
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#include "particle.h"
#include "r3vec.h"

particle::particle(){
    this->m = 1;    // default mass is 1
}

particle::particle(const r3vec& r, const r3vec& p, double m){
    this->r = r;
    this->p = p;
    this->m = m;
}

particle::~particle(){
    return;
}

r3vec particle::v(){
    return p*(1/m);
}

r3vec particle::get_disp(particle& a, particle& b){
    /* calculate the displacement of particle b relative to particle a */
    return b.r - a.r;
}

double particle::get_dist(particle& a, particle& b){
    /* calculate the distance of particle b relative to particle a */
    return (b.r - a.r).norm();
}

