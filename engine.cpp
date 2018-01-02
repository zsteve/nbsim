/*  Particle system simulation engine
    Models a system of point particles
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#include <iostream>
#include <cmath>
#include <tuple>
#include <random>

#include "particle.h"
#include "universe.h"

#include "minifb/include/MiniFB.h"

using namespace std;

std::random_device rd;
std::mt19937 gen(rd());

const int screen_width = 640;
const int screen_height = 480;
const int pi = 3.14159265359;

extern double grav_const;

const double dt = 0.005; // time step

typedef unsigned int uint32_t;

uint32_t buffer[screen_width*screen_height];

// scale our buffer by factor 2
uint32_t buffer_big[screen_width*screen_height*4];

#define MAKE_COLOR(r, g, b) (r << 16) | (g <<  8) | b

const double camera_x = 0;
const double camera_y = 0;
const double camera_z = 10;

// project point from 3d space to 2d space
std::tuple<double, double> project(double x, double y, double z){
    return std::make_tuple(-(x-camera_x)/(z-camera_z),\
                 -(y-camera_y)/(z-camera_z));
}

std::tuple<double, double> project(r3vec& v){
    return project(v.x, v.y, v.z);
}

void init_universe(universe& u, int n){
    static const r3vec k(1, 1, 1);
    static const int center_mass = 1;

    // initialise our universe
    // for the time being, randomly generate points on unit sphere
    uniform_real_distribution<double> dist(0, 1);   // uniform distribution
    particle p(r3vec(0, 0, 0), r3vec(), center_mass);
    p.c = MAKE_COLOR(255, 0, 0);
    u.add_particle(p);

    // for fun, generate some colours
    uniform_int_distribution<int> color_dist(0, 255);

    
    for(int i = 1 ; i < n; i++){
        // want to put other particles in orbit about the central particle
        double orbit_r = 1;
        double orbit_v = sqrt(grav_const*center_mass/orbit_r)*1.0;

        // generate polar and azimuthal angles randomly
        double theta = dist(gen)*pi*2;
        double phi = dist(gen)*pi;
        double x = orbit_r*cos(theta)*sin(phi);
        double y = orbit_r*sin(theta)*sin(phi);
        double z = orbit_r*cos(phi);

        // generate uniform momentum
        r3vec m(x, y, z);
        m = m.crossp(k).unit()*orbit_v*0;

        particle p(r3vec(x, y, z), m);
        p.c = MAKE_COLOR(color_dist(gen), color_dist(gen), color_dist(gen));
        u.add_particle(p);
    }
}

#define INTEGRATION_VERLET

void simulate_step(universe& u){

#ifdef INTEGRATION_VERLET
static int first_exec = 0;
#endif

#ifdef INTEGRATION_VERLET
    if(!first_exec){
#endif
    // set all net forces to zero first 
    for(int i = 0; i < u.get_size(); i++) u[i].f = r3vec(0, 0, 0);
    // now compute net force on each particle
    for(int i = 0; i < u.get_size(); i++){
        r3vec f_i();    // net force on particle i
        for(int j = i+1; j < u.get_size(); j++){
            // loop through all particles not previously considered, j
            r3vec f_ij = get_grav_atr(u[i], u[j]);
            u[i].f+=f_ij;
            u[j].f-=f_ij;
        }
    }
#ifdef INTEGRATION_VERLET
    first_exec = 1;
    }
#endif

    // computed resultant forces, now numerically integrate to
    // solve for momentum and position
    for(int i = 0; i < u.get_size(); i++){
#ifdef INTEGRATION_EULER
        // p(t + dt) = p(t) + F_tot(t) * dt
        u.get_next()[i].p = u[i].p + u[i].f*dt;
        // r(t + dt) = r(t) + 1/m*p(t) * dt
        u.get_next()[i].r = u[i].r + 1/(u[i].m)*u[i].p*dt;
#elif defined INTEGRATION_VERLET
        vector<particle>& v = u.get_next();

        // r(t + dt) = r(t) + p(t)dt/m + 0.5F(t)*dt^2/m
        v[i].r = u[i].r + u[i].p*dt*(1/u[i].m) + 0.5*u[i].f*dt*dt*(1/u[i].m);

        // need to compute forces for t+dt
        // set all net forces to zero first 
        for(int i = 0; i < v.size(); i++)
            v[i].f = r3vec(0, 0, 0);
        // now compute net force on each particle
        for(int i = 0; i < v.size(); i++){
            r3vec f_i();    // net force on particle i
            for(int j = i+1; j < v.size(); j++){
                // loop through all particles not previously considered, j
                r3vec f_ij = get_grav_atr(v[i], v[j]);
                v[i].f+=f_ij;
                v[j].f-=f_ij;
            }
        }
        // p(t + dt) = p(t) + 0.5(F(t) + F(t + dt))dt
        u.get_next()[i].p = u[i].p + 0.5*(v[i].f + u[i].f)*dt;
#endif
    }

    u.swap();

    return;
}

void plot_pixel(int x, int y, uint32_t c){
    // plot pixel using raw system coordinates
    if(x >= 0 && x < screen_width && y >= 0 && y < screen_height) 
        buffer[y*screen_width+x] = c;
}

void plot_pixel2(int x, int y, uint32_t c){
    // plot pixel using mathematical coordinate system centered (0, 0)
    plot_pixel(x+screen_width/2, y+screen_width/2, c);
}

void scale_up(){
    for(int i = 0; i < screen_height<<1; i++){
        for(int j = 0; j < screen_width<<1; j++){
            buffer_big[i*(screen_width<<1)+j] = buffer[(i >> 1)*screen_width+(j >> 1)];
        }
    }
}

const int num_particles = 5;
const double zoom_factor_x = 1000;
const double zoom_factor_y = 1000;

int main(){
    universe u(num_particles);
    init_universe(u, num_particles);

    if(!mfb_open("nbsim", screen_width<<1, screen_height<<1))
        return -1;
    
    while(1){
        int state;

        simulate_step(u);

        for(int i = 0; i < num_particles; i++){
            // only plot pixel if its visible
            if(u[i].r.z > camera_z) continue;
            std::tuple<double, double> pos = project(u[i].r);
            plot_pixel2(std::get<0>(pos)*zoom_factor_x, \
                        std::get<1>(pos)*zoom_factor_y, \
                        u[i].c);
        }
        scale_up();
        state = mfb_update(buffer_big);
        if(state < 0) break;
    }

    mfb_close();

    return 0;
}
