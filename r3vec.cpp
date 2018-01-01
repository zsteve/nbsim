/*  R3 vector class
    Implementation of vector/linear algebra operations
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#include "r3vec.h"
#include <cmath>
#include <iostream>

using namespace std;

r3vec::r3vec(){
    x = 0; y = 0; z = 0;
}

r3vec::r3vec(double x, double y, double z){
    this->x = x;
    this->y = y;
    this->z = z;
}

r3vec::~r3vec(){
    return;
}

r3vec operator*(const r3vec& lhs, double rhs){
    // scalar multiplication
    return r3vec(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
} 

r3vec operator*(double lhs, const r3vec& rhs){
    // scalar multiplication
    return r3vec(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
}

r3vec r3vec::operator+(const r3vec& rhs){
    // addition
    return r3vec(this->x+rhs.x, this->y+rhs.y, this->z+rhs.z);
}

r3vec r3vec::operator+=(const r3vec& rhs){
    this->x+=rhs.x; this->y+=rhs.y; this->z+=rhs.z;
    return *this;
}

r3vec r3vec::operator-(const r3vec& rhs){
    // subtraction
    return r3vec(this->x-rhs.x, this->y-rhs.y, this->z-rhs.z);
}

r3vec r3vec::operator-=(const r3vec& rhs){
    this->x-=rhs.x; this->y-=rhs.y; this->z-=rhs.z;
    return *this;
}

double r3vec::norm(){
    // return Euclidean norm
    return sqrt(x*x + y*y + z*z);
}

double r3vec::dotp(const r3vec& rhs){
    // return Euclidean inner product
    // TODO: maybe implement generalised inner products in future?
    return (x*rhs.x + y*rhs.y + z*rhs.z);
}

r3vec r3vec::crossp(const r3vec& rhs){
    // returns cross product (return) = *this x rhs
    return r3vec(this->y*rhs.z - this->z*rhs.y, 
                this->z*rhs.x - this->x*rhs.z, 
                this->x*rhs.y - this->y*rhs.x);
}

r3vec r3vec::unit(){
    // returns unit vector
    return (*this)*(1/this->norm());
}

ostream& operator<<(ostream& os, const r3vec& v){
    os << '(' << v.x << ' ' << v.y << ' ' << v.z << ')';
}
