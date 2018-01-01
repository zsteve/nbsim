/*  R3 vector class
    Implementation of vector/linear algebra operations, optimised for speed
    Author: Stephen Zhang; stephenz@student.unimelb.edu.au
*/

#ifndef __R3VEC__
#define __R3VEC__

#include <cmath>
#include <ostream>

using namespace std;

class r3vec{
public:
    double x, y, z;

    r3vec();
    r3vec(double x, double y, double z);
    ~r3vec();

    friend r3vec operator*(const r3vec& lhs, double rhs);
    friend r3vec operator*(double lhs, const r3vec& rhs);

    r3vec operator+(const r3vec& rhs);
    r3vec operator+=(const r3vec& rhs);
    r3vec operator-(const r3vec& rhs);
    r3vec operator-=(const r3vec& rhs);

    friend ostream& operator<<(ostream& os, const r3vec& v);

    double norm();
    double dotp(const r3vec& rhs);
    r3vec crossp(const r3vec& rhs);
    r3vec unit();
};

#endif