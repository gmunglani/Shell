#ifndef FRICTION_H_
#define FRICTION_H_

#include "point.h"

class Friction
{
private:
  Friction() {} // this class cannot be instatiated

public:
  static void apply(Real& Efric, Point& f_contact, Real f_normal, const Point& n, const Point& v,
                    Real abs_v_normal, Real v_stop, Real spring_stiffness, Real mu_s, Real mu_d,
                    Real damping, bool& stick, Point& spring_origin, const Point& contact_point);
};

#endif /* FRICTION_H_ */
