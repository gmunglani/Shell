#include "Friction.h"

// generic function to apply Margolis's slip-stick friction model (Proc. Inst. Mech. Eng., Part 1 219, 199 (2005))
void Friction::apply(Real& Efric, Point& f_contact, Real f_normal, const Point& n, const Point& v,
                     Real abs_v_normal, Real v_stop, Real spring_stiffness, Real mu_s, Real mu_d,
                     Real damping, bool& stick, Point& spring_origin, const Point& contact_point)
{
  libmesh_assert(f_normal >= 0);
  
  // compute the tangential velocity
  const Point v_normal = abs_v_normal * n;
  const Point v_tangential = v - v_normal;
  const Real abs_v_tangential = v_tangential.size();
  
  const Real v_eps = 1e-8 * v_stop;

  // static part
  if (mu_s > 0 && abs_v_tangential <= v_stop)
  { 
    // if there is no static friction spring yet, create one with initial amplitude
    // such that the force transition from dynamic to static is continuous
    if (!stick)
    {
      stick = true;

      // -k_s u - c_s k_s v_t = -mu_d f_n v_t / (|v_t| + eps), i.e.
      // u = (mu_d f_n / (k_s (|v_t| + eps)) - c_s) * v_t
      const Real prefactor = mu_d * f_normal / (spring_stiffness * (abs_v_tangential + v_eps)) - damping;
      spring_origin = contact_point - prefactor * v_tangential;
    }

    // the static friction spring is the current contact point minus the stored initial contact point
    Point u_spring = contact_point - spring_origin;

    // the spring may not be tangential to the contact, so we need to find the tangential
    // component of the spring direction (i.e. subtract the normal component)
    u_spring -= (u_spring * n) * n;

    // the potential energy contained in the static friction spring
    Efric = 0.5 * spring_stiffness * u_spring.size_sq();

    // Hooke's law (eq. (18))
    f_contact = -spring_stiffness * u_spring;
    
    // add the tangential damping force
    f_contact -= damping * spring_stiffness * v_tangential;
  }

  // dynamic part
  if (abs_v_tangential > v_stop || f_contact.size_sq() >= std::pow(mu_s * f_normal, 2))
  {
    // remove the static friction spring
    spring_origin.zero();
    stick = false;
    Efric = 0;

    // the elements are slipping, so apply dynamic friction (eq. (16))
    f_contact = -mu_d * f_normal / (abs_v_tangential + v_eps) * v_tangential;
  }
}
