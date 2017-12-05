#include "SmoothShellShellContact.h"
#include "ShellSystem.h" 

#include "face_tri3_sd.h"


// precompute the 64 static projection matrices (64 = 2^6, 6 being the number of constraints)
SmallMatrix<Real, SmoothShellShellContact::n_constraints, SmoothShellShellContact::n_constraints> SmoothShellShellContact::inv_AJTAJ [1 << SmoothShellShellContact::n_constraints] = { RECURSIVE_DEFINE_LIST_6(SmoothShellShellContact, inv_AJTAJ, 0) };
SmallMatrix<Real, SmoothShellShellContact::n_coords,      SmoothShellShellContact::n_coords>      SmoothShellShellContact::minus_PJ  [1 << SmoothShellShellContact::n_constraints] = { RECURSIVE_DEFINE_LIST_6(SmoothShellShellContact, minus_PJ , 0) };


SmoothShellShellContact::SmoothShellShellContact(unsigned int eid1, unsigned int eid2)
: ShellShellContact(eid1, eid2)
{
  // start searching the closests points of approach from the middle of the elements
  v[0] = w[0] = v[1] = w[1] = 1./3;
  constraint_bitmap = 0;
}

// computes the relevant data (distance, points etc.) based on the current nodal positions
void SmoothShellShellContact::compute()
{
  /*
  find min f(x) with A*x >= b, x = (v1,w1,v2,w2)^T, f(x) = |p2(v2,w2)-p1(v1,w1)|^2
      | 1  0  0  0|      | 0|
      | 0  1  0  0|      | 0|
  A = | 0  0  1  0|  b = | 0|
      | 0  0  0  1|      | 0|
      |-1 -1  0  0|      |-1|
      | 0  0 -1 -1|      |-1|
  */
  
  Real x [n_coords] = {v[0], w[0], v[1], w[1]};

  d2 = minimize_f(x);
  dc = p2 - p1;

  v[0] = x[0];
  w[0] = x[1];
  v[1] = x[2];
  w[1] = x[3];

  // this shouldn't happen if the time step is small enough
  if (d2 <= 0)
  {
#ifdef SHELL_OMP
#pragma omp critical(shell_register_problem)
#endif
    {
      system->register_problem(ShellSystem::CONTACT_PROBLEM);
      const Elem* e1 = system->get_mesh().elem(eid[0]);
      const Elem* e2 = system->get_mesh().elem(eid[1]);
      system->problem_info.set<unsigned int>("eid1") = eid[0];
      system->problem_info.set<unsigned int>("eid2") = eid[1];
      system->problem_info.set<unsigned int>("nid11") = e1->node(0);
      system->problem_info.set<unsigned int>("nid12") = e1->node(1);
      system->problem_info.set<unsigned int>("nid13") = e1->node(2);
      system->problem_info.set<unsigned int>("nid21") = e2->node(0);
      system->problem_info.set<unsigned int>("nid22") = e2->node(1);
      system->problem_info.set<unsigned int>("nid23") = e2->node(2);
      system->problem_info.set<Real>("d2") = d2;
    }
  }
}

// returns the constraint matrix A
SmallMatrix<Real, SmoothShellShellContact::n_constraints, SmoothShellShellContact::n_coords> SmoothShellShellContact::compute_A()
{
  SmallMatrix<Real, n_constraints, n_coords> A;
  
  A(0,0) = 1;
  A(1,1) = 1;
  A(2,2) = 1;
  A(3,3) = 1;
  A(4,0) = A(4,1) = -1;
  A(5,2) = A(5,3) = -1;
  
  return A;
}

// computes A * grad(f) efficiently
void SmoothShellShellContact::compute_Ag(Real* Ag, Real* g)
{
  Ag[0] = g[0];
  Ag[1] = g[1];
  Ag[2] = g[2];
  Ag[3] = g[3];
  Ag[4] = -g[0] - g[1];
  Ag[5] = -g[2] - g[3];
}

// determines the maximum allowed line search step size
Real SmoothShellShellContact::find_max_alpha(Real* x, Real* d)
{
  Real alpha = std::numeric_limits<Real>::max();

  if (d[0] < 0)
    alpha = std::min(alpha, -x[0] / d[0]);
  if (d[1] < 0)
    alpha = std::min(alpha, -x[1] / d[1]);
  if (d[2] < 0)
    alpha = std::min(alpha, -x[2] / d[2]);
  if (d[3] < 0)
    alpha = std::min(alpha, -x[3] / d[3]);
  if (d[0] + d[1] > 0)
    alpha = std::min(alpha, (1 - x[0] - x[1]) / (d[0] + d[1]));
  if (d[2] + d[3] > 0)
    alpha = std::min(alpha, (1 - x[2] - x[3]) / (d[2] + d[3]));
  
  // this probably shouldn't happen normally
  if (alpha == std::numeric_limits<Real>::max())
  {
    system->register_problem(ShellSystem::CONTACT_PROBLEM);
    const Elem* e1 = system->get_mesh().elem(eid[0]);
    const Elem* e2 = system->get_mesh().elem(eid[1]);
    system->problem_info.set<unsigned int>("eid1") = eid[0];
    system->problem_info.set<unsigned int>("eid2") = eid[1];
    system->problem_info.set<unsigned int>("nid11") = e1->node(0);
    system->problem_info.set<unsigned int>("nid12") = e1->node(1);
    system->problem_info.set<unsigned int>("nid13") = e1->node(2);
    system->problem_info.set<unsigned int>("nid21") = e2->node(0);
    system->problem_info.set<unsigned int>("nid22") = e2->node(1);
    system->problem_info.set<unsigned int>("nid23") = e2->node(2);
    system->problem_info.set<Real>("x[0]") = x[0];
    system->problem_info.set<Real>("x[1]") = x[1];
    system->problem_info.set<Real>("x[2]") = x[2];
    system->problem_info.set<Real>("x[3]") = x[3];
    
    alpha = 1;
  }

  return alpha;
}

// recalculates the bitmap of currently active constraints
void SmoothShellShellContact::update_constraint_bitmap(Real* x)
{
  // the tolerance for determining constraint activity
  static const Real constraint_tol = 100 * std::numeric_limits<Real>::epsilon();
  
  constraint_bitmap = 0;
  if (x[0]        <     constraint_tol) constraint_bitmap |= (1u << 0);
  if (x[1]        <     constraint_tol) constraint_bitmap |= (1u << 1);
  if (x[2]        <     constraint_tol) constraint_bitmap |= (1u << 2);
  if (x[3]        <     constraint_tol) constraint_bitmap |= (1u << 3);
  if (x[0] + x[1] > 1 - constraint_tol) constraint_bitmap |= (1u << 4);
  if (x[2] + x[3] > 1 - constraint_tol) constraint_bitmap |= (1u << 5);
}

// ensures that the solution doesn't even slightly violate the constraints
void SmoothShellShellContact::apply_constraints(Real* x)
{
  x[0] = std::max(x[0], (Real)0);
  x[1] = std::max(x[1], (Real)0);
  x[2] = std::max(x[2], (Real)0);
  x[3] = std::max(x[3], (Real)0);
  
  const Real u1 = 1 - x[0] - x[1];
  if (u1 < 0)
  {
    if (x[0] > x[1])
      x[0] += u1 - std::numeric_limits<Real>::epsilon();
    else
      x[1] += u1 - std::numeric_limits<Real>::epsilon();
  }
  
  const Real u2 = 1 - x[2] - x[3];
  if (u2 < 0)
  {
    if (x[2] > x[3])
      x[2] += u2 - std::numeric_limits<Real>::epsilon();
    else
      x[3] += u2 - std::numeric_limits<Real>::epsilon();
  }
}

// evaluates the square distance function using the current x = (v1,w1,v2,w2)^T
Real SmoothShellShellContact::evaluate_f(Real* x)
{
  // evaluate the shell positions at the barycentric coordinates
  evaluate(system, eid[0], x[0], x[1], p1);
  evaluate(system, eid[1], x[2], x[3], p2);

  // return the squared distance between the closest points of approach
  return (p2 - p1).size_sq();
}

// evaluates the gradient g of the square distance function
void SmoothShellShellContact::evaluate_gradf(Real* x, Real* g)
{
  // evaluate the shell tangents at the barycentric coordinates
  Point dpdxi1, dpdeta1;
  Point dpdxi2, dpdeta2;
  evaluate_deriv(system, eid[0], x[0], x[1], dpdxi1, dpdeta1);
  evaluate_deriv(system, eid[1], x[2], x[3], dpdxi2, dpdeta2);
  
  // assume that evaluate_f(x) has already computed p2 and p1
  const Point dp = p2 - p1;
  
  g[0] = -2 * (dp * dpdxi1);
  g[1] = -2 * (dp * dpdeta1);
  g[2] = 2 * (dp * dpdxi2);
  g[3] = 2 * (dp * dpdeta2);
}
