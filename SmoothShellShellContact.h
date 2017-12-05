#ifndef SMOOTHSHELLSHELLCONTACT_H_
#define SMOOTHSHELLSHELLCONTACT_H_

#include "ShellShellContact.h"
#include "SubdivEvaluation.h"
#include "ConstrainedMinimizer.h"

class ShellSystem;


class SmoothShellShellContact : public ShellShellContact, public SubdivEvaluation, public ConstrainedMinimizer<SmoothShellShellContact, 6, 4>
{
  friend class ShellSystem;
  
  // needed for the CRTP
  friend class ConstrainedMinimizer<SmoothShellShellContact, n_constraints, n_coords>;

public:
  SmoothShellShellContact(unsigned int eid1, unsigned int eid2);
  
  virtual void compute();
  
protected:
  Point p1, p2;

  static SmallMatrix<Real, n_constraints, n_constraints> inv_AJTAJ [1 << n_constraints];
  static SmallMatrix<Real, n_coords, n_coords> minus_PJ [1 << n_constraints];
  
  static SmallMatrix<Real, n_constraints, n_coords> compute_A();
  
  virtual Real evaluate_f(Real* x);
  virtual void evaluate_gradf(Real* x, Real* g);
  
  virtual void compute_Ag(Real* Ag, Real* g);
  virtual Real find_max_alpha(Real* x, Real* d);
  virtual void update_constraint_bitmap(Real* x);
  virtual void apply_constraints(Real* x);
};

#endif /* SMOOTHSHELLSHELLCONTACT_H_ */
