#ifndef SHELLSHELLCONTACT_H_
#define SHELLSHELLCONTACT_H_

#include "Contact.h"
#include "DataVector.h"


class ShellSystem; // forward declaration


class ShellShellContact : public Contact
{
  friend class ShellSystem;
  
  friend std::ostream& operator<<(std::ostream& os, const ShellShellContact& c);

protected:
  Point* p[2][3];
  Real v[2], w[2];

public:
  static ShellSystem* system;

  ShellShellContact(unsigned int eid1, unsigned int eid2);

  virtual void compute();
  void compute_on_mesh();
  virtual Point point();

protected:
  void edge_edge_test();
  void vertex_face_test();
};

#endif /* SHELLSHELLCONTACT_H_ */
