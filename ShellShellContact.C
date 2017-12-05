#include "ShellShellContact.h"
#include "ShellSystem.h"

ShellSystem* ShellShellContact::system = 0;

ShellShellContact::ShellShellContact(unsigned int eid1, unsigned int eid2)
: Contact(eid1, eid2)
{
  libmesh_assert(system);
  libmesh_assert(eid1 < system->get_mesh().n_elem());
  libmesh_assert(eid2 < system->get_mesh().n_elem());
  
  // store pointers to the nodal positions
  for (unsigned int e = 0; e < 2; ++e)
  {
    const Elem* elem = system->get_mesh().elem(eid[e]);
    for (unsigned int n = 0; n < 3; ++n)
      p[e][n] = &(system->node_cache[elem->node(n)].p);
  }
}

// prints a brief version of the contact data
std::ostream& operator<<(std::ostream& os, const ShellShellContact& c)
{
  const Elem* e1 = c.system->get_mesh().elem(c.eid[0]);
  const Elem* e2 = c.system->get_mesh().elem(c.eid[1]);
  os << "eid = [" << c.eid[0] << "," << c.eid[1] << "], "
     << "nid1 = [" << e1->node(0) << "," << e1->node(1) << "," << e1->node(2) << "], "
     << "nid2 = [" << e2->node(0) << "," << e2->node(1) << "," << e2->node(2) << "], "
     << "d = " << std::sqrt(c.d2) << ", v1 = " << c.v[0] << ", w1 = " << c.w[0] << ", v2 = " << c.v[1] << ", w2 = " << c.w[1] << std::endl;
  return os;
}

// computes the relevant data (distance, points etc.)
void ShellShellContact::compute()
{
  d2 = std::numeric_limits<Real>::max();
  
  edge_edge_test();
  vertex_face_test();

  // this shouldn't happen if the time step is small enough
  if (d2 <= 0)
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

// computes the distance using a the undeformed mesh
void ShellShellContact::compute_on_mesh()
{
  // use the nodal mesh positions instead of p
  for (unsigned int e = 0; e < 2; ++e)
  {
    const Elem* elem = system->get_mesh().elem(eid[e]);
    for (unsigned int n = 0; n < 3; ++n)
      p[e][n] = elem->get_node(n);
  }

  d2 = std::numeric_limits<Real>::max();

  edge_edge_test();
  vertex_face_test();

  // don't assert that d2 > 0 here, because compute_on_mesh() is called to find
  // the topological neighborhood of elements, in which case we might have d2 == 0
}

// returns the point of contact
Point ShellShellContact::point()
{
  return (1 - v[0] - w[0]) * *p[0][0] + v[0] * *p[0][1] + w[0] * *p[0][2] + 0.5 * dc;
}

// tests if any of the two edges of the two elements are close
void ShellShellContact::edge_edge_test()
{
  static const unsigned int mod3 [3] = {1, 2, 0};
  Real st[2];
  Point this_dc;

  // test each edge of element 1 with each edge of element 2
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      segment_segment_distance(*p[0][i], *p[0][mod3[i]], *p[1][j], *p[1][mod3[j]], this_dc, st);

      const Real this_d2 = this_dc.size_sq();

      // store the found contact data if the edges are touching and closer than
      // all previously found contacting edge pairs
      if (this_d2 < d2)
      {
        d2 = this_d2;
        dc = this_dc;
        v[0] = (i==0) * st[0] + (i==1) * (1-st[0]);
        w[0] = (i==1) * st[0] + (i==2) * (1-st[0]);
        v[1] = (j==0) * st[1] + (j==1) * (1-st[1]);
        w[1] = (j==1) * st[1] + (j==2) * (1-st[1]);
      }
    }
  }
}

// tests if any vertex of one element is close to the face of the other
void ShellShellContact::vertex_face_test()
{
  Real vw[2];
  Point this_dc;

  // test each vertex of each element with the face of the other
  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      point_triangle_distance(*p[i][j], *p[1-i][0], *p[1-i][1], *p[1-i][2], this_dc, vw);

      const Real this_d2 = this_dc.size_sq();

      // store the found contact data if the vertex is touching the face
      // and is closer than all previously found contacts
      if (this_d2 < d2)
      {
        d2 = this_d2;
        dc = (i == 0 ? this_dc : -this_dc); // dc must always points from element 1 to element 2
        v[0] = (i==0) * (j==1) + (i==1) * vw[0];
        w[0] = (i==0) * (j==2) + (i==1) * vw[1];
        v[1] = (i==1) * (j==1) + (i==0) * vw[0];
        w[1] = (i==1) * (j==2) + (i==0) * vw[1];
      }
    }
  }
}
