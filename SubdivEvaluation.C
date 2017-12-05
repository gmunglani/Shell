#include "SubdivEvaluation.h"
#include "ShellSystem.h"

#include "face_tri3_sd.h"

// precompute the subdivision matrices P*A^n for irregular vertices
PAnArray SubdivEvaluation::PAn;

// renumbering of the shape functions of a regular patch
const unsigned int SubdivEvaluation::cvi [12] = {3,6,2,0,1,4,7,10,9,5,11,8};


// evaluates the subdivision surface of the given element at the given barycentric coordinates
void SubdivEvaluation::evaluate(ShellSystem* shell_system, unsigned int eid, Real xi, Real eta, Point& p_shell)
{
  const Tri3SD* sd_elem = static_cast<const Tri3SD*>(shell_system->get_mesh().elem(eid));
  
  p_shell.zero();

  // the 1-ring of the shell element
  const std::vector<Node*>& patch = shell_system->elem_cache[eid]->patch;
  
  // the valence of the shell element
  const unsigned int valence = sd_elem->get_ordered_valence(0);
  
  libmesh_assert(valence + 6 == patch.size());
  
  if (valence == 6) // This means that all vertices are regular, i.e. we have 12 shape functions
  {
    const Point transformed_p(xi, eta);
    for (unsigned int i = 0; i < patch.size(); ++i)
    {
      const Point& p = shell_system->node_cache[patch[i]->id()].p;
      p_shell.add_scaled(p, FE<2,SUBDIV>::shape(sd_elem, FOURTH, cvi[i], transformed_p));
    }
  }
  else // vertex 0 is irregular by construction of the mesh
  {
    const int n = rescale_irregular_coords(xi, eta);

    // find out in which subdivided patch we are and set up the "selection matrix" P and the transformation Jacobian
    unsigned int k;
    if (xi > 0.5)
    {
      k = 0;
      xi  = 2*xi - 1;
      eta = 2*eta;
    }
    else if (eta > 0.5)
    {
      k = 2;
      xi  = 2*xi;
      eta = 2*eta - 1;
    }
    else
    {
      k = 1;
      xi  = 1 - 2*xi;
      eta = 1 - 2*eta;
    }

    const Point transformed_p(xi, eta);
    
    // temporary values
    std::vector<Real> phi(12);
    for (unsigned int i = 0; i < 12; ++i)
      phi[i] = FE<2,SUBDIV>::shape(sd_elem, FOURTH, i, transformed_p);

    const DenseMatrix<Real>& this_PAn = PAn(k, valence-3, n-1);
    for (unsigned int j = 0; j < patch.size(); ++j)
    {
      Real sum = 0;
      for (unsigned int i = 0; i < 12; ++i)
        sum += this_PAn(i,j) * phi[i];
      const Point& p = shell_system->node_cache[patch[j]->id()].p;
      p_shell.add_scaled(p, sum);
    }
    
  } // end irregular vertex
}

// evaluates the subdivision surface tangents of the given element at the given barycentric coordinates
void SubdivEvaluation::evaluate_deriv(ShellSystem* shell_system, unsigned int eid, Real xi, Real eta, Point& dpdxi, Point& dpdeta)
{ 
  const Tri3SD* sd_elem = static_cast<const Tri3SD*>(shell_system->get_mesh().elem(eid));
  
  dpdxi.zero();
  dpdeta.zero();

  // the 1-ring of the shell element
  const std::vector<Node*>& patch = shell_system->elem_cache[eid]->patch;
  
  // the valence of the shell element
  const unsigned int valence = sd_elem->get_ordered_valence(0);
  
  libmesh_assert(valence + 6 == patch.size());
  
  if (valence == 6) // This means that all vertices are regular, i.e. we have 12 shape functions
  {
    const Point transformed_p(xi, eta);
    for (unsigned int i = 0; i < patch.size(); ++i)
    {
      const Point& p = shell_system->node_cache[patch[i]->id()].p;
      dpdxi.add_scaled  (p, FE<2,SUBDIV>::shape_deriv(sd_elem, FOURTH, cvi[i], 0, transformed_p));
      dpdeta.add_scaled (p, FE<2,SUBDIV>::shape_deriv(sd_elem, FOURTH, cvi[i], 1, transformed_p));
    }
  }
  else // vertex 0 is irregular by construction of the mesh
  {
    const int n = rescale_irregular_coords(xi, eta);

    // find out in which subdivided patch we are and set up the "selection matrix" P and the transformation Jacobian
    unsigned int k;
    Real jfac; // the additional factor per derivative order
    if (xi > 0.5)
    {
      k = 0;
      xi  = 2*xi - 1;
      eta = 2*eta;
      jfac = std::pow((Real)(2), n);
    }
    else if (eta > 0.5)
    {
      k = 2;
      xi  = 2*xi;
      eta = 2*eta - 1;
      jfac = std::pow((Real)(2), n);
    }
    else
    {
      k = 1;
      xi  = 1 - 2*xi;
      eta = 1 - 2*eta;
      jfac = std::pow((Real)(-2), n);
    }

    const Point transformed_p(xi, eta);
    
    // temporary values
    std::vector<Real> dphidxi(12);
    std::vector<Real> dphideta(12);

    for (unsigned int i = 0; i < 12; ++i)
    {
      dphidxi[i]  = FE<2,SUBDIV>::shape_deriv(sd_elem, FOURTH, i, 0, transformed_p);
      dphideta[i] = FE<2,SUBDIV>::shape_deriv(sd_elem, FOURTH, i, 1, transformed_p);
    }

    const DenseMatrix<Real>& this_PAn = PAn(k, valence-3, n-1);
    for (unsigned int j = 0; j < patch.size(); ++j)
    {
      Real sum_xi = 0, sum_eta = 0;
      for (unsigned int i = 0; i < 12; ++i)
      {
        sum_xi  += this_PAn(i,j) * dphidxi[i];
        sum_eta += this_PAn(i,j) * dphideta[i];
      }
      const Point& p = shell_system->node_cache[patch[j]->id()].p;
      dpdxi.add_scaled  (p, sum_xi);
      dpdeta.add_scaled (p, sum_eta);
    }
    
    dpdxi  *= jfac;
    dpdeta *= jfac;
    
  } // end irregular vertex
}

// returns the number of subdivisions required plus one, and rescales the barycentric coordinates into that subelement
int SubdivEvaluation::rescale_irregular_coords(Real& xi, Real& eta)
{
  // the limit surface can't be evaluated exactly at the irregular vertex
  if (xi + eta < irregular_tol)
  {
    xi  = std::max(xi, irregular_tol);
    eta = std::max(eta, (Real)0);
  }
  
  // transform the barycentric coordinates according to the number
  // of subdivisions required (from Stam, SIGGRAPH '99 course notes)
  const int n = (int)floor(1 - log(xi + eta) / log((Real)2));
  const Real pow2 = std::pow((Real)2, n-1);
  xi  *= pow2;
  eta *= pow2;
  
  return n;
}
