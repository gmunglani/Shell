#ifndef __shellsystem_h__
#define __shellsystem_h__

#include "equation_systems.h"
#include "system.h"
#include "fe_base.h"
#include "quadrature.h"
#include "dense_vector.h"
#include "dense_subvector.h"
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "getpot.h"

#include <list>
#include <map>

#include "macros.h"
#include "Function.h"
#include "DataVector.h"
#include "BoundingBox.h"
#include "MovingAverage.h"


class ShellShellContact;
class SubdivEvaluation;


class ShellSystem : public System
{
  friend class ShellShellContact;
  friend class SubdivEvaluation;
  friend class SmoothShellShellContact;

public:
  // these are for better performance (allows optimizations at compile time)
  enum {UX_VAR, UY_VAR, UZ_VAR};

  enum BoundaryCondition {CLAMPED, CLAMPED_WEAK, PINNED, PINNED_WEAK,
                          FREE, FREE_WEAK, INVALID_BC};

  enum CoordinateSystem {LOCAL, CARTESIAN, CYLINDRICAL, SPHERICAL, INVALID_CS};

  enum Cavity {BOX, CYLINDER, ELLIPSOID, INVALID_CAV};
  
  enum ContactType {LINEAR, SMOOTH, NONE, INVALID_CT};

  // problems that occur during execution, in order of increasing badness
  enum Problem {NO_PROBLEM, NO_DISP_CHANGE, LARGE_DISP,
                STOP_SIM, TOO_MANY_CELLS, SMALL_TIMESTEP, CONTACT_PROBLEM, NOT_FINITE};

#ifdef MEASURE_PERFORMANCE
  PerfLog perf_log;
#endif

protected:
  struct SurfaceMetric
  {
    SurfaceMetric() : jac(1) {}
    RealVectorValue xyz; // the location of the quadrature points
    RealVectorValue vec[3], dvec[3]; // a_1 a_2 a_3 and a_{1,1} a_{1,2} a_{2,2}
    RealVectorValue ff1, ff2; // first and second fundamental form in Voigt notation (without factor 2)
    Real jac, sz; // surface Jacobian (Jacobian matrix determinant)
  };

  struct NodeCache
  {
    NodeCache() : is_ghost(true) {}
    bool is_ghost; // whether or not the node is ghosted
    Point p; // the nodal position (mesh coordinate plus displacement)
    Point spring_origin; // the origin of the static cavity friction spring
    std::vector<const Node*> neighbors;
  };

  struct ElemCache
  {
    ElemCache() : area(0),
                  JxW(NULL), phi(NULL), dphi(NULL), d2phi(NULL),
                  a_ref(NULL), a_def(NULL), H(NULL), G(NULL), G_old(NULL),
                  GV(NULL), alpha(NULL), beta(NULL), qxyz(NULL), odxi(NULL), odeta(NULL) {}
    ~ElemCache() { delete [] JxW;
                   delete [] a_ref; delete [] a_def; delete [] H; delete [] G; delete [] G_old;
                   delete [] GV; delete [] alpha; delete [] beta; delete [] qxyz; delete [] odxi; delete [] odeta; }
    Real area; // area of the deformed limit surface element
    std::vector<unsigned int> dof_indices, dof_indices_ux, dof_indices_uy, dof_indices_uz; // the global dof indices
    Real* JxW; // jacobian times quadrature weight
    Real* phi; // the shape functions evaluated at the quadrature points
    RealGradient* dphi; // first derivatives
    RealTensor* d2phi; // second derivatives
    SurfaceMetric* a_ref, * a_def; // surface metrics on reference and deformed configuration
    RealTensorValue* H; // matrix H from Koiter shell theory in Voigt notation
    RealVectorValue* G, * G_old; // growth tensors in Voigt notation
    RealTensorValue* GV; // inverse growth tensor in Voigt notation
    RealVectorValue* alpha, * beta; // membrane and bending strain tensors in Voigt notation
    Point* qxyz;
	RealGradient* odxi, * odeta;
    BoundingBox<Real> bounds;
    std::set<const Elem*> neighbor_elems; // the direct topological neighbor elements
    std::set<const Elem*> close_elems; // elements topologically close enough to be excluded from contact
    std::list<unsigned int> in_contact_with; // element ids this element is in contact with
    std::vector<Node*> patch; // the 1-ring of the element
  };
  
  struct BoundaryCache
  {
    BoundaryCache(unsigned int _nids[]) { nids[0] = _nids[0]; nids[1] = _nids[1]; nids[2] = _nids[2]; nids[3] = _nids[3]; }
    unsigned int nids [4]; // the four node ids involved in each boundary edge
  };

  // a pointer to a stream that should be used for output
  std::ostream* os;
  
  int current_problem, last_problem;
  Parameters problem_info;

  // data vectors
  ScalarVector u, v, a, f, Minv, C, u_pred, v_pred, a_pred;
  
  // measurement vectors
  ScalarVector kGauss_vec, kMean_vec, ekin_vec, etens_vec, ebend_vec, ecav_vec, eself_vec, efric_vec, mass_vec;
  std::map<std::string, ScalarVector*> node_vectors, elem_vectors;
  
  // measurement scalars
  unsigned int Nc, Nc_cav; // number of contacts
  Real Eself, Ecav, Efric, Etens, Ebend; // static energies
  Real cavity_force; // forces

  // material properties
  Real shell_mass_density; // mass density
  Real shell_young_modulus; // Young's modulus
  Real shell_young_modulus2; // Young's modulus 2
  Real shell_poisson_ratio; // poisson ratio
  Real shell_shear_modulus; // shell shear modulus
  Real memb_stiff_mult; // membrane stiffness multiplier
  Real bend_stiff_mult; // bending stiffness multiplier
  Real bulk_modulus; // bulk modulus for enclosed volume

  // geometric properties
  Real thickness0, thickness; // initial and current shell thickness
  Cavity cavity;
  Real Rx, Ry, Rz; // semiaxes of the cavity
  Real cavity_surface; // surface area of the cavity
  Real cavity_volume; // volume of the cavity

  // boundary condition properties
  BoundaryCondition shell_boundary;
  unsigned int bc_coords; // bitmask indicating in which of the 3 dimensions the BCs are applied
  unsigned int bc_every; // apply boundary conditions every x timesteps
  Real bc_stiff_mult; // boundary condition stiffness multiplier
  Real bc_stiff; // the stiffness of weak BC springs
  Real bc_criterion;

  // dynamic properties
  Real damping_alpha; // damping mass prefactor
  Real damping_beta; // damping stiffness prefactor
  Real damping_gamma; // viscous damping constant
  Real damping_rotation; // prefactor to a force opposed to the angular velocity
  Real perturbation_rate; // the rate at which the shell is perturbed
  Real perturbation_energy; // the amount of thermal energy added at the perturbation rate
  Real gravity; // the gravity constant
  unsigned int gravity_dir; // direction in which gravity is applied
  Real pressure; // surface pressure in a certain direction
  CoordinateSystem pressure_base; // basis in which the pressure is applied
  unsigned int pressure_dir; // direction in pressure_base in which the pressure is applied
  Real growth1, growth2; // growth factors in the two in-plane directions
  bool grow_dynamically; // whether or not non-prescribed growth is applied
  bool grow_in_plane; // whether or not in-plane growth is applied
  CoordinateSystem growth_base; // basis in which the growth tensor is diagonal
  unsigned int growth_dir [3]; // coordinate order in growth_base
  Real metric_growth_rate; // proportionality constant for growth in direction of the metric tensor
  Real curvature_growth_rate; // proportionality constant for growth in direction of the curvature tensor
  bool grow_mass; // whether mass should be recalculated in volumetric growth
  Real area; // the total surface area of the shell
  Real volume, volume0; // current and initial volume enclosed by the shell
  Real volume_growth; // growth factor for the enclosed volume
  bool stand_still; // whether or not the shell should stand still

  // contact-related parameters
  ContactType shell_shell_contact; // the type of self-contact
  bool shell_cavity_contact; // enable/disable shell-cavity interaction
  Real stiff_density; // spring stiffness density (k/area) for contacts
  Real mu_s; // static friction coefficient for self-contacts
  Real mu_d; // dynamic friction coefficient for self-contacts
  Real mu_s_cavity; // static friction coefficient for cavity contacts
  Real mu_d_cavity; // dynamic friction coefficient for cavity contacts
  bool self_friction, cavity_friction;
  Real self_damping; // damping of self-interaction normal forces
  Real cavity_damping; // damping of cavity interaction normal forces
  Real self_friction_damping; // viscous damping coefficient for the auxiliary static friction spring in self-contact
  Real cavity_friction_damping; // viscous damping coefficient for the auxiliary static friction spring in cavity contact
  Real friction_stop_velocity; // below which velocity dynamic friction should become static
  Real fric_stiff_mult; // static friction stiffness multiplier
  Real stiff_static_friction; // stiffness coefficient of the auxiliary spring for static friction
  Real static_friction_tol; // how far apart two contacting objects can be for keeping the spring alive
  Real friction_v_eps; // epsilon in Margolis' eq. (9) for the normalization of small tangential velocities
  Real d2_touch, d2_friction; // thickness^2, i.e. the squared distance at which two elements touch (plus friction margin)
  Real one_plus_delta_friction; // 1 + Delta s.t. depth = - static_friction_tol * thickness
  Real neighbor_margin; // multiple of the shell thickness to add to the topological neighborhood size
  Real Rx_eff, Ry_eff, Rz_eff; // effective cavity semiaxes
  Real avg_R_eff; // average effective cavity semiaxis
  Real inv_R2_eff [3]; // inverse squared effective cavity semiaxes
  std::vector<ShellShellContact*> shell_shell_contact_candidates;

  // constants for the plate theory
  Real K; // membrane stiffness or "stretching stiffness"
  Real D; // bending stiffness or "flexural rigidity"

  // cached data
  unsigned int n_qp; // number of quadrature points per element
  std::vector<Real> qw; // quadrature weights
  std::vector<std::vector<const Elem*> > nodes_to_elem_map; // a map from nodes to incident elements
  std::vector<NodeCache> node_cache;
  std::vector<ElemCache*> elem_cache;
  std::vector<BoundaryCache> boundary_cache;
  std::map<unsigned int, Real*> phi_all;
  std::map<unsigned int, RealGradient*> dphi_all;
  std::map<unsigned int, RealTensor*> d2phi_all;

  // time integration
  Real dt; // the current time increment
  Real newmark_beta, newmark_gamma;
  Real newmark_a [10];
  Real stopping_criterion;
  bool stop_on_problem;

  // adaptive stepsize control
  unsigned int adapt_dt_every; // try to adapt the timestep every x steps
  Real dt_min, dt_max; // minimum and maximum time increment
  Real rel_loc_err_target; // the target relative local error per timestep
  Real rel_loc_err_max; // upper bound
  Real rel_loc_err_min; // lower bound

  // numeric stuff
  Real distortion; // the fraction of the shell thickness for initial distortion

  // linked cell stuff
  Real du_sum; // the cumulative maximum displacement (for checking when to rebuild the linked cells)
  Real cell_size_mult; // multiple of the average element size for the size of the linked cells
  Real cell_margin; // prefactor to the additional margin given to AABBs to reuse them for a while
  Real avg_elem_size; // used for estimating the best cell size
  BoundingBox<Real> bounds; // the minimal bounding box containing the shell
  unsigned int n_cells_max; // maximum number of cells expected

  // function parsers
  std::map<std::string, Function*> functions;
  Point xyz; // temporary object for evaluating functions at coordinate (x,y,z)

  // other stuff
  Real increment_factor;
  unsigned int time_step;
  int extra_order; // > 0 if more than a single quadrature point per element is desired
  Real density; // packing density
  unsigned int n_ghost_nodes; // number of ghost nodes
  unsigned int n_ghost_elem; // number of ghost element
  timespec start_time; // the wall clock time when the program started
  Real wall_time; // no. of seconds since the program started
  bool const_derived, const_mass_damping; // whether or not some variables are constant during the simulation
  static const RealTensorValue Id; // 3x3 identity matrix

  // moving averages
  MovingAverage<Real> avg_Ekin, avg_Etens, avg_Ebend, avg_Ecav, avg_Eself, avg_Efric,
                      avg_area, avg_vol, avg_pcav, avg_kGauss, avg_kMean, avg_Nc, avg_Nc_cav, avg_dt;

public:
  ShellSystem(EquationSystems& es, const std::string& name, const unsigned int number);
  virtual ~ShellSystem();

  void set_ostream(std::ostream& o) { os = &o; }
  void set_parameters(GetPot& config_file, GetPot& command_line);
  void update_derived();
  void build_cache();
  void build_topological_neighborhood();
  void set_newmark_parameters(Real new_dt);
  void apply_initial_conditions(std::ifstream & loadfile);
  void compute_mass_damping();

  bool check_stopping_criterion();
  void prepare_next_step();
  bool advance(unsigned int next_time_step);
  bool check_for_problems();
  void register_problem(Problem p) { if (p > current_problem) current_problem = p; }

  void reset_properties();
  void measure_properties();
  inline void print_header();
  void print_properties();
  void save(std::ofstream & savefile);
  void write_vtk(std::string filename);
  void write_basis_vtk1(std::string filename);
  void write_basis_vtk2(std::string filename);

  const timespec& get_start_time() { return start_time;    }
  unsigned int get_n_ghost_nodes() { return n_ghost_nodes; }
  unsigned int get_n_ghost_elem()  { return n_ghost_elem;  }

protected:
  void compute_position();
  void compute_enclosed_volume();
  void compute_bounding_box(BoundingBox<Real>& bb, Real& avg_elem_size);
  void apply_boundary_conditions();
  void compute_reference_metrics(const unsigned int eid, const AutoPtr<FEBase>& fe);
  void compute_deformed_metrics();
  void compute_external_force();
  void subtract_internal_force();
  void compute_shell_shell_contact_candidates(bool force = false);
  void add_shell_shell_contact_force();
  void perturb();
  void reinit_mod(const Elem* elem, const AutoPtr<FEBase>& fe, std::vector<Real> qw);
  void compute_control_map(const unsigned int dim, const std::vector<Real>& qw, const Elem* elem);  

  Real avg_area_node(unsigned int nid);
  void check_increment(bool require_change = true);
  void compute_packing_density();
  void apply_growth();
  void reject_growth();

  static inline void cross_product_matrix(RealTensorValue& dest, const RealVectorValue& z);
  static inline RealTensorValue inverse(const RealTensorValue& A);

  inline Real cavity_isoparameter(const unsigned int nid);
  inline void get_nodal_props(const ScalarVector& vector, unsigned int eid, Point& p1, Point& p2, Point& p3);
  inline bool elems_are_close(unsigned int eid1, unsigned int eid2);

  // contact stuff
  inline Real shell_shell_contact_force(Real stiffness, Real distance);
  inline Real shell_shell_contact_energy(Real stiffness, Real distance);
  static inline Real contact_force(Real stiffness, Real depth);
  static inline Real contact_energy(Real stiffness, Real depth);

  void load(std::ifstream & loadfile);
};


// inlined member functions:

// calculates the cavity isoparameter 1+delta
Real ShellSystem::cavity_isoparameter(const unsigned int nid)
{
  const Point& p = node_cache[nid].p;
  Real x_red = p(0) / Rx_eff;
  Real y_red = p(1) / Ry_eff;
  Real z_red = p(2) / Rz_eff;
  x_red = x_red * x_red;
  y_red = y_red * y_red;
  z_red = z_red * z_red;

  switch (cavity)
  {
    case BOX:
    {
      Real max_red = x_red;
      if (y_red > max_red) max_red = y_red;
      if (z_red > max_red) max_red = z_red;
      return max_red;
    }
    case CYLINDER:
      return std::max(x_red, y_red + z_red);
    case ELLIPSOID:
      return x_red + y_red + z_red;
    default:
    {
      ASSERT(cavity != INVALID_CAV)
      return 0;
    }
  }
}

// copies nodal values for the specified element id from the specified vector into p1, p2, p3
void ShellSystem::get_nodal_props(const ScalarVector& vector, unsigned int eid, Point& p1, Point& p2, Point& p3)
{
  const Elem* elem = get_mesh().elem(eid);
  for (unsigned int var = 0; var < 3; ++var)
  {
    p1(var) = vector[elem->get_node(0)->dof_number(0, var, 0)];
    p2(var) = vector[elem->get_node(1)->dof_number(0, var, 0)];
    p3(var) = vector[elem->get_node(2)->dof_number(0, var, 0)];
  }
}

// returns true if two elements are topologically close enough to need orientation testing
bool ShellSystem::elems_are_close(unsigned int eid1, unsigned int eid2)
{
  const Elem* e2 = get_mesh().elem(eid2);
  return (elem_cache[eid1]->close_elems.find(e2) != elem_cache[eid1]->close_elems.end());
}

// calculates the normal force of a contact
Real ShellSystem::shell_shell_contact_force(Real stiffness, Real distance)
{
  libmesh_assert(distance > 0);
  libmesh_assert(stiffness > 0);
  return stiffness * (thickness - distance) * 0.5 * (1 + thickness / distance);
}

// calculates the potential energy of a contact
Real ShellSystem::shell_shell_contact_energy(Real stiffness, Real distance)
{
  libmesh_assert(distance > 0);
  libmesh_assert(stiffness > 0);
  const Real depth = thickness - distance;
  return 0.5 * stiffness * (0.5 * depth * depth - depth * thickness + thickness * thickness * std::log(thickness / distance));
}

// calculates the normal force of a contact
Real ShellSystem::contact_force(Real stiffness, Real depth)
{
  libmesh_assert(depth >= 0);
  libmesh_assert(stiffness > 0);
  return stiffness * depth;
}

// calculates the potential energy of a contact
Real ShellSystem::contact_energy(Real stiffness, Real depth)
{
  libmesh_assert(depth >= 0);
  libmesh_assert(stiffness > 0);
  return 0.5 * stiffness * depth * depth;
}

// computes the matrix representing a cross product from the left
void ShellSystem::cross_product_matrix(RealTensorValue& dest, const RealVectorValue& z)
{
  libmesh_assert(dest(0,0) == 0);
  libmesh_assert(dest(1,1) == 0);
  libmesh_assert(dest(2,2) == 0);
                      dest(0,1) = -z(2);  dest(0,2) =  z(1);
  dest(1,0) =  z(2);                      dest(1,2) = -z(0);
  dest(2,0) = -z(1);  dest(2,1) =  z(0);
}

// returns the inverse of the 3x3 matrix A
RealTensorValue ShellSystem::inverse(const RealTensorValue& A)
{
  const Real det00 = A(1,1)*A(2,2) - A(1,2)*A(2,1);
  const Real det01 = A(1,2)*A(2,0) - A(1,0)*A(2,2);
  const Real det02 = A(1,0)*A(2,1) - A(1,1)*A(2,0);

  const Real det10 = A(0,2)*A(2,1) - A(0,1)*A(2,2);
  const Real det11 = A(0,0)*A(2,2) - A(0,2)*A(2,0);
  const Real det12 = A(0,1)*A(2,0) - A(0,0)*A(2,1);

  const Real det20 = A(0,1)*A(1,2) - A(0,2)*A(1,1);
  const Real det21 = A(0,2)*A(1,0) - A(0,0)*A(1,2);
  const Real det22 = A(0,0)*A(1,1) - A(0,1)*A(1,0);

  const Real inv_det = 1 / (A(0,0) * det00 + A(0,1) * det01 + A(0,2) * det02);
  libmesh_assert(std::isfinite(inv_det));

  return RealTensorValue(inv_det * det00, inv_det * det10, inv_det * det20,
                         inv_det * det01, inv_det * det11, inv_det * det21,
                         inv_det * det02, inv_det * det12, inv_det * det22);
}

// prints the table header for print_properties()
void ShellSystem::print_header()
{
  *os
    << "  timestamp    step       dt          time        wall_time "
    << "area        volume      density     pressure    Nc_cav      Nc_self     "
    << "E_kin       E_tens      E_bend      E_cav       E_self      E_fric      kGauss      kMean      \n"
    << "------------------------------------------------------------"
    << "------------------------------------------------------------------------"
    << "-----------------------------------------------------------------------------------------------"
    << std::endl;
}

#endif
