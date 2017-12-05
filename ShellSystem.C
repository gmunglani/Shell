#include "ShellSystem.h"
#include "SmoothShellShellContact.h"
#include "Friction.h"
#include "tools.h"

#include "mesh_subdiv_support.h"
#include "dof_map_subdiv.h"
#include "mesh_tools.h"
#include "face_tri3_sd.h"

#include <ctime>

#include "fe_base.h"
#include <sstream>

std::fstream file, file2;
int press[4];
Real fulltime;
Real volume_growth_red;

const RealTensorValue ShellSystem::Id(1,0,0,0,1,0,0,0,1); // identity matrix

ShellSystem::ShellSystem(EquationSystems& es, const std::string& name, const unsigned int number)
: System(es, name, number),

#ifdef MEASURE_PERFORMANCE
  perf_log("ShellSystem"),
#endif

  os(&std::cout),
  current_problem(NO_PROBLEM),
  last_problem(NO_PROBLEM),
  
  // data vectors
  u(3 * get_mesh().n_nodes(), 0), // displacement
  v(u.size(), 0),                 // velocity
  a(u.size(), 0),                 // acceleration
  f(u.size(), 0),                 // force
  Minv(u.size(), 0),              // inverse masses
  C(u.size(), 0),                 // viscous damping coefficients
  u_pred(u.size(), 0),            // predicted displacement
  v_pred(u.size(), 0),            // predicted velocity
  a_pred(u.size(), 0),            // predicted acceleration

  // geometric properties
  cavity(INVALID_CAV),

  // boundary condition properties
  shell_boundary(INVALID_BC),
  bc_coords(7), // 1+2+4 = 7
  bc_every(0),

  // dynamic properties
  gravity_dir(2), // by default, apply gravity in z direction
  pressure_base(INVALID_CS), // base in which the pressure is applied
  pressure_dir(2),
  growth_base(INVALID_CS), // base in which the growth tensor is diagonal
  grow_mass(true),
  area(0),
  volume(0),
  volume0(volume),
  stand_still(false),

  // contact properties
  shell_cavity_contact(true),
  friction_stop_velocity(0.03), // velocity below which static friction sets in
  static_friction_tol(0.1), // how far apart two contacting objects can be (multiple of thickness) for keeping the spring alive
  neighbor_margin(0.4), // the additional size of the topological neighborhood (multiple of thickness)

  // time integration
  dt(1e-6),
  newmark_beta(0.25),
  newmark_gamma(0.5),
  stop_on_problem(false),

  // adaptive stepsize control
  adapt_dt_every(1),
  dt_min(0),
  dt_max(0),
  rel_loc_err_max(3e-4),
  rel_loc_err_min(3e-5),

  // numeric stuff
  distortion(1e-4), // the fraction of the shell thickness for initial distortion

  // linked cell stuff
  du_sum(0),
  cell_size_mult(4),
  cell_margin(0.02),
  avg_elem_size(0), // to be set when the first time compute_bounding_box() gets called

  // other stuff
  time_step(0),
  extra_order(0), // use a single quadrature point per element per default
  density(0),
  n_ghost_nodes(0),
  n_ghost_elem(0),
  wall_time(0)
{
  clock_gettime(CLOCK_MONOTONIC, &start_time);
  
  // this system must be the first
  ASSERT(number == 0)

  // let the contact class know where to get the data from
  ShellShellContact::system = this;

  // subdivision elements have function support on the 1-neighborhood,
  // so we need a custom Dof Map
  replace_dof_map(AutoPtr<DofMap>(new DofMapSubdiv(number)));

  add_variable("ux", FOURTH, SUBDIV); // displacement in x direction
  add_variable("uy", FOURTH, SUBDIV); // displacement in y direction
  add_variable("uz", FOURTH, SUBDIV); // displacement in z direction
  ASSERT(variable_number("ux") == UX_VAR);
  ASSERT(variable_number("uy") == UY_VAR);
  ASSERT(variable_number("uz") == UZ_VAR);

  set_newmark_parameters(dt);
  
  // register the measurement vectors
  elem_vectors["kGauss"] = &kGauss_vec;
  elem_vectors["kMean" ] = &kMean_vec;
  elem_vectors["etens" ] = &etens_vec;
  elem_vectors["ebend" ] = &ebend_vec;
  node_vectors["ekin"  ] = &ekin_vec;
  node_vectors["ecav"  ] = &ecav_vec;
  node_vectors["eself" ] = &eself_vec;
  node_vectors["efric" ] = &efric_vec;
  node_vectors["mass"  ] = &mass_vec;
  
  // resize all measurement vectors
  std::map<std::string, ScalarVector*>::iterator vit = elem_vectors.begin();
  std::map<std::string, ScalarVector*>::iterator end_vit = elem_vectors.end();
  for (; vit != end_vit; ++vit)
    vit->second->resize(get_mesh().n_elem());
  
  vit = node_vectors.begin();
  end_vit = node_vectors.end();
  for (; vit != end_vit; ++vit)
    vit->second->resize(get_mesh().n_nodes());
}

ShellSystem::~ShellSystem()
{
  std::map<std::string, Function*>::iterator func = functions.begin();
  for (; func != functions.end(); ++func)
    delete func->second;

  for (unsigned int eid = 0; eid < elem_cache.size(); ++eid)
    delete elem_cache[eid];

  std::map<unsigned int, Real*>::iterator phi_it = phi_all.begin();
  for (; phi_it != phi_all.end(); ++phi_it)
    delete [] phi_it->second;

  std::map<unsigned int, RealGradient*>::iterator dphi_it = dphi_all.begin();
  for (; dphi_it != dphi_all.end(); ++dphi_it)
    delete [] dphi_it->second;

  std::map<unsigned int, RealTensor*>::iterator d2phi_it = d2phi_all.begin();
  for (; d2phi_it != d2phi_all.end(); ++d2phi_it)
    delete [] d2phi_it->second;

  for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
    delete shell_shell_contact_candidates[i];
}

// sets system parameters by parsing the config file and the command line arguments
void ShellSystem::set_parameters(GetPot& config_file, GetPot& command_line)
{
  PARSE(extra_order) // > 0 if more than a single quadrature point per element is desired
  ASSERT(extra_order >= 0)

  // material properties
  MAKE_FUNCTION(shell_mass_density , 1e3); // mass density
  MAKE_FUNCTION(shell_young_modulus, 1e9); // Young's modulus
  MAKE_FUNCTION(shell_young_modulus2, 1e9); // Young's modulus2
  MAKE_FUNCTION(shell_poisson_ratio, 1/3); // Poisson ratio
  MAKE_FUNCTION(shell_shear_modulus, 0); // Shell shear modulus
  MAKE_FUNCTION(bulk_modulus, 0); // bulk modulus for enclosed volume
  MAKE_FUNCTION(memb_stiff_mult, 1);
  MAKE_FUNCTION(bend_stiff_mult, 1);

  // geometric properties
  MAKE_FUNCTION(thickness, 1e-3);

  STR_PARSE(cavity, ELLIPSOID);
  STR_TO_ENUM(cavity, BOX);
  STR_TO_ENUM(cavity, CYLINDER);
  STR_TO_ENUM(cavity, ELLIPSOID);
  ASSERT(cavity != INVALID_CAV)

  MAKE_FUNCTION(Rx, 0.1005);
  MAKE_FUNCTION(Ry, 0.1005);
  MAKE_FUNCTION(Rz, 0.1005);

  // boundary condition properties
  STR_PARSE(shell_boundary, FREE);
  STR_TO_ENUM(shell_boundary, CLAMPED);
  STR_TO_ENUM(shell_boundary, CLAMPED_WEAK);
  STR_TO_ENUM(shell_boundary, PINNED);
  STR_TO_ENUM(shell_boundary, PINNED_WEAK);
  STR_TO_ENUM(shell_boundary, FREE);
  STR_TO_ENUM(shell_boundary, FREE_WEAK);
  ASSERT(shell_boundary != INVALID_BC)

  PARSE(bc_coords)
  PARSE(bc_every)
  MAKE_FUNCTION(bc_stiff_mult, 1);
  MAKE_FUNCTION(bc_criterion, 1);
  ASSERT(bc_coords <= 7)

  // dynamic properties
  MAKE_FUNCTION(damping_alpha, 0);
  MAKE_FUNCTION(damping_beta, 0);
  MAKE_FUNCTION(damping_gamma, 0);
  MAKE_FUNCTION(damping_rotation, 0);
  MAKE_FUNCTION(perturbation_rate, 0);
  MAKE_FUNCTION(perturbation_energy, 0);

  MAKE_FUNCTION(gravity, 0);
  PARSE(gravity_dir)
  MAKE_FUNCTION(pressure, 0);
  PARSE(pressure_dir)
  MAKE_FUNCTION(growth1, 1);
  MAKE_FUNCTION(growth2, 1);
  MAKE_FUNCTION(metric_growth_rate, 0);
  MAKE_FUNCTION(curvature_growth_rate, 0);
  PARSE(grow_mass)
  ASSERT(pressure_dir < 3)

  STR_PARSE(pressure_base, LOCAL);
  STR_TO_ENUM(pressure_base, LOCAL);
  STR_TO_ENUM(pressure_base, CARTESIAN);
  STR_TO_ENUM(pressure_base, CYLINDRICAL);
  STR_TO_ENUM(pressure_base, SPHERICAL);
  ASSERT(pressure_base != INVALID_CS)

  STR_PARSE(growth_base, LOCAL);
  STR_TO_ENUM(growth_base, LOCAL);
  STR_TO_ENUM(growth_base, CARTESIAN);
  STR_TO_ENUM(growth_base, CYLINDRICAL);
  ASSERT(growth_base != INVALID_CS)

  MAKE_PARSE(unsigned int, thickness_dir, 2)
  ASSERT(thickness_dir < 3)
  growth_dir[0] = (thickness_dir < 1 ? 1 : 0);
  growth_dir[1] = (thickness_dir < 2 ? 2 : 1);
  growth_dir[2] =  thickness_dir;

  MAKE_FUNCTION(volume_growth, 1);
  PARSE(stand_still)

  // contact-related parameters
  STR_PARSE(shell_shell_contact, NONE);
  STR_TO_ENUM(shell_shell_contact, NONE);
  STR_TO_ENUM(shell_shell_contact, LINEAR);
  STR_TO_ENUM(shell_shell_contact, SMOOTH);
  ASSERT(shell_shell_contact != INVALID_CT)
  
  PARSE(shell_cavity_contact) // whether the shell interacts with the cavity or not
  MAKE_FUNCTION(mu_s, 0);
  MAKE_FUNCTION(mu_d, 0);
  MAKE_FUNCTION(mu_s_cavity, 0);
  MAKE_FUNCTION(mu_d_cavity, 0);
  MAKE_FUNCTION(self_damping, 0);
  MAKE_FUNCTION(cavity_damping, 0);
  MAKE_FUNCTION(self_friction_damping, 0);
  MAKE_FUNCTION(cavity_friction_damping, 0);
  MAKE_FUNCTION(fric_stiff_mult, 1);
  PARSE(friction_stop_velocity)
  PARSE(static_friction_tol)
  PARSE(cell_size_mult)
  PARSE(cell_margin)
  PARSE(neighbor_margin) // multiple of the shell thickness to add to the topological neighborhood size

  // the expected maximum number of cells
  n_cells_max = static_cast<unsigned int>(std::ceil(std::sqrt(get_mesh().n_elem()) / cell_size_mult));
  n_cells_max = n_cells_max * n_cells_max * n_cells_max;

  ASSERT(friction_stop_velocity > 0)

  // time integration
  PARSE(dt)
  PARSE(newmark_beta)
  PARSE(newmark_gamma)
  MAKE_FUNCTION(stopping_criterion, 0);
  PARSE(stop_on_problem)
  ASSERT(dt > 0)
  ASSERT(newmark_beta > 0)

  // timestep adaptivity
  PARSE(adapt_dt_every)
  PARSE(dt_min)
  PARSE(dt_max)
  PARSE(rel_loc_err_max)
  PARSE(rel_loc_err_min)

  if (adapt_dt_every > 0)
  {
    ASSERT(rel_loc_err_min >= 0)
    ASSERT(rel_loc_err_max > rel_loc_err_min)
    ASSERT(newmark_beta != 1./6) // error estimation not possible for this value
  }

  // take as the target relative local error the logarithmic mean of the upper and lower bound
  rel_loc_err_target = std::sqrt(rel_loc_err_min * rel_loc_err_max);

  dt_min = std::max(dt_min, std::numeric_limits<Real>::epsilon());
  if (dt_max <= 0) dt_max = std::numeric_limits<Real>::max();
  ASSERT(dt_min <= dt)
  ASSERT(dt_max >= dt)

  // numeric stuff
  MAKE_PARSE(unsigned int, random_seed, std::time(0))
  srand48(random_seed);
  PARSE(distortion) // the fraction of the shell thickness for initial distortion

  // compute system parameters the first time, for assertions below
  std::map<std::string, Function*>::iterator func = functions.begin();
  for (; func != functions.end(); ++func)
    func->second->evaluate();

  // assert logical integrity of initial conditions of function parameters
  ASSERT(shell_mass_density > 0)
  ASSERT(shell_young_modulus > 0)
  ASSERT(shell_young_modulus2 > 0)
  ASSERT(shell_shear_modulus > 0)
  ASSERT(shell_poisson_ratio >= 0 && shell_poisson_ratio < 1)

  ASSERT(thickness > 0)
//  ASSERT(Rx > thickness)
//  ASSERT(Ry > thickness)
//  ASSERT(Rz > thickness)
  thickness0 = thickness;

  ASSERT(bc_stiff_mult > 0)

  ASSERT(growth1 > 0)
  ASSERT(growth2 > 0)
  ASSERT(metric_growth_rate >= 0)
  ASSERT(curvature_growth_rate >= 0)
  ASSERT(volume_growth > 0)
  
  const bool grow_statically = !FUNCTION_IS(growth1, 1) || !FUNCTION_IS(growth2, 1);
  grow_dynamically = !FUNCTION_IS(metric_growth_rate, 0) || !FUNCTION_IS(curvature_growth_rate, 0);

  // combination of growth types is not implemented
  ASSERT(!grow_statically || !grow_dynamically)

  // is any kind of in-plane growth enabled?
  grow_in_plane = grow_statically || grow_dynamically;

  ASSERT(damping_alpha >= 0)
  ASSERT(damping_beta >= 0)
  ASSERT(damping_gamma >= 0)
  ASSERT(damping_rotation >= 0)
  ASSERT(perturbation_rate >= 0)
  ASSERT(perturbation_energy >= 0)

  if (shell_shell_contact != NONE)
  {
    ASSERT(mu_s >= 0)
    ASSERT(mu_d >= 0)
    ASSERT(self_damping >= 0)
    ASSERT(self_friction_damping >= 0)
    ASSERT(cell_size_mult > 0)
    ASSERT(cell_margin >= 0)
    ASSERT(neighbor_margin > 0)
  }

  if (shell_cavity_contact)
  {
    ASSERT(mu_s_cavity >= 0)
    ASSERT(mu_d_cavity >= 0)
    ASSERT(cavity_damping >= 0)
    ASSERT(cavity_friction_damping >= 0)
    ASSERT(0 <= static_friction_tol)
  }

  self_friction   = shell_shell_contact != NONE && (!FUNCTION_IS(mu_s, 0) || !FUNCTION_IS(mu_d, 0));
  cavity_friction = shell_cavity_contact && (!FUNCTION_IS(mu_s_cavity, 0) || !FUNCTION_IS(mu_d_cavity, 0));
  
  // determine which quantities will need to be recomputed during the simulation
  const_derived = (functions["thickness"]->is_constant() &&
                   functions["Rx"]->is_constant() &&
                   functions["Ry"]->is_constant() &&
                   functions["Rz"]->is_constant() &&
                   functions["shell_mass_density" ]->is_constant() &&
                   functions["shell_young_modulus"]->is_constant() &&
                   functions["shell_young_modulus2"]->is_constant() &&
                   functions["shell_poisson_ratio"]->is_constant() &&
                   functions["shell_shear_modulus"]->is_constant() &&
                   functions["memb_stiff_mult"]->is_constant() &&
                   functions["bend_stiff_mult"]->is_constant() &&
                   functions["fric_stiff_mult"]->is_constant() &&
                   functions[  "bc_stiff_mult"]->is_constant());
  
  const_mass_damping = (functions["thickness"]->is_constant() &&
                        functions["shell_mass_density" ]->is_constant() &&
                        functions["shell_young_modulus"]->is_constant() &&
                        functions["shell_young_modulus2"]->is_constant() &&
   			functions["shell_poisson_ratio"]->is_constant() &&
   			functions["shell_shear_modulus"]->is_constant() &&
                        functions["memb_stiff_mult"]->is_constant() &&
                        functions["bend_stiff_mult"]->is_constant() &&
                        functions["damping_alpha"  ]->is_constant() &&
                        functions["damping_beta"   ]->is_constant() &&
                        functions["damping_gamma"  ]->is_constant() &&
                        !(grow_in_plane && grow_mass));

  set_newmark_parameters(dt);
  update_derived();
  build_cache();
  // compute M and C
  compute_mass_damping();
}

// calculates/updates some properties that are depending on others
void ShellSystem::update_derived()
{
  // effective semi-axes for cavity interaction
  const Real h_half = 0.5 * thickness;
  Rx_eff = Rx - h_half;
  Ry_eff = Ry - h_half;
  Rz_eff = Rz - h_half;

  // the average of them, for a fast approximation of the indentation depth
  avg_R_eff = (Rx_eff + Ry_eff + Rz_eff) / 3;

  // precompute inverse squared effective cavity semiaxes for efficieny
  inv_R2_eff[0] = 1 / (Rx_eff * Rx_eff);
  inv_R2_eff[1] = 1 / (Ry_eff * Ry_eff);
  inv_R2_eff[2] = 1 / (Rz_eff * Rz_eff);

  d2_touch = thickness * thickness; // touching distance squared
  d2_friction = std::pow((1 + static_friction_tol) * thickness, 2); // touching distance squared including friction scope
  one_plus_delta_friction = std::pow(1 - static_friction_tol * thickness / avg_R_eff, 2);

  switch (cavity)
  {
    case BOX:
    {
      cavity_surface = 8 * (Rx*Ry + Rx*Rz + Ry*Rz);
      cavity_volume = 8 * Rx * Ry * Rz;
      break;
    }
    case CYLINDER:
    {
      // approximation for the circumference of an ellipse from http://de.wikipedia.org/wiki/Ellipse#Umfang
      // (the symmetry axis of the cylinder is x)
      const Real lambda = (Ry - Rz) / (Ry + Rz);
      const Real lambda32 = 3 * lambda * lambda;
      cavity_surface = 2 * pi * ((Ry + Rz) * (1 + lambda32 / (10 + std::sqrt(4 - lambda32))) + Ry * Rz);
      cavity_volume  = 2 * pi * Rx * Ry * Rz;
      break;
    }
    case ELLIPSOID:
    {
      // Knud Thomsen's approximation formula for the surface area of an arbitrary ellipsoid
      cavity_surface = 4 * pi * std::pow((std::pow(Rx * Ry, 1.6)
                                        + std::pow(Rx * Rz, 1.6)
                                        + std::pow(Ry * Rz, 1.6)) / 3, 0.625);
      cavity_volume = 4./3 * pi * Rx * Ry * Rz;
      break;
    }
    default:
      ASSERT(cavity != INVALID_CAV)
  }

  // stiffness constants for the thin shell theory
  const Real shell_poisson_ratio2 = (shell_young_modulus2 * shell_poisson_ratio) / shell_young_modulus;
  const Real one_nu2 = 1 - shell_poisson_ratio * shell_poisson_ratio2;
  K = memb_stiff_mult * thickness / one_nu2; // membrane stiffness
  D = bend_stiff_mult * std::pow(thickness, 3) / (12 * one_nu2); // bending stiffness
  const Real shell_young_modulus_avg = (shell_young_modulus+shell_young_modulus2) * 0.5; 

  // weak boundary condition stiffness
  bc_stiff = bc_stiff_mult * shell_young_modulus_avg  * thickness / one_nu2;

  // Hertz force density, to be multiplied by the element area
  stiff_density = 0.5 * shell_young_modulus_avg / (one_nu2 * thickness);

  // stiffness of static friction springs
  stiff_static_friction = fric_stiff_mult * shell_young_modulus_avg * thickness / one_nu2;
}

// caches the dof incides for each element, to avoid recalculation
void ShellSystem::build_cache()
{
  const MeshBase& mesh = get_mesh();
  const DofMap& dof_map = get_dof_map();
  
  node_cache.resize(mesh.n_nodes());
  elem_cache.resize(mesh.n_elem(), NULL);
 
  // initialize the quadrature rule using element 0
  FEType fe_type = variable_type(UX_VAR);
  AutoPtr<QBase> qrule = fe_type.default_quadrature_rule(2, extra_order);
  AutoPtr<FEBase> fe(FEBase::build(2, fe_type));
  fe->attach_quadrature_rule(qrule.get());
  fe->reinit(mesh.elem(0));
  qw = qrule->get_weights();
  n_qp = qw.size();
  
  // in the loop below, we also count the number of ghost elements
  unsigned int n_real_elem = 0;
#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:n_real_elem)
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    // cache non-ghost finite elements only
    const Elem* elem = mesh.elem(eid);
    if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;
    
    ++n_real_elem;

    // create a cache entry
    elem_cache[eid] = new ElemCache;
    ElemCache& cache = *elem_cache[eid];

    // cache dof maps
    dof_map.dof_indices(elem, cache.dof_indices);
    dof_map.dof_indices(elem, cache.dof_indices_ux, UX_VAR);
    dof_map.dof_indices(elem, cache.dof_indices_uy, UY_VAR);
    dof_map.dof_indices(elem, cache.dof_indices_uz, UZ_VAR);

    // allocate enough space to hold the element caches
    cache.JxW   = new Real           [n_qp];
    cache.a_ref = new SurfaceMetric  [n_qp];
    cache.a_def = new SurfaceMetric  [n_qp];
    cache.H     = new RealTensorValue[n_qp];
    cache.alpha = new RealVectorValue[n_qp];
    cache.beta  = new RealVectorValue[n_qp];

    cache.qxyz = new Point           [n_qp];
    if (grow_in_plane)
    {
      // identity in-plane growth tensors
      cache.G = new RealVectorValue[n_qp];
      std::fill_n(cache.G, n_qp, RealVectorValue(1,1,0));
      
      if (grow_dynamically)
        cache.G_old = new RealVectorValue[n_qp];
      
      cache.GV  = new RealTensorValue[n_qp];
    }
    // we will need the 1-ring if limit surface contact is activated
    MeshTools::Subdiv::find_one_ring(static_cast<const Tri3SD*>(elem), cache.patch);
    
    // initialize the finite element
    AutoPtr<FEBase> fe(FEBase::build(2, fe_type));
    
#ifdef SHELL_OMP
#pragma omp critical(shell_fe)
#endif
    {
      // these are not thread-safe
      fe->attach_quadrature_rule(qrule.get());
      fe->reinit(elem);
      reinit_mod(elem, fe, qw);
	  //std::cout << "**************************************************" << std::endl;
	  //if (eid == 1) EXIT;	
	}
    // precalculate stuff for the reference configuration
    compute_reference_metrics(eid, fe);

    // here, the finite element goes out of scope and it is destroyed
  }

  // count the number of ghost elements
  n_ghost_elem = mesh.n_elem() - n_real_elem;

  // precompute maps and which nodes are ghost nodes
  MeshTools::build_nodes_to_elem_map(mesh, nodes_to_elem_map);
  Real n_ghost_nodes_sum = 0;
#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:n_ghost_nodes_sum)
#endif
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    std::vector<const Node*>& neighbors = node_cache[nid].neighbors;
    MeshTools::find_nodal_neighbors(mesh, mesh.node(nid), nodes_to_elem_map, neighbors);
    
    node_cache[nid].is_ghost = true;
    const std::vector<const Elem*>& elems = nodes_to_elem_map[nid];
    for (unsigned int eid = 0; eid < elems.size(); ++eid)
    {
      if (!static_cast<const Tri3SD*>(elems[eid])->is_ghost())
      {
        node_cache[nid].is_ghost = false;
        break;
      }
    }
    
    n_ghost_nodes_sum += node_cache[nid].is_ghost;
  }
  n_ghost_nodes = n_ghost_nodes_sum;
  
  // determine the set of node ids involved in boundary conditions
  if (n_ghost_nodes > 0 && (bc_coords & 7) && bc_every > 0)
  {
    // iterate over all ghost elements
#ifdef SHELL_OMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      const Elem* elem = mesh.elem(eid);
      const Tri3SD* gh_elem = static_cast<const Tri3SD*>(elem);
      if (!gh_elem->is_ghost()) continue;

      // for all sides on which a non-ghost element lies
      for (unsigned int s = 0; s < elem->n_sides(); ++s)
      {
        const Tri3SD* nb_elem = static_cast<const Tri3SD*>(elem->neighbor(s));
        if (nb_elem == NULL || nb_elem->is_ghost()) continue;

/*
    n4
   /  \
  / gh \
n2 ---- n3
  \ nb /
   \  /
    n1
*/

        unsigned int nids [4]; // n1, n2, n3, n4
        nids[1] = gh_elem->node(s); // n2
        nids[2] = gh_elem->node(MeshTools::Subdiv::next[s]); // n3
        nids[3] = gh_elem->node(MeshTools::Subdiv::prev[s]); // n4

        // find n1
        nids[0] = nb_elem->node(0);
        for (unsigned int n = 1; nids[0] == nids[1] || nids[0] == nids[2]; ++n)
          nids[0] = nb_elem->node(n);

#ifdef SHELL_OMP
#pragma omp critical(shell_add_boundary)
#endif
        {
          // add this boundary edge only if its middle point meets the BC criterion
          xyz = 0.5 * (mesh.node(nids[1]) + mesh.node(nids[2]));
          functions["bc_criterion"]->evaluate();
          if (bc_criterion == 1)
            boundary_cache.push_back(BoundaryCache(nids));
        }
      }
    }
  }

  if (shell_shell_contact == NONE) return; // everything below is needed only for shell-shell contacts

  // determine the direct topological neighbors for all non-ghost elements
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    const Elem* elem = mesh.elem(eid);
    if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;

    std::set<const Elem*>& neighbor_elems = elem_cache[eid]->neighbor_elems;

    for (unsigned int n = 0; n < elem->n_nodes(); ++n)
    {
      std::vector<const Elem*>& node_elements = nodes_to_elem_map[elem->node(n)];
      for (unsigned int i = 0; i < node_elements.size(); ++i)
      {
        // add the neighbor element unless it is the element itself or a ghost
        if (node_elements[i]->id() != eid && !static_cast<const Tri3SD*>(node_elements[i])->is_ghost())
          neighbor_elems.insert(node_elements[i]);
      }
    }
  }
  
  build_topological_neighborhood();
}

// determines the set of topologically close elements
void ShellSystem::build_topological_neighborhood()
{
  PERFLOG_START("build_topological_neighborhood()");

  const MeshBase& mesh = get_mesh();
  
  // neighbor_margin should be chosen large if the shell
  // is expected to be much complessed in plane
  const Real cutoff = (1 + neighbor_margin) * thickness / std::min(growth1, growth2);
  const Real cutoff2 = cutoff * cutoff;

#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

    std::set<const Elem*>& close_elems = elem_cache[eid]->close_elems;
    close_elems.clear();
    const std::set<const Elem*>& neighbor_elems = elem_cache[eid]->neighbor_elems;
    std::set<const Elem*> added = close_elems = neighbor_elems;

    do // test the neighbor's neighbors until none of them are close enough any more
    {
      // add the neighbors of the previously added elements to the set of new candidates
      std::set<const Elem*> candidates;
      std::set<const Elem*>::iterator a_it = added.begin();
      for (; a_it != added.end(); ++a_it)
      {
        const std::set<const Elem*>& neighbors_of_added = elem_cache[(*a_it)->id()]->neighbor_elems;
        std::set<const Elem*>::iterator n_it = neighbors_of_added.begin();
        const std::set<const Elem*>::iterator n_end = neighbors_of_added.end();
        for (; n_it != n_end; ++n_it)
          if ((*n_it)->id() != eid) // unless it is the element itself
            candidates.insert(*n_it);
      }

      // forget the previously added elements
      added.clear();

      // add a candidate if it is close enough and hasn't been added yet
      std::set<const Elem*>::iterator c_it = candidates.begin();
      for (; c_it != candidates.end(); ++c_it)
      {
        ShellShellContact contact(eid, (*c_it)->id());
        contact.compute_on_mesh();
        if (contact.d2 < cutoff2 && close_elems.find(*c_it) == close_elems.end())
        {
          added.insert(*c_it);
          close_elems.insert(*c_it);
        }
      }
    } while (added.size() > 0);
  }

  PERFLOG_STOP("build_topological_neighborhood()");
}

// sets the newmark integration coefficients
void ShellSystem::set_newmark_parameters(Real new_dt)
{
  dt = std::max(new_dt, dt_min);
  dt = std::min(    dt, dt_max);

  newmark_a[0] = newmark_beta * dt * dt;
  newmark_a[1] = newmark_gamma / (newmark_beta * dt);
  newmark_a[2] = 1 / (newmark_beta * dt);
  newmark_a[3] = 1 / (2 * newmark_beta) - 1;
  newmark_a[4] = newmark_gamma / newmark_beta - 1;
  newmark_a[5] = dt * (newmark_gamma / (2 * newmark_beta) - 1);
  newmark_a[6] = (1 - newmark_gamma) * dt;
  newmark_a[7] = newmark_gamma * dt;
  newmark_a[8] = dt * dt * (0.5 - newmark_beta);
  newmark_a[9] = dt * dt * std::fabs(newmark_beta - 1./6);
}

// computes M and C
void ShellSystem::compute_mass_damping()
{
  const MeshBase& mesh = get_mesh();

  Minv.zero();
  C.zero();

#ifdef SHELL_OMP
#pragma omp parallel
#endif
  {
    // the element mass and damping matrices and submatrices (lumped to vectors)
    DenseVector<Real> Me, Ce;
    DenseSubVector<Real> Mx(Me), My(Me), Mz(Me);

    // iterate over all non-ghost elements
#ifdef SHELL_OMP
#pragma omp for
#endif
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      const Elem* elem = mesh.elem(eid);
      if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;

      const ElemCache& cache = *elem_cache[eid];
      const unsigned int n_var_dofs = cache.dof_indices_ux.size();

      Me.resize(cache.dof_indices.size());
      Mx.reposition(UX_VAR * n_var_dofs, n_var_dofs);
      My.reposition(UY_VAR * n_var_dofs, n_var_dofs);
      Mz.reposition(UZ_VAR * n_var_dofs, n_var_dofs);
      Ce.resize(Me.size());

      Real avg_growth = 0;
      for (unsigned int qp = 0; qp < n_qp; ++qp)
      {
        Real JxW;
        
        // account for change of mass during growth
        if (grow_in_plane && grow_mass)
        {
          const RealVectorValue& G = cache.G[qp];
          const Real det = G(0) * G(1) - 0.25 * G(2) * G(2); // G contains the Voigt factor 2
          JxW = cache.JxW[qp] * det;
          avg_growth += det;
        }
        else
        {
          JxW = cache.JxW[qp];
          avg_growth += 1;
        }
        
        for (unsigned int i = 0; i < n_var_dofs; ++i)
        {
          Real lumped_mass = 0;
          for (unsigned int j = 0; j < n_var_dofs; ++j)
            lumped_mass += cache.phi[j * n_qp + qp];
          lumped_mass *= JxW * cache.phi[i * n_qp + qp];
          Mx(i) += lumped_mass;
          My(i) += lumped_mass;
          Mz(i) += lumped_mass;
        }
      }
      avg_growth /= n_qp;
      
      Me.scale(shell_mass_density * (grow_mass ? thickness : thickness0));

      // viscous damping factors
      // Note: an alternative model would be C_i = crit_damp_frac * 2 * sqrt(M_i * K_i)
      const Real avg_stiff = 2 * (K + D / (avg_growth * elem->hmin() * elem->hmax())) / 3;
      for (unsigned int i = 0; i < Me.size(); ++i)
        Ce(i) = damping_alpha * Me(i) + damping_beta * avg_stiff;
     
//      Me.print(std::cout);
//      Ce.print(std::cout);
//	  EXIT
#ifdef SHELL_OMP
#pragma omp critical(shell_assemble_MC)
#endif
      {
        Minv.add_vector(Me, cache.dof_indices);
        C.add_vector   (Ce, cache.dof_indices);
      }
    }
  }
  
  C.add(damping_gamma);
  
  // invert the mass matrix so that Minv stores the inverse lumped masses
  Minv.reciprocal();
}

// applies initial conditions (using a loadfile, if given) and precomputes stuff
void ShellSystem::apply_initial_conditions(std::ifstream & loadfile)
{
  if (loadfile.is_open())
  {
    load(loadfile);
  }
  else
  {
    const MeshBase& mesh = get_mesh();
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      const Node& node = mesh.node(nid);
      for (unsigned int var = 0; var < 3; ++var)
        u[node.dof_number(0, var, 0)] += distortion * (2 * drand48() - 1) * thickness;
    }
  }

  // make sure the initial force calculation etc. uses the initial u
  u_pred = u;
  compute_position();
  
  reset_properties();

  compute_bounding_box(bounds, avg_elem_size);
  
  // precompute the the deformed surface metrics
  compute_deformed_metrics();

  // compute the enclosed volume
  compute_enclosed_volume();
  volume0 = volume;

  // always compute the lists of contact pairs at the beginning
  compute_shell_shell_contact_candidates(true);

  // compute the initial acceleration a = M^{-1}*f such that the error estimate
  // in adaptive time stepping is small at the beginning
  if (!loadfile.is_open())
  {
    const Real old_time = this->time;
    const unsigned int old_time_step = time_step;
    ++time_step;
    prepare_next_step();
    compute_external_force();
    subtract_internal_force();
//	EXIT
	a.pointwise_mult(Minv, f);
    this->time = old_time;
    time_step = old_time_step;
  }
}

// checks if the current state meets a stopping citerion, returning true if the simulation should be terminated
bool ShellSystem::check_stopping_criterion()
{
  functions["stopping_criterion"]->evaluate();
  
  if (stopping_criterion != 0)
  {
#ifdef MEASURE_PERFORMANCE
    perf_log.print_log();
#endif
    return true;
  }

  return false;
}

// prepares some stuff for the next timestep
void ShellSystem::prepare_next_step()
{
  reset_properties();
  
  // advance the time
  this->time += dt;

  const Real old_thickness = thickness;

  // update the system parameters, which can be functions of time, density, etc.
  std::map<std::string, Function*>::iterator func = functions.begin();
  for (; func != functions.end(); ++func)
    func->second->evaluate();

  if (thickness > old_thickness)
    du_sum += thickness - old_thickness;

  // recompute derived properties unless nothing relevant changed
  if (!const_derived) update_derived();

  // let the shell grow (if enabled)
  if (grow_in_plane) apply_growth();

  // recompute mass and damping matrices unless nothing relevant changed
  if (!const_mass_damping) compute_mass_damping();
}

// integrates one step in time
bool ShellSystem::advance(unsigned int next_time_step)
{
  PERFLOG_START("advance()");
  
  time_step = next_time_step;
  const bool adapt_dt_now = (adapt_dt_every && time_step % adapt_dt_every == 0);

  // repeat until the local error is small enough
  while (true)
  {
    prepare_next_step();

#ifdef SHELL_OMP
#pragma omp parallel sections
#endif
    {
#ifdef SHELL_OMP
#pragma omp section
#endif
      {
        // predict the displacement
        u_pred = u;
        u_pred.add(dt, v);
        u_pred.add(newmark_a[8], a);
    
        // now that we have a new u_pred, compute the total nodal positions
        check_increment(false);
        compute_position();
      }
#ifdef SHELL_OMP
#pragma omp section
#endif
      {
        // predict the velocity
        v_pred = v;
        v_pred.add(newmark_a[6], a);
        perturb();
      }
    }
    
    // update the contact lists if necessary
    compute_shell_shell_contact_candidates();
    
    // compute the the predicted deformed surface metrics
    compute_deformed_metrics();

    // compute the generalized forces
    compute_external_force();
    subtract_internal_force();

    // apply boundary conditions (weak ones are external forces)
    apply_boundary_conditions();

    // a_pred is unused here, so we can use it as a temporary for the damping force
    a_pred.pointwise_mult(C, v_pred);
    f -= a_pred;
    a_pred.pointwise_mult(Minv, f);

    Real rel_loc_err = rel_loc_err_target;
    Real new_dt = dt;
    if (adapt_dt_now)
    {
      // Compute the relative local error according to Zienkiewicz and Xie
      // with the modification that the reference length is not the recorded
      // maximum displacement but the shell thickness. We can use f as a
      // temporary object because it is no longer needed below.
      f = a_pred;
      f -= a;
      const Real norm_a = f.linfty_norm();
      
      if (std::isfinite(norm_a))
        rel_loc_err = newmark_a[9] * norm_a / thickness;
      else
      {
        register_problem(NOT_FINITE);
        problem_info.set<Real>("norm_a") = norm_a;
      }
      
      if (rel_loc_err == 0)
        new_dt = 2 * dt;
      else if (rel_loc_err > rel_loc_err_max || rel_loc_err < rel_loc_err_min)
        new_dt = dt * std::pow(rel_loc_err_target/rel_loc_err, 1./3);

      // if any type of contact is enabled, make sure the time step is small
      // enough to resolve the contact depth with at least 10 steps
      if (shell_cavity_contact || shell_shell_contact != NONE)
      {
        static const unsigned int contact_resolution = 10;
        const Real vmax = v_pred.linfty_norm();
        if (new_dt * vmax * contact_resolution * 2 > thickness)
          new_dt = thickness / (2 * contact_resolution * vmax);
      }
    }
    
    if (new_dt >= dt)
    {
      // apply hard boundary conditions (weak ones have no effect here)
      // Note: the fact that hard BCs also affect a_pred allows this to be
      // called before the displacement and velocity corrections below
      apply_boundary_conditions();
      
      if (stand_still)
      {
        a_pred.zero();
        v_pred.zero();
      }
      
#ifdef SHELL_OMP
#pragma omp parallel sections
#endif
      {
#ifdef SHELL_OMP
#pragma omp section
#endif
        {
          // correct the displacement
          u_pred.add(newmark_a[0], a_pred);
          check_increment();
          
          // now that we have a new u_pred, compute the total nodal positions
          compute_position();

          // accept the new displacement
          u.swap(u_pred);
        }
#ifdef SHELL_OMP
#pragma omp section
#endif
        {
          // correct the velocity
          v_pred.add(newmark_a[7], a_pred);
          
          // accept the new velocity
          v.swap(v_pred);
        }
      }

      // accept the new acceleration
      a.swap(a_pred);

      // update the contact lists if necessary
      compute_shell_shell_contact_candidates();
      
      // measure some properties of the new configuration
      // (must be before setting the new dt to have the old dt available for output)
      measure_properties();

      set_newmark_parameters(new_dt);

      break; // exit the time step repetition loop
    }
    
    // if we made it through here, we will need to repeat the time step

    if (grow_dynamically) reject_growth();
    this->time -= dt;
    set_newmark_parameters(new_dt);

    // going below the machine precision doesn't make sense
    if (dt <= std::numeric_limits<Real>::epsilon())
    {
      register_problem(SMALL_TIMESTEP);
      problem_info.set<Real>("dt") = dt;
      break;
    }
  }

  PERFLOG_STOP("advance()");
  
  return check_for_problems();
}

// subtracts the internal force vector from f
void ShellSystem::subtract_internal_force()
{
  PERFLOG_START("subtract_internal_force()");
  
  const MeshBase& mesh = get_mesh();

  Real Etens_sum = 0, Ebend_sum = 0;
#ifdef SHELL_OMP
#pragma omp parallel reduction(+:Etens_sum,Ebend_sum)
#endif
  {
    DenseVector<Real> qe; // element internal force
    DenseSubVector<Real> qx(qe), qy(qe), qz(qe); // element internal force contributions from each variable
    RealTensorValue BmT, BbT; // membrane and bending strain matrices (T for transposed)
    RealVectorValue n, m; // membrane and bending stress resultants in Voigt notation
    RealVectorValue p11, p21, p12, p22, p13, p23; // intermediate parameters
    RealVectorValue column, a2xa3, a3xa1, fint_density;
    RealTensorValue GVT; // transposed inverse growth tensor in Voigt notation

    // iterate over all non-ghost elements
#ifdef SHELL_OMP
#pragma omp for
#endif
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

      ElemCache& cache = *elem_cache[eid];

      const std::vector<unsigned int>& dof_inds = cache.dof_indices;
      const unsigned int n_var_dofs = cache.dof_indices_ux.size();

      qe.resize(dof_inds.size());
      qx.reposition(UX_VAR * n_var_dofs, n_var_dofs);
      qy.reposition(UY_VAR * n_var_dofs, n_var_dofs);
      qz.reposition(UZ_VAR * n_var_dofs, n_var_dofs);

      Real Etens_elem = 0, Ebend_elem = 0;

      for (unsigned int qp = 0; qp < n_qp; ++qp)
      {
        const Real& JxW = cache.JxW[qp];
        
        // the surface metrics
        const SurfaceMetric& a_ref = cache.a_ref[qp];
        const SurfaceMetric& a_def = cache.a_def[qp];

        // membrane strain in Voigt notation, eqs. (34,15)
        // bending strain in Voigt notation, eqs. (34,19), not (16)! (different sign)
        RealVectorValue& alpha = cache.alpha[qp];
        RealVectorValue& beta  = cache.beta[qp];
        if (grow_in_plane)
        {
          alpha = 0.5 * (cache.GV[qp] * a_def.ff1 - a_ref.ff1);
          beta  = a_ref.ff2 - cache.GV[qp] * a_def.ff2;
          
          // also precompute the Voigt-corrected inverse transposed growth tensor
          GVT = cache.GV[qp].transpose();
          GVT(0,2) *= 2;
          GVT(1,2) *= 2;
          GVT(2,0) *= 0.5;
          GVT(2,1) *= 0.5;
        }
        else
        {
          alpha = 0.5 * (a_def.ff1 - a_ref.ff1);
          beta  = a_ref.ff2 - a_def.ff2;
        }
        alpha(2) *= 2;
        beta (2) *= 2;

        RealTensorValue& H = cache.H[qp];
	
		// poisson ratio in direction 2
        const Real shell_poisson_ratio2 = (shell_young_modulus * shell_poisson_ratio) / shell_young_modulus2;

	    // transformed elasticity tensor
        H(0,0)          = shell_young_modulus * a_ref.ff1(1) * a_ref.ff1(1);
        H(0,1) = H(1,0) = shell_young_modulus * (shell_poisson_ratio * a_ref.ff1(1) * a_ref.ff1(0)
        	                + (1 - shell_poisson_ratio) * a_ref.ff1(2) * a_ref.ff1(2));
        H(0,2) = H(2,0) = shell_young_modulus * -a_ref.ff1(1) * a_ref.ff1(2);
        H(1,1)          = shell_young_modulus2 * a_ref.ff1(0) * a_ref.ff1(0);
        H(1,2) = H(2,1) = shell_young_modulus2 * -a_ref.ff1(0) * a_ref.ff1(2);
        H(2,2)          = shell_shear_modulus * ((1 - shell_poisson_ratio * shell_poisson_ratio2) * a_ref.ff1(1) * a_ref.ff1(0)
                            + (1 - shell_poisson_ratio * shell_poisson_ratio2) * ((1 + shell_poisson_ratio) / (1 - shell_poisson_ratio)) * a_ref.ff1(2) * a_ref.ff1(2));

	    const Real invdet = 1 / (a_ref.ff1(0) * a_ref.ff1(1) - a_ref.ff1(2) * a_ref.ff1(2));
        libmesh_assert(std::isfinite(invdet)); 
        H *= invdet * invdet;	
	
		// stress resultants in Voigt notation
        n = K * (H * alpha); // eq. (35)
        m = D * (H * beta);  // eq. (36)

//		std::cout << "D*H " << D*H << std::endl;
//		std::cout << "Beta " << beta << std::endl;

		// precalculate parameters for bending strain tensor, eq. (A6)
        a2xa3 = a_def.vec[1].cross(a_def.vec[2]);
        a3xa1 = a_def.vec[2].cross(a_def.vec[0]);

        const Real invjac = 1 / a_def.jac;
        p11 = invjac * (a_def.dvec[0].cross(a_def.vec[1]) + a_def.ff2(0) * a2xa3);
        p12 = invjac * (a_def.dvec[1].cross(a_def.vec[1]) + a_def.ff2(1) * a2xa3);
        p13 = invjac * (a_def.dvec[2].cross(a_def.vec[1]) + a_def.ff2(2) * a2xa3);

        p21 = invjac * (a_def.vec[0].cross(a_def.dvec[0]) + a_def.ff2(0) * a3xa1);
        p22 = invjac * (a_def.vec[0].cross(a_def.dvec[1]) + a_def.ff2(1) * a3xa1);
        p23 = invjac * (a_def.vec[0].cross(a_def.dvec[2]) + a_def.ff2(2) * a3xa1);

        // iterate over all basis functions
        for (unsigned int i = 0; i < n_var_dofs; ++i)
        {
          const unsigned int idx = i * n_qp + qp;
          const RealGradient& dphi = cache.dphi [idx];
          const RealTensor&  d2phi = cache.d2phi[idx];

          // membrane part, eq. (A5) transposed
          column = a_def.vec[0] * dphi(0);
          BmT(0,0) = column(0);
          BmT(1,0) = column(1);
          BmT(2,0) = column(2);
          column = a_def.vec[1] * dphi(1);
          BmT(0,1) = column(0);
          BmT(1,1) = column(1);
          BmT(2,1) = column(2);
          column = a_def.vec[0] * dphi(1) + a_def.vec[1] * dphi(0);
          BmT(0,2) = column(0);
          BmT(1,2) = column(1);
          BmT(2,2) = column(2);

          // bending part, eq. (A6) transposed (note: Cirak is missing the factor 2!)
          column = a_def.vec[2] * -d2phi(0,0) + (dphi(0) * p11 + dphi(1) * p21);
          BbT(0,0) = column(0);
          BbT(1,0) = column(1);
          BbT(2,0) = column(2);
          column = a_def.vec[2] * -d2phi(1,1) + (dphi(0) * p12 + dphi(1) * p22);
          BbT(0,1) = column(0);
          BbT(1,1) = column(1);
          BbT(2,1) = column(2);
          column = 2 * (a_def.vec[2] * -d2phi(0,1) + (dphi(0) * p13 + dphi(1) * p23));
          BbT(0,2) = column(0);
          BbT(1,2) = column(1);
          BbT(2,2) = column(2);
          
          // growth correction
          if (grow_in_plane)
          {
            BmT = BmT * GVT;
            BbT = BbT * GVT;
          }

          // nonlinear form of eq. (45)
          fint_density = BmT * n + BbT * m;
          qx(i) += JxW * fint_density(0);
          qy(i) += JxW * fint_density(1);
          qz(i) += JxW * fint_density(2);
        }

        // measure elastic energies (non-variational form of eq. (38))
        Etens_elem += 0.5 * JxW * (n * alpha);
        Ebend_elem += 0.5 * JxW * (m * beta);
 
//		std::cout << "BmT and n " << BmT*n << std::endl;
//		std::cout << "BbT and m " << BbT*m << std::endl;
//		std::cout << "n and alpha " << n*alpha << std::endl;
//		std::cout << "m and beta " << m*beta << std::endl;
//		std::cout << "Ebend " << Ebend_elem << std::endl;
//		std::cout << "m " << m << std::endl;
//		std::cout << "Bmt " << BmT << std::endl;
     }
      
      Etens_sum += Etens_elem;
      Ebend_sum += Ebend_elem;
      etens_vec[eid] += Etens_elem / cache.area;
      ebend_vec[eid] += Ebend_elem / cache.area;

      // insert the element contribution into the vector
#ifdef SHELL_OMP
#pragma omp critical(shell_assemble_vec)
#endif
      f.subtract_vector(qe, dof_inds);
    }
  }
  Etens += Etens_sum;
  Ebend += Ebend_sum;
//  EXIT

  PERFLOG_STOP("subtract_internal_force()");
}

// computes the reference surface metrics and H matrices for the given element
void ShellSystem::compute_reference_metrics(const unsigned int eid, const AutoPtr<FEBase>& fe)
{
  const std::vector<Point>&        xyz          = fe->get_xyz();
  const std::vector<RealGradient>& dxyzdxi      = fe->get_dxyzdxi();
  const std::vector<RealGradient>& dxyzdeta     = fe->get_dxyzdeta();
  const std::vector<RealGradient>& d2xyzdxi2    = fe->get_d2xyzdxi2();
  const std::vector<RealGradient>& d2xyzdeta2   = fe->get_d2xyzdeta2();
  const std::vector<RealGradient>& d2xyzdxideta = fe->get_d2xyzdxideta();

  for (unsigned int qp = 0; qp < n_qp; ++qp)
  {
    SurfaceMetric& a_ref = elem_cache[eid]->a_ref[qp];
    RealTensorValue& H = elem_cache[eid]->H[qp];

    a_ref.xyz = elem_cache[eid]->qxyz[qp];

    // eq. (3) in Cirak et al.
    a_ref.vec[0]  = dxyzdxi[qp];
    a_ref.vec[1]  = dxyzdeta[qp];
    // eq. (9)
    a_ref.vec[2]  = a_ref.vec[0].cross(a_ref.vec[1]);
    a_ref.jac     = a_ref.vec[2].size();
    libmesh_assert(a_ref.jac != 0);
    a_ref.vec[2] *= 1 / a_ref.jac;

    // derivative
    a_ref.dvec[0] = d2xyzdxi2[qp];
    a_ref.dvec[1] = d2xyzdeta2[qp];
    a_ref.dvec[2] = d2xyzdxideta[qp];

    // first fundamental form of reference surface (eq. (4), left)
    a_ref.ff1(0) = a_ref.vec[0] * a_ref.vec[0];
    a_ref.ff1(1) = a_ref.vec[1] * a_ref.vec[1];
    a_ref.ff1(2) = a_ref.vec[0] * a_ref.vec[1];

    // second fundamental form of reference surface (first term of rhs of eq. (19))
    a_ref.ff2(0) = a_ref.dvec[0] * a_ref.vec[2];
    a_ref.ff2(1) = a_ref.dvec[1] * a_ref.vec[2];
    a_ref.ff2(2) = a_ref.dvec[2] * a_ref.vec[2];

//	if (eid < 2)
//	{
//		std::cout << "A_Ref.ff2 " << a_ref.vec[2] << std::endl;
//	    EXIT	
//	}
    // H matrix (eq. (37) via eq. (5) in Cirak)
    H(0,0)          = a_ref.ff1(1) * a_ref.ff1(1);
    H(0,1) = H(1,0) =    shell_poisson_ratio  * a_ref.ff1(1) * a_ref.ff1(0)
                    + (1-shell_poisson_ratio) * a_ref.ff1(2) * a_ref.ff1(2);
    H(0,2) = H(2,0) = -a_ref.ff1(1) * a_ref.ff1(2);
    H(1,1)          =  a_ref.ff1(0) * a_ref.ff1(0);
    H(1,2) = H(2,1) = -a_ref.ff1(0) * a_ref.ff1(2);
    H(2,2)          = 0.5 * ((1-shell_poisson_ratio) * a_ref.ff1(1) * a_ref.ff1(0)
                           + (1+shell_poisson_ratio) * a_ref.ff1(2) * a_ref.ff1(2));
    const Real invdet = 1 / (a_ref.ff1(0) * a_ref.ff1(1) - a_ref.ff1(2) * a_ref.ff1(2));
    libmesh_assert(std::isfinite(invdet));
    H *= invdet * invdet;
    
    // copy the reference metric to the initial deformed metric
    // so that the mass and damping matrices can be computed using the deformed state
    elem_cache[eid]->a_def[qp] = a_ref;
  }
}

// computes the deformed surface metrics based on the current displacement
void ShellSystem::compute_deformed_metrics()
{
  PERFLOG_START("compute_deformed_metrics()");

  const MeshBase& mesh = get_mesh();
  
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

    ElemCache& cache = *elem_cache[eid];
    cache.area = 0;

    const Real*          phi = cache.phi;
    const RealGradient* dphi = cache.dphi;
    const RealTensor*  d2phi = cache.d2phi;
    const std::vector<unsigned int>& dof_inds_ux = cache.dof_indices_ux;
    const std::vector<unsigned int>& dof_inds_uy = cache.dof_indices_uy;
    const std::vector<unsigned int>& dof_inds_uz = cache.dof_indices_uz;
    
    RealVectorValue u_qp; // current displacement at the quadrature point

    for (unsigned int qp = 0; qp < n_qp; ++qp)
    {
      const SurfaceMetric& a_ref = cache.a_ref[qp];
            SurfaceMetric& a_def = cache.a_def[qp];

      // find the current displacement approximant at the quadrature point
      RealGradient u, dudxi, dudeta, d2udxi2, d2udeta2, d2udxideta;
      for (unsigned int i = 0; i < dof_inds_ux.size(); ++i)
      {
        const unsigned int idx = i * n_qp + qp;
        u_qp(0) = u_pred[dof_inds_ux[i]];
        u_qp(1) = u_pred[dof_inds_uy[i]];
        u_qp(2) = u_pred[dof_inds_uz[i]];
        u.add_scaled         (u_qp, phi[idx]);
        dudxi.add_scaled     (u_qp, dphi[idx](0));
        dudeta.add_scaled    (u_qp, dphi[idx](1));
        d2udxi2.add_scaled   (u_qp, d2phi[idx](0,0));
        d2udeta2.add_scaled  (u_qp, d2phi[idx](1,1));
        d2udxideta.add_scaled(u_qp, d2phi[idx](0,1));
      }

      // compute the deformed surface metric
      a_def.xyz = a_ref.xyz + u;

      // eqs. (3,20)
      a_def.vec[0]  = a_ref.vec[0] + dudxi;
      a_def.vec[1]  = a_ref.vec[1] + dudeta;
      // eq. (18)
      a_def.vec[2]  = a_def.vec[0].cross(a_def.vec[1]);
      a_def.jac     = a_def.vec[2].size(); // Ciarlet (2000) proposed to use a_ref.jac here
      libmesh_assert(a_def.jac > 0);
      a_def.vec[2] *= 1 / a_def.jac;

      // derivative
      a_def.dvec[0] = a_ref.dvec[0] + d2udxi2;
      a_def.dvec[1] = a_ref.dvec[1] + d2udeta2;
      a_def.dvec[2] = a_ref.dvec[2] + d2udxideta;

      // first fundamental form of deformed surface (eq. (4), right)
      a_def.ff1(0) = a_def.vec[0] * a_def.vec[0];
      a_def.ff1(1) = a_def.vec[1] * a_def.vec[1];
      a_def.ff1(2) = a_def.vec[0] * a_def.vec[1];

      // second fundamental form of deformed surface (second term of rhs of eq. (19))
      a_def.ff2(0) = a_def.dvec[0] * a_def.vec[2];
      a_def.ff2(1) = a_def.dvec[1] * a_def.vec[2];
      a_def.ff2(2) = a_def.dvec[2] * a_def.vec[2];

//	  if (eid < 2)
//	  	  std::cout << "A_def.ff2 " << a_def.ff2 << std::endl;

	// add the contribution to the deformed element area
      cache.area += qw[qp] * a_def.jac;
    }
  }

  PERFLOG_STOP("compute_deformed_metrics()");
}

// computes all external forces
void ShellSystem::compute_external_force()
{
  PERFLOG_START("compute_external_force()");

  f.zero();
  const MeshBase& mesh = get_mesh();
  
  volume_growth_red = std::max((1.0 - static_cast<int>(time*5.0)*0.01),0.5);
  if (time_step == 2) volume0 = volume;	
  const Real volume_pressure = (bulk_modulus == 0 ? 0 : bulk_modulus * (volume0*volume_growth_red/volume - 1));
//  std::cout << "Volume_print " << volume_pressure << " " << volume0 << " " << volume << std::endl;
//  if (time_step > 20) EXIT  
  
  pressure *= fulltime;
  if (pressure != 0 || volume_pressure != 0)
  {
#ifdef SHELL_OMP
#pragma omp parallel
#endif
    {
      DenseVector<Real> qe; // element internal force
      DenseSubVector<Real> qx(qe), qy(qe), qz(qe); // element internal force contributions from each variable
      
      // iterate over all non-ghost elements
#ifdef SHELL_OMP
#pragma omp for
#endif
      for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
      {
        if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

        const ElemCache& cache = *elem_cache[eid];
        const Real* JxW = cache.JxW;
        const Real* phi = cache.phi;
        const std::vector<unsigned int>& dof_inds = cache.dof_indices;
        const unsigned int n_var_dofs = cache.dof_indices_ux.size();
        
        qe.resize(dof_inds.size());
        qx.reposition(UX_VAR * n_var_dofs, n_var_dofs);
        qy.reposition(UY_VAR * n_var_dofs, n_var_dofs);
        qz.reposition(UZ_VAR * n_var_dofs, n_var_dofs);

        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
          const SurfaceMetric& a_def = cache.a_def[qp];
          Point pressure_vec;

          if (pressure != 0)
          {
            switch (pressure_base)
            {
              case LOCAL:
              {
                // pressure in surface normal direction
                pressure_vec = a_def.vec[2];
                break;
              }
              case CARTESIAN:
              {
                pressure_vec(pressure_dir) = 1;
                break;
              }
              case CYLINDRICAL:
              {
                // pressure either in axial or radial direction
                if (pressure_dir == 0)
                  pressure_vec(0) = 1;
                else
                {
                  const Real inv_r = 1 / std::sqrt(a_def.xyz(1) * a_def.xyz(1) + a_def.xyz(2) * a_def.xyz(2));
                  if (std::isfinite(inv_r))
                  {
                    pressure_vec(1) = inv_r * a_def.xyz(1);
                    pressure_vec(2) = inv_r * a_def.xyz(2);
                  }
                }
                break;
              }
              case SPHERICAL:
              {
                // pressure in radial direction
                pressure_vec = a_def.xyz;
                const Real radius = pressure_vec.size();
                if (radius > 0)
                  pressure_vec *= 1 / radius;
                else
                  pressure_vec.zero();
                break;
              }
              default:
                ASSERT(pressure_base != INVALID_CS)
            }

            pressure_vec *= pressure;
          }

          if (volume_pressure != 0)
		  {
            pressure_vec.add_scaled(a_def.vec[2], volume_pressure);
//			std::cout << "P " << pressure_vec(0) << " " << pressure_vec(1) << " " << pressure_vec(2) << std::endl;
//		    HERE
		  }

          for (unsigned int i = 0; i < n_var_dofs; ++i)
          {
            const Real prefactor = JxW[qp] * phi[i * n_qp + qp];
            qx(i) += prefactor * pressure_vec(0);
            qy(i) += prefactor * pressure_vec(1);
            qz(i) += prefactor * pressure_vec(2);
          }
		}

#ifdef SHELL_OMP
#pragma omp critical(shell_assemble_vec)
#endif
        f.add_vector(qe, dof_inds);
      }
    }
  }

  if (gravity != 0)
  {
    // iterate over all nodes, including ghosts
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      const unsigned int dof = mesh.node(nid).dof_number(0, gravity_dir, 0);
      f[dof] -= gravity / Minv[dof];
    }
  }

  // interaction with cavity: forces that keep the shell inside the cavity
  if (shell_cavity_contact)
  {
    unsigned int Nc_cav_sum = 0;
    Real Efric_sum = 0, Ecav_sum = 0;
    Real cavity_force_sum = 0;
  
#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:Nc_cav_sum,Efric_sum,Ecav_sum,cavity_force_sum)
#endif
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      NodeCache& cache = node_cache[nid];
      if (cache.is_ghost) continue;
      
      const Node& node = mesh.node(nid);

      // check if the node is close to the cavity
      const Real one_plus_delta = cavity_isoparameter(nid);
      if (one_plus_delta > one_plus_delta_friction)
      {
        const Real depth = avg_R_eff * (std::sqrt(one_plus_delta) - 1); // penetration depth
        const Real area_node = avg_area_node(nid);
        const Real stiffness = stiff_density * area_node;
        Real f_normal = 0;
        if (depth > 0)
        {
          ++Nc_cav_sum;
          const Real Ecav_node = contact_energy(stiffness, depth);
          Ecav_sum += Ecav_node;
          ecav_vec[nid] += Ecav_node / area_node;
          f_normal = contact_force(stiffness, depth);
          cavity_force_sum += f_normal;
        }
        Point n;

        switch (cavity)
        {
          case BOX:
          {
            const Real x_red = cache.p(0) * cache.p(0) * inv_R2_eff[0];
            const Real y_red = cache.p(1) * cache.p(1) * inv_R2_eff[1];
            const Real z_red = cache.p(2) * cache.p(2) * inv_R2_eff[2];
            if (x_red > y_red && x_red > z_red)
              n(0) = sgn(cache.p(0));
            else if (y_red > z_red)
              n(1) = sgn(cache.p(1));
            else
              n(2) = sgn(cache.p(2));
            break;
          }
          case CYLINDER:
          {
            // treat radial and longitudinal directions separately for simplicity
            const Real x_red = cache.p(0) * cache.p(0) * inv_R2_eff[0];
            const Real y_red = cache.p(1) * cache.p(1) * inv_R2_eff[1];
            const Real z_red = cache.p(2) * cache.p(2) * inv_R2_eff[2];
            if (x_red > y_red + z_red)
              n(0) = sgn(cache.p(0));
            else
            {
              n(1) = cache.p(1) * inv_R2_eff[1];
              n(2) = cache.p(2) * inv_R2_eff[2];
            }
            break;
          }
          case ELLIPSOID:
          {
            n(0) = cache.p(0) * inv_R2_eff[0];
            n(1) = cache.p(1) * inv_R2_eff[1];
            n(2) = cache.p(2) * inv_R2_eff[2];
            break;
          }
          default:
            ASSERT(cavity != INVALID_CAV)
        }
        libmesh_assert(n.size() > 0);
        n *= 1 / n.size();

        // compute the normal velocity
        const unsigned int ux_dof = node.dof_number(0, UX_VAR, 0);
        const unsigned int uy_dof = node.dof_number(0, UY_VAR, 0);
        const unsigned int uz_dof = node.dof_number(0, UZ_VAR, 0);
        const Point vel(v[ux_dof], v[uy_dof], v[uz_dof]);
        const Real abs_v_normal = vel * n;

        Point f_contact; // zeroed by default

        if (cavity_friction)
        {
          Real Efric_node = 0;
          
          // compute friction force and energy, update static friction spring if necessary
          bool stick = (cache.spring_origin.size_sq() != 0);
          Friction::apply(Efric_node, f_contact, f_normal, const_cast<const Point&>(n), vel, abs_v_normal, friction_stop_velocity, stiff_static_friction,
                          mu_s_cavity, mu_d_cavity, cavity_friction_damping, stick, cache.spring_origin, cache.p);
          
          // measure static friction energy (total, and nodal density)
          Efric_sum += Efric_node;
          efric_vec[nid] += Efric_node / area_node;
        }
        
        // add the normal damping force
        f_normal += cavity_damping * stiffness * abs_v_normal;

        f_contact -= f_normal * n;
        f[ux_dof] += f_contact(0);
        f[uy_dof] += f_contact(1);
        f[uz_dof] += f_contact(2);
      }
      else
      {
        // discard the static friction spring only if the contact distance is larger
        // than a certain tolerance value, because D fluctuates around zero when a contact happens,
        // so we need to sustain the spring also when the two elements are slightly apart
        cache.spring_origin.zero();
      }
    }
    
    Nc_cav += Nc_cav_sum;
    Efric += Efric_sum;
    Ecav += Ecav_sum;
    cavity_force += cavity_force_sum;
  }

  // prevent the accumulation of angular momentum
  // (e.g. induced by boundary conditions)
  if (damping_rotation > 0)
  {
    RealVectorValue L; // angular momentum about the origin
    RealTensorValue I; // moment of inertia tensor with respect to the origin
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      const Node& node = mesh.node(nid);

      const unsigned int ux_dof = node.dof_number(0, UX_VAR, 0);
      const unsigned int uy_dof = node.dof_number(0, UY_VAR, 0);
      const unsigned int uz_dof = node.dof_number(0, UZ_VAR, 0);

      const Point& ri = node_cache[nid].p;
      const Point vi(v[ux_dof], v[uy_dof], v[uz_dof]);
      const Real mi = 1 / Minv[ux_dof];
      const Point rxv = ri.cross(vi);
#ifdef SHELL_OMP
#pragma omp critical(shell_angular_momentum)
#endif
      {
        L += mi * rxv;

        I(0,0) += mi * (ri(1) * ri(1) + ri(2) * ri(2));
        I(1,1) += mi * (ri(0) * ri(0) + ri(2) * ri(2));
        I(2,2) += mi * (ri(0) * ri(0) + ri(1) * ri(1));

        I(0,1) -= mi * ri(0) * ri(1);
        I(0,2) -= mi * ri(0) * ri(2);
        I(1,2) -= mi * ri(1) * ri(2);
      }
    }

    // fill in the lower half
    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);

    // angular velocity
    const RealVectorValue omega = inverse(I) * L;
    RealTensorValue W;
    cross_product_matrix(W, omega);

#ifdef SHELL_OMP
#pragma omp parallel for
#endif
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      const Node& node = mesh.node(nid);
      const Point vi = W * node_cache[nid].p;
      f[node.dof_number(0, UX_VAR, 0)] -= damping_rotation * vi(0);
      f[node.dof_number(0, UY_VAR, 0)] -= damping_rotation * vi(1);
      f[node.dof_number(0, UZ_VAR, 0)] -= damping_rotation * vi(2);
    }
  }

  if (shell_shell_contact != NONE)
    add_shell_shell_contact_force();

  PERFLOG_STOP("compute_external_force()");
}

// adds the external forces resulting from self-contact to the vector
void ShellSystem::add_shell_shell_contact_force()
{
  PERFLOG_START("add_shell_shell_contact_force()");
  
  const MeshBase& mesh = get_mesh();

  unsigned int Nc_sum = 0;
  Real Efric_sum = 0, Eself_sum = 0;
#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:Nc_sum,Efric_sum,Eself_sum)
#endif
  for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
  {
    ShellShellContact* contact = shell_shell_contact_candidates[i];
    if (!contact) continue;

    // check if the elements are close
    contact->compute();
    if (contact->d2 < d2_friction)
    {
      const Real& area1 = elem_cache[contact->eid[0]]->area;
      const Real& area2 = elem_cache[contact->eid[1]]->area;
      const Real area_min = std::min(area1, area2);
      libmesh_assert(area_min > 0);

      const Real abs_dc = std::sqrt(contact->d2);
      const Real depth = thickness - abs_dc; // penetration depth
      const Point n = contact->dc / abs_dc; // normal vector
      const Real stiffness = stiff_density * area_min;
      Real f_normal = 0;
      Real Eself_contact = 0, Efric_contact = 0;
      if (depth > 0)
      {
        ++Nc_sum;
        f_normal = shell_shell_contact_force(stiffness, abs_dc);
        Eself_contact = shell_shell_contact_energy(stiffness, abs_dc);
      }
      Eself_sum += Eself_contact;

      // compute the relative normal velocity
      Point v10, v11, v12, v20, v21, v22;
      get_nodal_props(v, contact->eid[0], v10, v11, v12);
      get_nodal_props(v, contact->eid[1], v20, v21, v22);
      const Real u[2] = { 1 - contact->v[0] - contact->w[0], 1 - contact->v[1] - contact->w[1] };
      const Point v_rel = u[1] * v20 + contact->v[1] * v21 + contact->w[1] * v22
                       - (u[0] * v10 + contact->v[0] * v11 + contact->w[0] * v12);
      const Real abs_v_normal = v_rel * n;

      Point f_contact; // zeroed by default

      if (self_friction)
      {
        // compute friction force and energy, update static friction spring if necessary
        Friction::apply(Efric_contact, f_contact, f_normal, n, v_rel, abs_v_normal, friction_stop_velocity, stiff_static_friction,
                        mu_s, mu_d, self_friction_damping, contact->stick, contact->spring_origin, contact->point());
        
        // sum up the static friction energy
        Efric_sum += Efric_contact;
      }
      
      // add the normal damping force
      f_normal -= self_damping * stiffness * abs_v_normal;

      f_contact += f_normal * n;

      // distribute the energies and forces among the six nodes, using linear interpolation
#ifdef SHELL_OMP
#pragma omp critical(shell_assemble_vec)
#endif
      {
        for (unsigned int e = 0; e < 2; ++e)
        {
          const Elem* elem = mesh.elem(contact->eid[e]);
          const Real& v = contact->v[e];
          const Real& w = contact->w[e];
          const Real u = 1 - v - w;
          const Real prefactor0 = 0.5 * u / area_min;
          const Real prefactor1 = 0.5 * v / area_min;
          const Real prefactor2 = 0.5 * w / area_min;

          eself_vec[elem->node(0)] += prefactor0 * Eself_contact;
          eself_vec[elem->node(1)] += prefactor1 * Eself_contact;
          eself_vec[elem->node(2)] += prefactor2 * Eself_contact;
          efric_vec[elem->node(0)] += prefactor0 * Efric_contact;
          efric_vec[elem->node(1)] += prefactor1 * Efric_contact;
          efric_vec[elem->node(2)] += prefactor2 * Efric_contact;
       
          const Node* n0 = elem->get_node(0);
          const Node* n1 = elem->get_node(1);
          const Node* n2 = elem->get_node(2);
          const Real prefactor = (e == 0 ? -0.5 : 0.5);
          for (unsigned int var = 0; var < 3; ++var)
          {
            const Real fc = prefactor * f_contact(var);
            f[n0->dof_number(0, var, 0)] += u * fc;
            f[n1->dof_number(0, var, 0)] += v * fc;
            f[n2->dof_number(0, var, 0)] += w * fc;
          }
        }
      }
    }
    else
    {
      // discard the static friction spring only if the contact distance is larger
      // than a certain tolerance value, because D fluctuates around zero when a contact happens,
      // so we need to sustain the spring also when the two elements are slightly apart
      contact->spring_origin.zero();
      contact->stick = false;
    }
  }
  Nc += Nc_sum;
  Efric += Efric_sum;
  Eself += Eself_sum;

  PERFLOG_STOP("add_shell_shell_contact_force()");
}

// computes the self-contact pairs
void ShellSystem::compute_shell_shell_contact_candidates(bool force)
{
  if (shell_shell_contact == NONE) return;

  PERFLOG_START("compute_shell_shell_contact_candidates()");

  // rebuild the linked cell list only when necessary, unless we enforce it
  Real cell_size = cell_size_mult * avg_elem_size;
  Real tolerance = cell_margin * cell_size;
  if (force || du_sum >= tolerance)
  {
    // reuse the list of contact pairs computed below for a while
    du_sum = 0;

    // recompute the close neighborhood unless the shell is growing uniformly in all 3 directions
    const Real growth3 = thickness / thickness0;
    if (growth1 != growth3 || growth2 != growth3)
      build_topological_neighborhood();

    // compute the element bounding boxes
    compute_bounding_box(bounds, avg_elem_size);

    // determine size
    cell_size = cell_size_mult * avg_elem_size;
    tolerance = cell_margin * cell_size;
    cell_size += tolerance;

    // the maximum distance between two elements for considering contacts
    Real d_max = thickness + tolerance;
    
    // for smooth contact evaluation, the shell surfaces could be an unknown amount
    // closer than the distance between the faceted elements,
    // so we add another average shell element size as a safety margin
    if (shell_shell_contact == SMOOTH)
      d_max += avg_elem_size;
    
    const Real d2_max = d_max * d_max;

    // allocate memory
    unsigned int n_cells_dir[3];
    n_cells_dir[0] = static_cast<int>((bounds.size(0) + d_max) / cell_size + 1);
    n_cells_dir[1] = static_cast<int>((bounds.size(1) + d_max) / cell_size + 1);
    n_cells_dir[2] = static_cast<int>((bounds.size(2) + d_max) / cell_size + 1);
    const unsigned int n_cells = n_cells_dir[0] * n_cells_dir[1] * n_cells_dir[2];

    // when the system explodes it may happen that n_cells is extremely large, and
    // allocating enough memory for the cells fails
    if (n_cells > n_cells_max)
    {
      register_problem(TOO_MANY_CELLS);
      problem_info.set<unsigned int>("n_cells") = n_cells;
      problem_info.set<unsigned int>("n_cells_max") = n_cells_max;
      problem_info.set<Real>("cell_size") = cell_size;
      return;
    }

    const MeshBase& mesh = get_mesh();
    std::vector<std::vector<unsigned int> > elem_to_cell_map(mesh.n_elem());
    std::vector<std::vector<unsigned int> > cell_to_elem_map(n_cells);

    // Put the non-ghost elements into all cells their bounding box
    // reaches into. This loop is tough to parallelize.
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

      BoundingBox<Real>& elem_bounds = elem_cache[eid]->bounds;

      // find the cell indices in each direction
      int cell_ind_min[3], cell_ind_max[3];
      for (unsigned int var = 0; var < 3; ++var)
      {
        cell_ind_min[var] = static_cast<int>((elem_bounds.min[var] - bounds.min[var]        ) / cell_size);
        cell_ind_max[var] = static_cast<int>((elem_bounds.max[var] - bounds.min[var] + d_max) / cell_size);

        // assert that the indices are within the valid range
        libmesh_assert(cell_ind_min[var] >= 0);
        libmesh_assert(cell_ind_max[var] >= 0);
        libmesh_assert(cell_ind_min[var] < (int)n_cells_dir[var]);
        libmesh_assert(cell_ind_max[var] < (int)n_cells_dir[var]);
      }

      // add the element to all those cells
      for (int cx = cell_ind_min[0]; cx <= cell_ind_max[0]; ++cx)
      {
        for (int cy = cell_ind_min[1]; cy <= cell_ind_max[1]; ++cy)
        {
          for (int cz = cell_ind_min[2]; cz <= cell_ind_max[2]; ++cz)
          {
            // the global cell index
            const int c = cx * n_cells_dir[1] * n_cells_dir[2]
                        + cy * n_cells_dir[2]
                        + cz;
            libmesh_assert(c < (int)n_cells);

            elem_to_cell_map[eid].push_back(c);
            cell_to_elem_map[c].push_back(eid);
          }
        }
      }
    }
    
    // Drop all contact candidates whose elements are too far apart
    // and which are not connected by a static friction spring, and at the same
    // time determine the first empty slot in the list of contact candidates.
    // This loop is tough to parallelize.
    unsigned int first_empty = shell_shell_contact_candidates.size();
    for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
    {
      ShellShellContact* contact = shell_shell_contact_candidates[i];
      if (contact && !contact->stick && contact->d2 > d2_max)
      {
        elem_cache[contact->eid[0]]->in_contact_with.remove(contact->eid[1]);
        delete contact;
        shell_shell_contact_candidates[i] = NULL; // mark this slot as emtpy
        first_empty = std::min(first_empty, i);
      }
    }

    // iterate over element pairs
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
    for (unsigned int eid1 = 0; eid1 < mesh.n_elem(); ++eid1)
    {
      if (static_cast<const Tri3SD*>(mesh.elem(eid1))->is_ghost()) continue;

      // iterate over the cell indices in which this element lies
      const std::vector<unsigned int>& cells = elem_to_cell_map[eid1];
      for (unsigned int cidx = 0; cidx < cells.size(); ++cidx)
      {
        // iterate over the elements contained in this cell
        std::vector<unsigned int>& elems_in_cell = cell_to_elem_map[cells[cidx]];
        for (unsigned int i = 0; i < elems_in_cell.size(); ++i)
        {
          const unsigned int eid2 = elems_in_cell[i];
          
          // prevent adding element pairs twice
          if (eid1 >= eid2) continue;

          // exclude topologically close elements
          if (elems_are_close(eid1, eid2)) continue;

          // exclude elements pairs that are in the list of contact candidates already
          const std::list<unsigned int>& contacts_with_low_id = elem_cache[eid1]->in_contact_with;
          if (std::find(contacts_with_low_id.begin(), contacts_with_low_id.end(), eid2) != contacts_with_low_id.end())
            continue;

          // exclude element pairs whose bounding boxes aren't close enough
          if (elem_cache[eid1]->bounds.apart(elem_cache[eid2]->bounds, d_max))
            continue;

          // create a new contact candidate
          ShellShellContact* contact_candidate;
          switch (shell_shell_contact)
          {
            case LINEAR:
            {
              contact_candidate = new ShellShellContact(eid1, eid2);
              break;
            }
            case SMOOTH:
            {
              contact_candidate = new SmoothShellShellContact(eid1, eid2);
              break;
            }
            default:
              ASSERT(shell_shell_contact != INVALID_CT)
          }
          
          // check if the elements are close enough
          contact_candidate->compute();
          if (contact_candidate->d2 > d2_max)
          {
            delete contact_candidate;
            continue;
          }

          // finally, if we made it through here, we have found an element pair
          // potentially making contact, thus add it to the list of contact candidates
#ifdef SHELL_OMP
#pragma omp critical(shell_cc_addition)
#endif
          {
            if (first_empty == shell_shell_contact_candidates.size())
            {
              ++first_empty;
              shell_shell_contact_candidates.push_back(contact_candidate);
            }
            else
            {
              libmesh_assert(shell_shell_contact_candidates[first_empty] == NULL);
              shell_shell_contact_candidates[first_empty] = contact_candidate;
              
              // find the next empty slot
              std::vector<ShellShellContact*>::iterator next_empty = std::find(shell_shell_contact_candidates.begin() + first_empty + 1, shell_shell_contact_candidates.end(), (ShellShellContact*)NULL);
              first_empty = next_empty - shell_shell_contact_candidates.begin();
            }
          }
          
          elem_cache[eid1]->in_contact_with.push_back(eid2);
        }
      }
    }
  }
/*
  if (shell_shell_contact_candidates.size() > 0)
  {
    *os << "Shell-shell contact candidates in timestep " << time_step << ":\n";
    for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
    {
      ShellShellContact* contact = shell_shell_contact_candidates[i];
      if (!contact) continue;
      
      contact->compute();
      const Elem* e1 = get_mesh().elem(contact->eid[0]);
      const Elem* e2 = get_mesh().elem(contact->eid[1]);
      *os << "elems = " << contact->eid[0] << "," << contact->eid[1] << ",\tnodes = "
          << e1->node(0) << "," << e1->node(1) << "," << e1->node(2) << ","
          << e2->node(0) << "," << e2->node(1) << "," << e2->node(2) << ",\t"
          << "dist = " << std::sqrt(contact->d2) << std::endl;
    }
  }
*/
  PERFLOG_STOP("compute_shell_shell_contact_candidates()");
}

// enforces the boundary conditions using ghosted nodes
void ShellSystem::apply_boundary_conditions()
{
  if (bc_every <= 0 || (time_step-1) % bc_every > 0) return;
  
  PERFLOG_START("apply_boundary_conditions()");
  
  const Real psh = 2e2;
  const MeshBase& mesh = get_mesh();
  	
/*
  const Real infty = std::numeric_limits<Real>::infinity();
  std::vector<RealVectorValue> points, targets, forces;
  
  points.push_back(RealVectorValue(10,0,0));
  points.push_back(RealVectorValue(0,10,0));
  points.push_back(RealVectorValue(-10,0,0));
  points.push_back(RealVectorValue(0,-10,0));

  const Real angle = pi/2 - pi/10;
  points.push_back(RealVectorValue(10*std::cos(angle),0,10*std::sin(angle)));
  points.push_back(RealVectorValue(0,10*std::cos(angle),10*std::sin(angle)));
  points.push_back(RealVectorValue(-10*std::cos(angle),0,10*std::sin(angle)));
  points.push_back(RealVectorValue(0,-10*std::cos(angle),10*std::sin(angle)));

  targets.push_back(RealVectorValue(infty,0,0));
  targets.push_back(RealVectorValue(0,infty,infty));
  targets.push_back(RealVectorValue(infty,0,0));
  targets.push_back(RealVectorValue(0,infty,infty));

  targets.push_back(RealVectorValue(infty,0,infty));
  targets.push_back(RealVectorValue(0,infty,infty));
  targets.push_back(RealVectorValue(infty,0,infty));
  targets.push_back(RealVectorValue(0,infty,infty));

  forces.push_back(RealVectorValue(-fulltime*psh,0,0));
  forces.push_back(RealVectorValue(0,fulltime*psh,0));
  forces.push_back(RealVectorValue(fulltime*psh,0,0));
  forces.push_back(RealVectorValue(0,-fulltime*psh,0));
  
 // if (time_step%1000 == 0) {std::cout << "Increment " << increment_factor << " " << fulltime << std::endl;}	

  forces.push_back(RealVectorValue(0,0,0));
  forces.push_back(RealVectorValue(0,0,0));
  forces.push_back(RealVectorValue(0,0,0));
  forces.push_back(RealVectorValue(0,0,0));

  for (unsigned int i = 0; i < points.size(); ++i)
  {
    unsigned nearest_nid = Node::invalid_id;
    Real nearest_d2 = std::numeric_limits<Real>::infinity();
    for (unsigned int nid = 0; nid < get_mesh().n_nodes(); ++nid)
    {
      const Real d2 = (points[i] - get_mesh().point(nid)).size_sq();
      if (d2 < nearest_d2)
      {
        nearest_d2 = d2;
        nearest_nid = nid;
      }
    }
    const Node *n = get_mesh().node_ptr(nearest_nid);
	if (i < 4) press[i] = n->id();
	
	if (n->processor_id() == processor_id())
    {
      for (unsigned int var = 0; var < 3; ++var)
      {
        const unsigned int dof = n->dof_number(number(), var, 0);
        const Real target = targets[i](var);
        if (std::isfinite(target))
        {
          u_pred[dof] = target;
          v_pred[dof] = 0;
		  a_pred[dof] = 0;
		}
        f[dof] -= -forces[i](var);  
	  }
	}
  }
*/
  // iterate over all boundary edges
  for (unsigned int i = 0; i < boundary_cache.size(); ++i)
  {
    const BoundaryCache& cache = boundary_cache[i];
    
    for (unsigned int var = 0; var < 3; ++var)
    {
      // only apply the boundary conditions in the specified directions
      if (!(bc_coords & (1 << var)))
        continue;
      
      const unsigned int dof [4] = { mesh.node(cache.nids[0]).dof_number(0, var, 0),
                                     mesh.node(cache.nids[1]).dof_number(0, var, 0),
                                     mesh.node(cache.nids[2]).dof_number(0, var, 0),
                                     mesh.node(cache.nids[3]).dof_number(0, var, 0) };

      switch (shell_boundary)
      {
        case CLAMPED: // Cirak et al.: u1 = u2 = u3 = u4 = 0
        {
          for (unsigned int n = 0; n < 4; ++n)
          {
            u_pred[dof[n]] = 0;
            v_pred[dof[n]] = 0;
            a_pred[dof[n]] = 0;
          }
          break;
        }
        case CLAMPED_WEAK: // penalize to u1 = u2 = u3 = u4 = 0
        {
          for (unsigned int n = 0; n < 4; ++n)
            f[dof[n]] -= bc_stiff * u_pred[dof[n]];
          break;
        }
        case PINNED: // u2 = u3 = 0; u4 = -u1
        {
          u_pred[dof[1]] = 0;
          v_pred[dof[1]] = 0;
          a_pred[dof[1]] = 0;
          
          u_pred[dof[2]] = 0;
          v_pred[dof[2]] = 0;
          a_pred[dof[2]] = 0;

          u_pred[dof[3]] = -u_pred[dof[0]];
          v_pred[dof[3]] = -v_pred[dof[0]];
          a_pred[dof[3]] = -a_pred[dof[0]];
          break;
        }
        case PINNED_WEAK: // penalize to u2 = u3 = 0; u4 = -u1
        {
          f[dof[1]] -= bc_stiff * u_pred[dof[1]];
          f[dof[2]] -= bc_stiff * u_pred[dof[2]];
          f[dof[3]] -= bc_stiff * (u_pred[dof[3]] + u_pred[dof[0]]);
          break;
        }
        case FREE: // u4 = u2 + u3 - u1
        {
          u_pred[dof[3]] = u_pred[dof[1]] + u_pred[dof[2]] - u_pred[dof[0]];
          v_pred[dof[3]] = v_pred[dof[1]] + v_pred[dof[2]] - v_pred[dof[0]];
          a_pred[dof[3]] = a_pred[dof[1]] + a_pred[dof[2]] - a_pred[dof[0]];
          break;
        }
        case FREE_WEAK: // penalize to u4 = u2 + u3 - u1
        {
          f[dof[3]] -= bc_stiff * (u_pred[dof[3]] + u_pred[dof[0]] - u_pred[dof[1]] - u_pred[dof[2]]);
          break;
        }
        default:
          break;
      }
    }
  }

  PERFLOG_STOP("apply_boundary_conditions()");
}

// computes the growth deformation gradient tensors G at the quadrature points
void ShellSystem::apply_growth()
{
  PERFLOG_START("apply_growth()");

  const MeshBase& mesh = get_mesh();
  
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    const Elem* elem = mesh.elem(eid);
    if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;

    const ElemCache& cache = *elem_cache[eid];

    for (unsigned int qp = 0; qp < n_qp; ++qp)
    {
      const SurfaceMetric& a_ref = cache.a_ref[qp];
      RealTensorValue& GV        = cache.GV[qp];
      RealTensorValue Ginv;

      if (grow_dynamically)
      {
        RealVectorValue& alpha = cache.alpha[qp];
        RealVectorValue& beta  = cache.beta[qp];
        RealVectorValue& G     = cache.G[qp];
        RealVectorValue& G_old = cache.G_old[qp];

        // increment the growth tensor by a growth rate proportional to the strains
        G_old = G;
        G.add_scaled(alpha, dt * metric_growth_rate);
        G.add_scaled(beta , dt * curvature_growth_rate);

        // cheap version of Ginv = inverse(G)
        const Real det = G(0) * G(1) - 0.25 * G(2) * G(2); // G contains the Voigt factor 2
        libmesh_assert(det != 0);
        Ginv(0,0) = G(1) / det;
        Ginv(1,1) = G(0) / det;
        Ginv(0,1) = Ginv(1,0) = -0.5 * G(2) / det;
        Ginv(2,2) = 1;
      }
      else
      {
        // compute the matrix B containing the base vectors w.r.t. which the growth tensor is diagonal
        RealTensorValue B;
        switch (growth_base)
        {
          case LOCAL:
          {
            // use the normalized reference surface tangent vectors as base
            for (unsigned int i = 0; i < 3; ++i)
            {
              libmesh_assert(a_ref.vec[i].size() > 0);
              const Real inv_length = 1 / a_ref.vec[i].size();
              B(0,i) = inv_length * a_ref.vec[i](0);
              B(1,i) = inv_length * a_ref.vec[i](1);
              B(2,i) = inv_length * a_ref.vec[i](2);
            }
            break;
          }
          case CARTESIAN:
          {
            B = Id;
            break;
          }
          case CYLINDRICAL:
          {
            // note: the cylinder axis is x, i.e. the (r,phi) plane is the (y,z) plane
            const Real inv_r = 1 / std::sqrt(a_ref.xyz(1) * a_ref.xyz(1) + a_ref.xyz(2) * a_ref.xyz(2));
            if (std::isfinite(inv_r))
            {
              B(0,0) =  1;
              B(1,1) =  inv_r * a_ref.xyz(1); //  cos(a)
              B(2,1) =  inv_r * a_ref.xyz(2); //  sin(a)
              B(1,2) = -inv_r * a_ref.xyz(2); // -sin(a)
              B(2,2) =  inv_r * a_ref.xyz(1); //  cos(a)
            }
            else
            {
              GV = Id;
              continue;
            }
            break;
          }
          default:
            ASSERT(growth_base != INVALID_CS)
        }

        // compute the matrix A containing the unnormalized reference surface base vectors a1, a2, a3
        RealTensorValue A;
        for (unsigned int i = 0; i < 3; ++i)
        {
          A(0,i) = a_ref.vec[i](0);
          A(1,i) = a_ref.vec[i](1);
          A(2,i) = a_ref.vec[i](2);
        }
        A(0,2) *= a_ref.jac;
        A(1,2) *= a_ref.jac;
        A(2,2) *= a_ref.jac;
        
        RealVectorValue& G = cache.G[qp];
#ifdef SHELL_OMP
#pragma omp critical(shell_use_xyz)
#endif
        {
          // re-evaluate the in-plane growth factor functions at the Gauss point (x,y,z)
          xyz = a_ref.xyz;
          functions["growth1"]->evaluate();
          functions["growth2"]->evaluate();

          // the inverse growth tensor in diagonal form (in basis B)
          libmesh_assert(growth1 > 0);
          libmesh_assert(growth2 > 0);
          G(0) = growth1;
          G(1) = growth2;
          Ginv(growth_dir[0],growth_dir[0]) = 1 / growth1;
          Ginv(growth_dir[1],growth_dir[1]) = 1 / growth2;
        }
        G(2) = 0;
        Ginv(growth_dir[2],growth_dir[2]) = 1;

        // change of basis B -> A: G^-1' = A^-1 B G^-1 B^-1 A
        A = inverse(B) * A;
        Ginv = inverse(A) * Ginv * A;
      }

      // G is the matrix in Voigt notation for the transformation G^-T X G^-1 for some matrix X.
      // This can be written as GV*XV, where XV is X as a matrix operator in Voigt notation and GV as below.
      // This is allowed becase the above transformation still results in a symmetric matrix.
      GV(0,0) =     Ginv(0,0) * Ginv(0,0);
      GV(0,1) =     Ginv(1,0) * Ginv(1,0);
      GV(0,2) = 2 * Ginv(0,0) * Ginv(1,0);

      GV(1,0) =     Ginv(0,1) * Ginv(0,1);
      GV(1,1) =     Ginv(1,1) * Ginv(1,1);
      GV(1,2) = 2 * Ginv(0,1) * Ginv(1,1);

      GV(2,0) = Ginv(0,0) * Ginv(0,1);
      GV(2,1) = Ginv(1,0) * Ginv(1,1);
      GV(2,2) = Ginv(0,1) * Ginv(1,0) + Ginv(0,0) * Ginv(1,1);
    }
  }

  PERFLOG_STOP("apply_growth()");
}

// reverses the proposed groth tensors (necessary for adaptive time stepping, where a step may be rejected)
void ShellSystem::reject_growth()
{
  const MeshBase& mesh = get_mesh();
  
#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    const Elem* elem = mesh.elem(eid);
    if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;

    ElemCache& cache = *elem_cache[eid];
    for (unsigned int qp = 0; qp < n_qp; ++qp)
      cache.G[qp] = cache.G_old[qp];
  }
}

// determine the average element area around a node
Real ShellSystem::avg_area_node(unsigned int nid)
{
  if (node_cache[nid].is_ghost) return 0;

  Real node_area = 0;
  unsigned int n_real_elems = 0;
  for (unsigned int i = 0; i < nodes_to_elem_map[nid].size(); ++i)
  {
    const Elem* elem = nodes_to_elem_map[nid][i];
    if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;
    
    ++n_real_elems;
    node_area += elem_cache[elem->id()]->area;
  }
  libmesh_assert(n_real_elems > 0);
  node_area /= n_real_elems;
  libmesh_assert(node_area > 0);

  return node_area;
}

// computes the total nodal positions (mesh + displacement)
void ShellSystem::compute_position()
{
  PERFLOG_START("compute_position()");

  const MeshBase& mesh = get_mesh();

#ifdef SHELL_OMP
#pragma omp parallel for
#endif
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    const Node& node = mesh.node(nid);
    node_cache[nid].p = node;
    for (unsigned int var = 0; var < 3; ++var)
      node_cache[nid].p(var) += u_pred[node.dof_number(0, var, 0)];
  }

  PERFLOG_STOP("compute_position()");
}

// computes the surface area and the volume currently enclosed by the shell
void ShellSystem::compute_enclosed_volume()
{
  // the enclosed volume can only be calculated for closed meshes
  if (n_ghost_elem > 0)
    return;

  PERFLOG_START("compute_enclosed_volume()");
  
  const MeshBase& mesh = get_mesh();
  
  area = volume = 0;
  Real area_sum = 0, volume_sum = 0;
#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:area_sum,volume_sum)
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    const Elem * elem = mesh.elem(eid);
    RealTensorValue triad;
    for (unsigned int n = 0; n < 3; ++n)
    {
      const unsigned int nid = elem->node(n);
      for (unsigned int var = 0; var < 3; ++var)
        triad(var,n) = node_cache[nid].p(var);
    }
    area_sum += elem_cache[eid]->area;
    volume_sum += triad.det();
  }
  area = area_sum;
  volume = volume_sum;

//  std::cout << "Volume_sum " << volume << std::endl;
  volume = std::fabs(volume/6) - area * 0.5 * thickness;
  libmesh_assert(volume > 0);

  PERFLOG_STOP("compute_enclosed_volume()");
}

// computes the bounding boxes for all elements and the shell in total
void ShellSystem::compute_bounding_box(BoundingBox<Real>& bb, Real& avg_elem_size)
{
  PERFLOG_START("compute_bounding_box()");

  bb.reset();
  Real avg_elem_size_sum = 0;

  const MeshBase& mesh = get_mesh();
#ifdef SHELL_OMP
#pragma omp parallel
#endif
  {
    BoundingBox<Real> bb_local;
    unsigned int processed = 0;
#ifdef SHELL_OMP
#pragma omp for nowait reduction(+:avg_elem_size_sum)
#endif
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      const Elem * elem = mesh.elem(eid);
      if (static_cast<const Tri3SD*>(elem)->is_ghost()) continue;
      
      ++processed;
      BoundingBox<Real>& elem_bounds = elem_cache[eid]->bounds;
      elem_bounds.reset();
      for (unsigned int n = 0; n < elem->n_nodes(); ++n)
        elem_bounds.include(node_cache[elem->node(n)].p);
      bb_local.include(elem_bounds);
      avg_elem_size_sum += elem_bounds.max_size();
    }
    
    if (processed > 0)
    {
#ifdef SHELL_OMP
#pragma omp critical(shell_include_bb)
#endif
      bb.include(bb_local);
    }
  }

  avg_elem_size = avg_elem_size_sum / mesh.n_elem();

  PERFLOG_STOP("compute_bounding_box()");
}

// checks if the max. displacement increment is small enough to guarantee a physical solution (without tunnelling)
// and accumulates the maximum displacement for knowing when the linked cell lists need to be rebuilt
void ShellSystem::check_increment(bool require_change)
{
  PERFLOG_START("check_increment()");

  const MeshBase& mesh = get_mesh();
  Real du_max = 0;

#ifdef SHELL_OMP
#pragma omp parallel
#endif
  {
    Real du_max_local = 0;
    
#ifdef SHELL_OMP
#pragma omp for nowait
#endif
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      if (node_cache[nid].is_ghost) continue;
      
      const Node& node = mesh.node(nid);
      const unsigned int ux_dof = node.dof_number(0, UX_VAR, 0);
      const unsigned int uy_dof = node.dof_number(0, UY_VAR, 0);
      const unsigned int uz_dof = node.dof_number(0, UZ_VAR, 0);

      const Real dx = u_pred[ux_dof] - u[ux_dof];
      const Real dy = u_pred[uy_dof] - u[uy_dof];
      const Real dz = u_pred[uz_dof] - u[uz_dof];

      const Real du = dx * dx + dy * dy + dz * dz;
      du_max_local = std::max(du_max_local, du);
    }
    
#ifdef SHELL_OMP
#pragma omp critical(shell_du_max)
#endif
    du_max = std::max(du_max, du_max_local);
  }

  du_max = std::sqrt(du_max);

  if (std::isfinite(du_max))
  {
    du_sum += 2 * du_max; // worst case: the two fastest nodes move towards each other, hence the factor 2

    if (require_change && du_max == 0)
      register_problem(NO_DISP_CHANGE);
    else if (shell_shell_contact != NONE && du_max > 0.5 * thickness)
    {
      register_problem(LARGE_DISP);
      problem_info.set<Real>("du_max") = du_max;
      problem_info.set<Real>("thickness") = thickness;
    }
  }
  else
  {
    register_problem(NOT_FINITE);
    problem_info.set<Real>("du_max") = du_max;
  }

  PERFLOG_STOP("check_increment()");
}

// adds kinetic energy to the system at a specified rate
void ShellSystem::perturb()
{
  if (perturbation_rate <= 0 || perturbation_energy <= 0)
    return;

  PERFLOG_START("perturb()");

  const MeshBase& mesh = get_mesh();
  const Real probability = perturbation_rate * dt;
  const unsigned int n_real_nodes = mesh.n_nodes() - n_ghost_nodes;
  const Real kT = 2 * perturbation_energy / (3 * n_real_nodes); // the "added temperature" per dof

  // this is not parallelized because drand48 and randn are not thread-safe
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    if (!node_cache[nid].is_ghost && drand48() < probability)
    {
      const Node& node = mesh.node(nid);
      const unsigned int ux_dof = node.dof_number(0, UX_VAR, 0);
      const Real sigma = std::sqrt(kT * Minv[ux_dof]);
      
      v_pred[ux_dof                       ] += randn<Real>(0, sigma);
      v_pred[node.dof_number(0, UY_VAR, 0)] += randn<Real>(0, sigma);
      v_pred[node.dof_number(0, UZ_VAR, 0)] += randn<Real>(0, sigma);
    }
  }

  PERFLOG_STOP("perturb()");
}

// resets some measured properties of the shell
void ShellSystem::reset_properties()
{
  const MeshBase& mesh = get_mesh();
  
  // reset all measurement vectors
  std::map<std::string, ScalarVector*>::iterator vit = elem_vectors.begin();
  std::map<std::string, ScalarVector*>::iterator end_vit = elem_vectors.end();
  for (; vit != end_vit; ++vit)
  {
    vit->second->resize(mesh.n_elem());
    vit->second->zero();
  }
  
  vit = node_vectors.begin();
  end_vit = node_vectors.end();
  for (; vit != end_vit; ++vit)
  {
    vit->second->resize(mesh.n_nodes());
    vit->second->zero();
  }
  
  // reset some measurement scalars
  Nc = 0, Nc_cav = 0; // number of contacts
  Eself = 0, Ecav = 0, Efric = 0, Etens = 0, Ebend = 0; // static energies
  cavity_force = 0;
}

// measures some properties of the shell
void ShellSystem::measure_properties()
{
  PERFLOG_START("measure_properties()");
  
  const MeshBase& mesh = get_mesh();
  
  // The properties measured from the element cache (area & curvatures) are computed on the predicted
  // configuration instead of the corrected configuration to avoid having to recompute the deformed
  // surface metrics for only measurement purposes. This induces a tiny measurement error.
  
  // measure curvature, see
  // http://en.wikipedia.org/wiki/Shape_operator
  // http://de.wikipedia.org/wiki/Mittlere_Krümmung#Eigenschaften
  // http://de.wikipedia.org/wiki/Gaußsche_Krümmung#Eigenschaften
  Real kGauss = 0, kMean = 0, area_sum = 0;
 
  Real xl = 0;
  Real yl = 0;
  Point p1, p2, p3, p4;

#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:kGauss,kMean,area_sum)
#endif
  for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
  {
    if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

    const ElemCache& cache = *elem_cache[eid];
    libmesh_assert(cache.area > 0);

    Real kGauss_elem = 0, kMean_elem = 0;
    for (unsigned int qp = 0; qp < n_qp; ++qp)
    {
      // growth-enhanced fundamental forms
      RealVectorValue ff1, ff2;
      if (grow_in_plane)
      {
        ff1 = cache.GV[qp] * cache.a_def[qp].ff1;
        ff2 = cache.GV[qp] * cache.a_def[qp].ff2;
      }
      else
      {
        ff1 = cache.a_def[qp].ff1;
        ff2 = cache.a_def[qp].ff2;
      }
      
      const Real w = cache.a_def[qp].jac * qw[qp];
      const Real invdet = 1 / (ff1(0) * ff1(1) - ff1(2) * ff1(2));
      kGauss_elem += w * invdet * (ff2(0) * ff2(1) - ff2(2) * ff2(2));
      kMean_elem  += w * invdet * 0.5 * (ff2(0) * ff1(1) + ff2(1) * ff1(0) - 2 * ff2(2) * ff1(2));
    }
 
 	if (time_step%1000 == 0)
	{
   		// To find the limit surface point for the edges of the elements
		if (eid == 0)
		{
			SubdivEvaluation::evaluate(this, eid, xl, yl, p1);
		} else if (eid == 512){
			SubdivEvaluation::evaluate(this, eid, xl, yl, p2);
		} else if (eid == 1024){
			SubdivEvaluation::evaluate(this, eid, xl, yl, p3);
    	} else if (eid == 1536){
			SubdivEvaluation::evaluate(this, eid, xl, yl, p4);
		}	
	}

    kGauss += kGauss_elem;
    kMean  += kMean_elem;
    kGauss_vec[eid] += kGauss_elem / cache.area;
    kMean_vec [eid] += kMean_elem  / cache.area;
    
    // also compute the total area here
    area_sum += cache.area;
  }
  area = area_sum;
  // measure kinetic properties
  Real Ekin = 0;

#ifdef SHELL_OMP
#pragma omp parallel for reduction(+:Ekin)
#endif
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    if (node_cache[nid].is_ghost) continue;
    
    // nodal mass
    mass_vec[nid] = 1 / Minv[mesh.node(nid).dof_number(0, 0, 0)];

    Real Ekin_node = 0;
    for (unsigned int var = 0; var < 3; ++var)
    {
      const unsigned int dof = mesh.node(nid).dof_number(0, var, 0);
      Ekin_node  += v[dof] * v[dof] / Minv[dof];
    }
    Ekin_node *= 0.5;
    Ekin += Ekin_node;
    ekin_vec[nid] += Ekin_node / avg_area_node(nid);
  }
  
  if (!std::isfinite(Ekin))
  {
    register_problem(NOT_FINITE);
    problem_info.set<Real>("Ekin") = Ekin;
  }
  
  // compute the enclosed volume
  compute_enclosed_volume();

  real sphericity = pi

  // the current packing density
  if (shell_cavity_contact)
    density = area * thickness / cavity_volume;

  // add the interesting quantities to the moving average calculators, weighted with dt
  avg_dt.add(dt);
  avg_area.add(dt * area);
  avg_vol.add(dt * volume);
  avg_pcav.add(dt * cavity_force / cavity_surface);
  avg_Nc.add(dt * Nc);
  avg_Nc_cav.add(dt * Nc_cav);
  avg_Ekin.add(dt * Ekin);
  avg_Etens.add(dt * Etens);
  avg_Ebend.add(dt * Ebend);
  avg_Ecav.add(dt * Ecav);
  avg_Eself.add(dt * Eself);
  avg_Efric.add(dt * Efric);
  avg_kGauss.add(dt * kGauss);
  avg_kMean.add(dt * kMean);

  PERFLOG_STOP("measure_properties()");
}

// prints a line with a state summary
void ShellSystem::print_properties()
{
//  if (fulltime > 0.9200 && fulltime < 0.9208){return;}
  char line[512];
  const Real adt = avg_dt.average();
  Real secs_running = wall_time + time_elapsed<Real>(start_time);
  if (secs_running < 1)
    secs_running = std::floor(secs_running * 1e7) * 1e-7; // for sprintf

  time_t rawtime;
  std::time(&rawtime);
  struct tm * timeinfo = localtime(&rawtime);
  char now[13];
  strftime(now, 13, "%y%m%d%H%M%S", timeinfo);

  sprintf(
      line,
      "%12s %10d %11.6g %11.6g %9.8g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g %11.6g % 11.5g % 11.5g",
      now,
      time_step,
      adt,
      this->time,
      secs_running,
      avg_area.average() / adt,
      avg_vol.average() / adt,
      density,
      avg_pcav.average() / adt,
      avg_Nc_cav.average() / adt,
      avg_Nc.average() / adt,
      avg_Ekin.average() / adt,
      avg_Etens.average() / adt,
      avg_Ebend.average() / adt,
      avg_Ecav.average() / adt,
      avg_Eself.average() / adt,
      avg_Efric.average() / adt,
      avg_kGauss.average() / adt,
      avg_kMean.average() / adt);

  *os << "# " << line << std::endl;

#ifdef MEASURE_PERFORMANCE
  perf_log.print_log();
#endif
  
    if (increment_factor < 1)
  	  increment_factor += 0.02;

	std::cout << "Volume growth red " << volume_growth_red << std::endl;
//	std::cout << "Fulltime " << fulltime << std::endl;
}

// returns true if there has been a problem that requires termination, and prints the worst one that occurred
bool ShellSystem::check_for_problems()
{
  // only print the problem if it's not the same as in the last step, to avoid spamming stdout
  if (current_problem != last_problem)
  {
    switch (current_problem)
    {
      case NO_PROBLEM:
        *os << "CLEAR: Problem resolved";
        break;
      case NO_DISP_CHANGE:
        *os << "WARNING: No change in displacement";
        break;
      case LARGE_DISP:
        *os << "WARNING: A physical solution is not guaranteed, because the increment is too large";
        break;
      case SMALL_TIMESTEP:
        *os << "ERROR: The timestep is too small";
        break;
      case TOO_MANY_CELLS:
        *os << "ERROR: The number of linked cells seems to have exploded";
        break;
      case CONTACT_PROBLEM:
        *os << "ERROR: Problem in shell-shell contact";
        break;
      case NOT_FINITE:
        *os << "ERROR: The solution is not finite";
        break;
      default:
        *os << "ERROR: An unknown problem (" << current_problem << ") occurred";
        current_problem = STOP_SIM;
    }

    *os << " in time step " << time_step << "!" << std::endl;

    // print some more helpful info
    if (problem_info.n_parameters() > 0)
    {
      Parameters::const_iterator it = problem_info.begin();
      for (; it != problem_info.end(); ++it)
      {
        *os << it->first << " = ";
        it->second->print(*os);
        *os << std::endl;
      }
    }
    problem_info.clear();
  }

  last_problem = current_problem;
  current_problem = NO_PROBLEM;
  return (last_problem >= STOP_SIM || (stop_on_problem && last_problem > NO_PROBLEM));
}

// writes the current simulation state to a file
void ShellSystem::save(std::ofstream & savefile)
{
  const MeshBase& mesh = get_mesh();

  // save global data
  savefile << dt << "\n";
  savefile << this->time << "\n";
  savefile << wall_time + time_elapsed<Real>(start_time) << "\n";

  // save the shell state
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    const Node& node = mesh.node(nid);
    for (unsigned int var = 0; var < 3; ++var)
    {
      const unsigned int dof = node.dof_number(0, var, 0);
      savefile << u[dof] << "\t"
               << v[dof] << "\t"
               << a[dof] << "\t";
    }
    
    const NodeCache& cache = node_cache[nid];
    savefile << cache.spring_origin(0) << "\t"
             << cache.spring_origin(1) << "\t"
             << cache.spring_origin(2) << "\n";
  }

  // save accumulated measurement data
  savefile << avg_dt.size() << "\n";
  for (unsigned int i = 0; i < avg_dt.size(); ++i)
  {
    savefile << avg_dt(i) << "\t"
             << avg_area(i) << "\t"
             << avg_vol(i) << "\t"
             << avg_pcav(i) << "\t"
             << avg_Nc_cav(i) << "\t"
             << avg_Nc(i) << "\t"
             << avg_Ekin(i) << "\t"
             << avg_Etens(i) << "\t"
             << avg_Ebend(i) << "\t"
             << avg_Ecav(i) << "\t"
             << avg_Eself(i) << "\t"
             << avg_Efric(i) <<"\t"
             << avg_kGauss(i) << "\t"
             << avg_kMean(i) << "\n";
  }

  // save the self-contacts that are currently in static friction mode,
  // because they hold the static friction springs
  unsigned int nc = 0;
  for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
  {
    const ShellShellContact* contact = shell_shell_contact_candidates[i];
    nc += (contact && contact->stick);
  }
  savefile << nc << "\n";
  if (nc > 0)
  {
    for (unsigned int i = 0; i < shell_shell_contact_candidates.size(); ++i)
    {
      const ShellShellContact* contact = shell_shell_contact_candidates[i];
      if (contact && contact->stick)
      {
        savefile << contact->eid[0] << "\t"
                 << contact->eid[1] << "\t"
                 << contact->spring_origin(0) << "\t"
                 << contact->spring_origin(1) << "\t"
                 << contact->spring_origin(2) << "\n";
      }
    }
  }

  // save growth tensors and strains for the computation of the next growth increment
  if (grow_dynamically)
  {
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      if (static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost()) continue;

      for (unsigned int qp = 0; qp < n_qp; ++qp)
      {
        RealVectorValue& alpha = elem_cache[eid]->alpha[qp];
        RealVectorValue& beta  = elem_cache[eid]->beta[qp];
        RealVectorValue& G     = elem_cache[eid]->G[qp];
        savefile << alpha(0) << "\t" << alpha(1) << "\t" << alpha(2) << "\t";
        savefile << beta (0) << "\t" << beta (1) << "\t" << beta (2) << "\t";
        savefile << G(0) << "\t" << G(1) << "\t" << G(2) << "\n";
      }
    }
  }
}

// loads a saved simulation state from a file
void ShellSystem::load(std::ifstream & loadfile)
{
  // load global data
  Real dt_loaded;
  loadfile >> dt_loaded >> this->time >> wall_time;

  // if adaptive time stepping is disabled, allow respecifying the time step
  // in the config file or command line, i.e. don't use the loaded one
  if (adapt_dt_every > 0)
    set_newmark_parameters(dt_loaded);

  // load the shell state
  const MeshBase& mesh = get_mesh();
  for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
  {
    const Node& node = mesh.node(nid);

    for (unsigned int var = 0; var < 3; ++var)
    {
      const unsigned int dof = node.dof_number(0, var, 0);
      loadfile >> u[dof] >> v[dof] >> a[dof];
      u_pred[dof] = u[dof];
    }
    
    NodeCache& cache = node_cache[nid];
    loadfile >> cache.spring_origin(0)
             >> cache.spring_origin(1)
             >> cache.spring_origin(2);
  }

  // load accumulated measurement data
  unsigned int n_measured;
  loadfile >> n_measured;
  for (unsigned int i = 0; i < n_measured; ++i)
  {
    loadfile
    >> avg_dt.add()
    >> avg_area.add()
    >> avg_vol.add()
    >> avg_pcav.add()
    >> avg_Nc.add()
    >> avg_Nc.add()
    >> avg_Ekin.add()
    >> avg_Etens.add()
    >> avg_Ebend.add()
    >> avg_Ecav.add()
    >> avg_Eself.add()
    >> avg_Efric.add()
    >> avg_kGauss.add()
    >> avg_kMean.add();
  }

  // load self-contacts
  unsigned int eid1, eid2, n_load_contacts;
  loadfile >> n_load_contacts;
  for (unsigned int i = 0; i < n_load_contacts; ++i)
  {
    loadfile >> eid1 >> eid2;
    ShellShellContact* contact;
    
    switch (shell_shell_contact)
    {
      case LINEAR:
      {
        contact = new ShellShellContact(eid1, eid2);
        break;
      }
      case SMOOTH:
      {
        contact = new SmoothShellShellContact(eid1, eid2);
        break;
      }
      default:
        ASSERT(shell_shell_contact != INVALID_CT)
    }
    
    loadfile >> contact->spring_origin(0) >> contact->spring_origin(1) >> contact->spring_origin(2);
    contact->stick = true;
    shell_shell_contact_candidates.push_back(contact);
    elem_cache[eid1]->in_contact_with.push_back(eid2);
  }

  // load growth tensors and strains for the computation of the next growth increment
  if (grow_dynamically)
  {
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      if (elem_cache[eid])
      {
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
          RealVectorValue& alpha = elem_cache[eid]->alpha[qp];
          RealVectorValue& beta  = elem_cache[eid]->beta[qp];
          RealVectorValue& G     = elem_cache[eid]->G[qp];
          loadfile >> alpha(0) >> alpha(1) >> alpha(2)
                   >> beta (0) >> beta (1) >> beta (2)
                   >> G(0) >> G(1) >> G(2);
        }
      }
    }
  }

  // don't count loading the program in the continued wall time measurement
  wall_time -= time_elapsed<Real>(start_time);
}

// writes the current state in "Serial vtkUnstructuredGrid" XML format
void ShellSystem::write_vtk(std::string filename)
{
  const MeshBase& mesh = get_mesh();
  unsigned char value_count; // deliberate overflow after 256 values have been written

 // if (fulltime > 0.9200 && fulltime < 0.9208){return;}

/*
  if (time_step > 1)
  {
  std::ostringstream s;
  s << "ux_poisson3_" << shell_young_modulus2 << "_" << press[3] << ".dat"; 
  std::string ux(s.str());


//  std::cout << "press " << press[3] << std::endl;
  file2.open(ux.c_str() , std::fstream::in | std::fstream::out | std::fstream::app);	
  file2 << u[mesh.node(press[0]).dof_number(0, 0, 0)] << "\t" << u[mesh.node(press[1]).dof_number(0, 1, 0)] << std::endl;
//  file2 << u[mesh.node(0).dof_number(0, 0, 0)] << "\t" << u[mesh.node(544).dof_number(0, 0, 0)] << "\t";
//  file2 << u[mesh.node(272).dof_number(0, 1, 0)] << "\t" << u[mesh.node(816).dof_number(0, 1, 0)] << std::endl;
//  file2 << u[mesh.node(0).dof_number(0, 0, 0)] << "\t" << u[mesh.node(529).dof_number(0, 0, 0)] << "\t";
//  file2 << u[mesh.node(257).dof_number(0, 1, 0)] << "\t" << u[mesh.node(801).dof_number(0, 1, 0)] << std::endl;	
//  file2 << u[mesh.node(0).dof_number(0, 0, 0)] << "\t" << u[mesh.node(144).dof_number(0, 0, 0)] << "\t";
//  file2 << u[mesh.node(72).dof_number(0, 1, 0)] << "\t" << u[mesh.node(216).dof_number(0, 1, 0)] << std::endl;
  file2.close();
  }
*/
  std::ofstream file(filename.c_str());
  if (file.is_open())
  {
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << mesh.n_nodes() - n_ghost_nodes << "\" NumberOfCells=\"" << mesh.n_elem() - n_ghost_elem << "\">\n";
    file << "      <PointData>\n";
    
    for (unsigned int var = 0; var < n_vars(); ++var)
    {
      value_count = 0;
      file << "        <DataArray type=\"Float64\" Name=\"" << variable_name(var) << "\" format=\"ascii\">\n";
      for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
      {
        if (!node_cache[nid].is_ghost)
        {
          file << u[mesh.node(nid).dof_number(0, var, 0)] << " ";
          if (++value_count == 0) file << "\n";
        }
      }
      file << "\n";
      file << "        </DataArray>\n";
    }
    
    std::map<std::string, ScalarVector*>::iterator it = node_vectors.begin();
    std::map<std::string, ScalarVector*>::iterator end_it = node_vectors.end();
    for (; it != end_it; ++it)
    {
      value_count = 0;
      file << "        <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"ascii\">\n";
      const ScalarVector& vec = *(it->second);
      for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
      {
        if (!node_cache[nid].is_ghost)
        {
          file << vec[nid] << " ";
          if (++value_count == 0) file << "\n";
        }
      }
      file << "\n";
      file << "        </DataArray>\n";
    }
    
    file << "      </PointData>\n";
    file << "      <CellData>\n";

    it = elem_vectors.begin();
    end_it = elem_vectors.end();
    for (; it != end_it; ++it)
    {
      value_count = 0;
      file << "        <DataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"ascii\">\n";
      const ScalarVector& vec = *(it->second);
      for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
      {
        if (!static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost())
        {
          file << vec[eid] << " ";
          if (++value_count == 0) file << "\n";
        }
      }
      file << "\n";
      file << "        </DataArray>\n";
    }

    file << "      </CellData>\n";
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    value_count = 0;
    for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
    {
      if (!node_cache[nid].is_ghost)
      {
        const Node& node = mesh.node(nid);
        for (unsigned int var = 0; var < 3; ++var)
        {
          file << node(var) << " ";
          if (++value_count == 0) file << "\n";
        }
      }
    }
    file << "\n";
    
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    file << "      <Cells>\n";
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    
    value_count = 0;
    for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
    {
      const Elem* elem = mesh.elem(eid);
      if (!static_cast<const Tri3SD*>(elem)->is_ghost())
      {
        for (unsigned int n = 0; n < 3; ++n)
        {
          file << elem->node(n) << " ";
          if (++value_count == 0) file << "\n";
        }
      }
    }
    file << "\n";

    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    
    value_count = 0;
    const unsigned int max_offset = 3 * (mesh.n_elem() - n_ghost_elem);
    for (unsigned int i = 3; i <= max_offset; i += 3)
    {
      file << i << " ";
      if (++value_count == 0) file << "\n";
    }
    file << "\n";
    
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    
    value_count = 0;
    for (unsigned int i = n_ghost_elem; i < mesh.n_elem(); ++i)
    {
      file << "5 "; // VTK_TRIANGLE
      if (++value_count == 0) file << "\n";
    }
    file << "\n";
    
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
  }
  file.close();
//  write_basis_vtk1(filename);
//  write_basis_vtk2(filename);
}
void ShellSystem::reinit_mod (const Elem* elem, const AutoPtr<FEBase>& fe, std::vector<Real> qw)
{
    const MeshBase& mesh = get_mesh();
    const DofMap& dof_map = get_dof_map();
   
	const Tri3SD* sd_elem = static_cast<const Tri3SD*>(elem);

	libmesh_assert(!sd_elem->is_ghost());
	libmesh_assert(sd_elem->is_subdiv_updated());
	
	// Find basis and transfomation from current basis to preferred basis 
	const std::vector<RealGradient>& dbasis1 = fe->get_dxyzdxi();
	const std::vector<RealGradient>& dbasis2 = fe->get_dxyzdeta();
    const std::vector<Point>&        qpc     = fe->get_xyz();

    ElemCache& cache = *elem_cache[elem->id()];
    std::copy(fe->get_JxW().begin(), fe->get_JxW().end(), cache.JxW);

	Real angle, anglejb1, kfac, jfac;
    RealVectorValue n, dir, j, k, b1, b2, v, pivot, pivotn, normal;
    RealTensorValue A, B, R, Rod, Rn;

	n = dbasis1[0].cross(dbasis2[0]);
	normal = n;
    b1 = dbasis1[0];
    b2 = dbasis2[0];

	// Direction from qp to reference point
//	dir(0) =   0 - qpc[0](0);
//	dir(1) =   0 - qpc[0](1);
//	dir(2) =  10 - qpc[0](2);

	dir(0) =   0;
	dir(1) =   0;
	dir(2) =  100;

	// Normalize direction
	dir = dir / dir.size();
	n = n / n.size();

	// Projection 
	j = dir - (dir*n)*n;
	k = -n.cross(j);

	// n = a * b * sin(a/b)
    anglejb1 = atan2((j.cross(b1)).size(), j*b1);
	jfac = normal.size() / (b1.size() * sin(anglejb1));
	kfac = normal.size() / jfac; 

    // Ensure that the new direction vectors are normalized
    j *= jfac / j.size();
    k *= kfac / k.size();

    // Align normal vector with [0 0 1]
    v(0) = v(1) = 0.0;
 	v(2)        = 1.0;
    pivot = n.cross(v);
    pivotn = pivot/pivot.size();
    angle = atan2(pivot.size(), n*v);
 
 	// Rodrigues 
 	Rod(0,0) = Rod(1,1) = Rod(2,2) = 0;
 	Rod(0,1)                       = -pivotn(2);
 	Rod(0,2)                       =  pivotn(1);
 	Rod(1,0)                       =  pivotn(2);
 	Rod(1,2)                       = -pivotn(0);
 	Rod(2,0)                       = -pivotn(1);
 	Rod(2,1)                       =  pivotn(0);

    // Final rotation matrix
	Rn = Id + sin(angle) * Rod + (1-cos(angle)) * Rod * Rod;

	// Set up A and B matrices
	for (unsigned int i = 0; i < 3; ++i)
	{
    	A(i,0) = b1(i);
    	A(i,1) = b2(i);
   		B(i,0) = k(i);
   		B(i,1) = j(i);
    	A(i,2) = B(i,2) = normal(i);
	}

    const int l = 210;
	if (elem->id() < 2)
	{
//		std::cout << "Qp " << qpc[0] << std::endl;
//		std::cout << "Fac " << jfac << "  " << kfac << std::endl;
//		std::cout << std::endl << "Jacobian " << *cache.JxW << "  " << n.size() << std::endl;
//	    std::cout << "Direction " << dir << std::endl;
//        std::cout << "*************************************************" << std::endl;
//    	std::cout << std::endl << "A and B @3D " << std::endl << A << std::endl << B << std::endl;
//        std::cout << "The angle is " << angle << std::endl;
//        std::cout << "Pivot " << pivotn <<  std::endl;
//		std::cout << "Old basis " << std::endl << b1 << std::endl << b2 << std::endl << normal << std::endl << std::endl;
//		std::cout << "Transformed basis " << std::endl << j << std::endl << k << std::endl << normal << std::endl << std::endl;		
//		std::cout << "Rotation " << Rn << std::endl;
	}

    // Rotate A and B in their 2D space
    A = Rn * A;
    B = Rn * B; 

	// Transformation A - B
    R = inverse(A) * B;

	if (elem->id() < 2)
	{
//		std::cout << std::endl << "A and B @2D " << std::endl << A << std::endl << B << std::endl;
//		std::cout << "Transform " << R(0,0) << " " << R(1,0) << " " << R(0,1) << " " << R(1,1) << std::endl;
	}
	

	// Elem properties
    const unsigned int sf = sd_elem->get_ordered_valence(0) + 6;
    const unsigned int n_qp = qw.size();   
	
	// allocate memory for phi, dphi and d2phi
    cache.phi   = new Real        [sf * n_qp];
    cache.dphi  = new RealGradient[sf * n_qp];
    cache.d2phi = new RealTensor  [sf * n_qp];

	// Get original shape functions
	const std::vector<std::vector<Real > >& phin            = fe->get_phi();
    const std::vector<std::vector<RealGradient > >& dphin   = fe->get_dphi(); 
    const std::vector<std::vector<RealTensor > >& d2phin    = fe->get_d2phi(); 

	// Create new vectors for the shape function transformation
	std::vector<std::vector<Real > > dphidxi(sf,std::vector<Real>(n_qp));   
    std::vector<std::vector<Real > > dphideta(sf,std::vector<Real>(n_qp));
    std::vector<std::vector<Real > > d2phidxi2(sf,std::vector<Real>(n_qp));
    std::vector<std::vector<Real > > d2phideta2(sf,std::vector<Real>(n_qp));
    std::vector<std::vector<Real > > d2phidxideta(sf,std::vector<Real>(n_qp));

	std::vector<std::vector<Real > >         phi(sf,std::vector<Real>(n_qp));
    std::vector<std::vector<RealGradient > > dphi(sf,std::vector<RealGradient>(n_qp));
    std::vector<std::vector<RealTensor > >   d2phi(sf,std::vector<RealTensor>(n_qp));

 	// copy JxW from a std::vector to a C array for higher efficiency
	std::copy(fe->get_xyz().begin(), fe->get_xyz().end(), cache.qxyz);

    // Assign new dphi and d2phi
	for (unsigned int p = 0; p < n_qp; ++p )
	{
    	for ( unsigned int w = 0; w < sf; ++w )
    	{
			dphidxi[w][p]      = dphin[w][p](0) * R(0,0) + dphin[w][p](1) * R(1,0);  
     		dphideta[w][p]     = dphin[w][p](0) * R(0,1) + dphin[w][p](1) * R(1,1);
 
        	dphi[w][p](0)      = dphidxi[w][p];
        	dphi[w][p](1)      = dphideta[w][p];	

			d2phidxi2[w][p]    = d2phin[w][p](0,0) * R(0,0) * R(0,0) + d2phin[w][p](1,1) * R(1,0) * R(1,0) + 2 * d2phin[w][p](1,0) * R(0,0) * R(1,0);
            d2phideta2[w][p]   = d2phin[w][p](0,0) * R(0,1) * R(0,1) + d2phin[w][p](1,1) * R(1,1) * R(1,1) + 2 * d2phin[w][p](0,1) * R(0,1) * R(1,1);
            d2phidxideta[w][p] = d2phin[w][p](0,0) * R(0,0) * R(0,1) + d2phin[w][p](1,1) * R(1,0) * R(1,1) + d2phin[w][p](0,1) * R(0,1) * R(1,0) + d2phin[w][p](1,0) * R(0,0) * R(1,1);

			d2phi[w][p](0,0)   = d2phidxi2[w][p];
			d2phi[w][p](1,1)   = d2phideta2[w][p];
        	d2phi[w][p](0,1)   = d2phi[w][p](1,0) = d2phidxideta[w][p];	
	
	        // copy the data from std::vectors to C arrays for higher efficiency
    	    std::copy(fe->get_phi ()[w].begin(), fe->get_phi ()[w].end(), &cache.phi [w*n_qp]);
       		std::copy(dphi [w].begin(), dphi [w].end(), &cache.dphi [w*n_qp]);
        	std::copy(d2phi[w].begin(), d2phi[w].end(), &cache.d2phi[w*n_qp]);
		}
	}
	
	// Let the FESubdivMap compute the new basis vectors
	fe->_fe_map->get_phi_map()          = phi;
	fe->_fe_map->get_dphidxi_map()      = dphidxi;
	fe->_fe_map->get_dphideta_map()     = dphideta;
	fe->_fe_map->get_d2phidxi2_map()    = d2phidxi2;
	fe->_fe_map->get_d2phideta2_map()   = d2phideta2;
	fe->_fe_map->get_d2phidxideta_map() = d2phidxideta;
  
//	#pragma omp critical
	// Compute basis map
    fe->_fe_map->compute_map (2, qw, elem);

}

// writes the current state in "Serial vtkUnstructuredGrid" XML format
void ShellSystem::write_basis_vtk1(std::string filename)
{
  const MeshBase& mesh = get_mesh();
  unsigned char value_count; // deliberate overflow after 256 values have been written

  filename.replace(filename.end()-4,filename.end(),"quad1.vtu");
  std::ofstream file(filename.c_str());


  const unsigned int max_offset = (mesh.n_elem() - n_ghost_elem);
  if (file.is_open())
  {
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << max_offset << "\" NumberOfCells=\"" << max_offset << "\">\n";
    file << "      <PointData>\n";

      unsigned int var = 0;	
      value_count = 0;
      file << "        <DataArray type=\"Float64\" Name=\"Basis" << var << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
      {
 	    const Elem* elem = mesh.elem(eid);
        if (!static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost())
        {
	      const ElemCache& cache = *elem_cache[eid];
          const SurfaceMetric& a_def = cache.a_def[0];
          if (eid%2 == 0)
    	  {
          	file << a_def.vec[var](0) << " " << a_def.vec[var](1) << " " << a_def.vec[var](2) << " ";
          } else {
          	file << a_def.vec[1](0) << " " << a_def.vec[1](1) << " " << a_def.vec[1](2) << " ";
		  } 
          if (++value_count == 0) file << "\n";		  
     	}
      }  
    file << "        </DataArray>\n";        

    file << "      </PointData>\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";    
    value_count = 0;
    for (unsigned int eid = 0; eid < max_offset; ++eid)
    {
      const Elem* elem = mesh.elem(eid);
      if (!static_cast<const Tri3SD*>(elem)->is_ghost())
      { 
      	const ElemCache& cache = *elem_cache[eid];
	  	const SurfaceMetric& a_def = cache.a_def[0];
      	for (unsigned int var = 0; var < 3; ++var)
        {
          file << a_def.xyz(var) << " ";
          if (++value_count == 0) file << "\n";
		}
      }    
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

	file << "      <Cells>\n";
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
   
   	value_count = 0;
    for (unsigned int j = 0; j < max_offset; ++j)
	{
	  file << j << " ";
      if (++value_count == 0) file << "\n";
	}
    file << "\n";

    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    value_count = 0;

    for (unsigned int i = 1; i < (max_offset+1); ++i)
	{
	  file << i << " ";
      if (++value_count == 0) file << "\n";   
	}
 
    file << "\n";
    value_count = 0;
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
 
 	for (unsigned int j = 0; j < max_offset; ++j)
	{
	  file << "1" << " ";
       if (++value_count == 0) file << "\n";   
	}
   
    file << "\n";
    
    file << "        </DataArray>\n";
    file << "      </Cells>\n";


    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
  }
  file.close();
}

// writes the current state in "Serial vtkUnstructuredGrid" XML format
void ShellSystem::write_basis_vtk2(std::string filename)
{
  const MeshBase& mesh = get_mesh();
  unsigned char value_count; // deliberate overflow after 256 values have been written

  filename.replace(filename.end()-4,filename.end(),"quad2.vtu");
  std::ofstream file(filename.c_str());


  const unsigned int max_offset = (mesh.n_elem() - n_ghost_elem);
  if (file.is_open())
  {
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << max_offset << "\" NumberOfCells=\"" << max_offset << "\">\n";
    file << "      <PointData>\n";

      unsigned int var = 0;	
      value_count = 0;
      file << "        <DataArray type=\"Float64\" Name=\"Basis" << var << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
      for (unsigned int eid = 0; eid < mesh.n_elem(); ++eid)
      {
 	    const Elem* elem = mesh.elem(eid);
        if (!static_cast<const Tri3SD*>(mesh.elem(eid))->is_ghost())
        {
	      const ElemCache& cache = *elem_cache[eid];
          const SurfaceMetric& a_def = cache.a_def[0];
          if (eid%2 == 1)
    	  {
          	file << a_def.vec[var](0) << " " << a_def.vec[var](1) << " " << a_def.vec[var](2) << " ";
          } else {
          	file << 0 << " " << 0 << " " << 0 << " ";
		  } 
          if (++value_count == 0) file << "\n";		  
     	}
      }  
    file << "        </DataArray>\n";        

    file << "      </PointData>\n";

    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";    
    value_count = 0;
    for (unsigned int eid = 0; eid < max_offset; ++eid)
    {
      const Elem* elem = mesh.elem(eid);
      if (!static_cast<const Tri3SD*>(elem)->is_ghost())
      { 
      	const ElemCache& cache = *elem_cache[eid];
	  	const SurfaceMetric& a_def = cache.a_def[0];
      	for (unsigned int var = 0; var < 3; ++var)
        {
          file << a_def.xyz(var) << " ";
          if (++value_count == 0) file << "\n";
		}
      }    
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

	file << "      <Cells>\n";
    file << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
   
   	value_count = 0;
    for (unsigned int j = 0; j < max_offset; ++j)
	{
	  file << j << " ";
      if (++value_count == 0) file << "\n";
	}
    file << "\n";

    file << "        </DataArray>\n";
    file << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    value_count = 0;

    for (unsigned int i = 1; i < (max_offset+1); ++i)
	{
	  file << i << " ";
      if (++value_count == 0) file << "\n";   
	}
 
    file << "\n";
    value_count = 0;
    file << "        </DataArray>\n";
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
 
 	for (unsigned int j = 0; j < max_offset; ++j)
	{
	  file << "1" << " ";
       if (++value_count == 0) file << "\n";   
	}
   
    file << "\n";
    
    file << "        </DataArray>\n";
    file << "      </Cells>\n";


    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
  }
  file.close();
}

