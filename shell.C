#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <csignal>

#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_modification.h"
#include "mesh_refinement.h"
#include "mesh_subdiv_support.h"
#include "boundary_info.h"
#include "equation_systems.h"
#include "fe_base.h"
#include "getpot.h"

#include "tools.h"
#include "ShellSystem.h"

//====================================================================================================================

enum ExitCode { UNFINISHED = -1, SUCCESS = 0, TIME_LIMIT_REACHED, SIGNAL_CAUGHT, ERROR };

ExitCode exit_code = UNFINISHED;

void catchSignal(int sig)
{
  char now [20];
  timeNow(now);
  std::cout << "Received signal " << sig << ": " << strsignal(sig) << " (" << now << ")" << std::endl;
  exit_code = (SIGNAL_CAUGHT > exit_code ? SIGNAL_CAUGHT : exit_code);
}

//====================================================================================================================

void print_help(char* prog_name, std::ostream& os = std::cout)
{
  os << "Usage: " << prog_name << " [-c configfile] [-l loadfile] [-m meshfile] [-d datadir] [options]\n";
}

//====================================================================================================================

int main(int argc, char** argv)
{
  time_t start_time(time(NULL));
  
  // parse command line arguments
  GetPot command_line(argc, argv);
  
  // print help and exit if requested
  if (command_line.search("-h") || command_line.search("--help"))
  {
    print_help(argv[0], std::cout);
    return SUCCESS;
  }

  LibMeshInit init(argc, argv);

  { // make sure stack is emptied before LibMesh::close() is called

    // print some info for restorability (should be after libMesh initialization, for MPI)
    char now[20];
    timeNow(now);
    std::cout << "Binary was compiled on " << __DATE__ << " at " << __TIME__ << std::endl;
    std::cout << "Invocation at " << now << ": " << argv[0];
    for (int i = 1; i < argc; ++i)
      std::cout << " " << argv[i];
    std::cout << std::endl;
    
    std::ostream* os = &std::cout;
    ASSERT(n_processors() == 1)

    // bind signal handlers
    signal(SIGHUP, SIG_IGN); // for easy running through ssh
    signal(SIGQUIT, catchSignal);
    signal(SIGTERM, catchSignal); // important for job killing on a cluster

    // parse config file
    std::string cfgname = "config.cfg";
    if (command_line.search("-c"))
      cfgname = command_line.next(cfgname);
    
    // make sure the config file can be read
    std::ifstream cfgfile(cfgname.c_str());
    if (!cfgfile)
    {
      std::cout << "ERROR: Can't open config file \"" << cfgname << "\" for reading!\n";
      return ERROR;
    }
    cfgfile.close();
    
    GetPot config_file(cfgname);

    // create the data directory in case it doesn't exist
    std::string data_dir = "data";
    if (command_line.search("-d"))
      data_dir = command_line.next(data_dir);
    std::system((std::string("mkdir -p '") + data_dir + "'; cp '" + cfgname + "' '" + data_dir + "'").c_str());

    // build data file names
    std::string vtk_file_path_format = data_dir + "/shell%010u.vtu";
    std::string state_filename_format = "shell_state%010u.dat";
    std::string latest_state_filename = data_dir + "/shell_state_latest.dat";
    char file_path_buffer [1024];

    // determine the name of the load file, if any
    std::string load_filename;
    if (command_line.search("-l"))
      load_filename = command_line.next(load_filename);

    // read parameters
    MAKE_PARSE(unsigned int, time_steps, 0)
    MAKE_PARSE(double, walltime_limit, 0)
    MAKE_PARSE(double, data_time,  0.01)
    MAKE_PARSE(double, vtk_time,   0.01)
    MAKE_PARSE(double, state_time, 0   )
    ASSERT(walltime_limit >= 0)
    ASSERT(data_time  >= 0)
    ASSERT(vtk_time   >= 0)
    ASSERT(state_time >= 0)

    unsigned int time_step = 0;

    // .. and also load some data if a loadfile was specified
    std::ifstream loadfile;

    if (!load_filename.empty())
    {
      loadfile.open(load_filename.c_str());
      if (!loadfile.is_open())
        std::cout << "WARNING: Can't open load file \"" << load_filename << "\" for reading! Starting new simulation.\n";
    }

    if (loadfile.is_open())
    {
      loadfile >> time_step;
      if (time_steps && time_step >= time_steps)
      {
        std::cout << "ERROR: Loaded state is already finished!\n";
        loadfile.close();
        return SUCCESS;
      }
    }
    
    // create the mesh
    Mesh mesh(2);
    
    MAKE_PARSE(unsigned int, refinements, 0)
    
    if (command_line.search("-m"))
    {
      // make sure the mesh file can be read
      std::string mesh_filename = "";
      mesh_filename = command_line.next(mesh_filename);
      std::ifstream meshfile(mesh_filename.c_str());
      if (!meshfile)
      {
        std::cout << "ERROR: Can't open mesh file \"" << mesh_filename << "\" for reading!\n";
        return ERROR;
      }
      meshfile.close();
      
      // copy the mesh file into the data folder for later reproducibility
      std::system((std::string("cp '") + mesh_filename + "' '" + data_dir + "'").c_str());
      
      mesh.read(mesh_filename);
      
      // refine the mesh 'refinements' times
      if (refinements > 0)
      {
        MeshRefinement mesh_refinement(mesh);
        mesh_refinement.uniformly_refine(refinements);
        MeshTools::Modification::flatten(mesh); // remove refinement data
      }
    }
    else
    {
      // if no shell mesh file is given, build a delaunay-triangulated square with unit half diagonal
      const unsigned int nx = (1u << refinements);
      const Real L_half = std::sqrt(0.5);
      MeshTools::Generation::build_delaunay_square(mesh, nx, nx, -L_half, L_half, -L_half, L_half, TRI3);
      mesh.boundary_info->clear();
      
      // refine the mesh once to remove triangles with more than one irregular vertex
      MeshRefinement mesh_refinement(mesh);
      mesh_refinement.uniformly_refine(1);
      MeshTools::Modification::flatten(mesh); // remove refinement data
    }
    
    MAKE_PARSE(Real, randomness, 0)
    if (randomness != 0)
      MeshTools::Modification::distort(mesh, randomness);
    
    // project the mesh onto an ellipsoid, (elliptic) cylinder or plane if the projection vector is nonzero
    MAKE_VECTOR_PARSE(project, 0);
    if (project[0] != 0 || project[1] != 0 || project[2] != 0)
    {
      for (unsigned int nid = 0; nid < mesh.n_nodes(); ++nid)
      {
        Point n;
        for (unsigned int var = 0; var < 3; ++var)
          if (project[var] != 0)
            n(var) = mesh.node(nid)(var) / project[var];
        const Real factor = 1 / n.size();
        for (unsigned int var = 0; var < 3; ++var)
          if (project[var] != 0)
            mesh.node(nid)(var) *= factor;
      }
    }

    MAKE_VECTOR_PARSE(scale, 1);
    ASSERT(scale[0] != 0)
    ASSERT(scale[1] != 0)
    ASSERT(scale[2] != 0)
    if (scale[0] != 1 || scale[1] != 1 || scale[2] != 1)
      MeshTools::Modification::scale(mesh, scale[0], scale[1], scale[2]);

    MAKE_VECTOR_PARSE(rotate, 0);
    if (fmod(rotate[0],1) != 0 || fmod(rotate[1],1) != 0 || fmod(rotate[2],1) != 0)
      MeshTools::Modification::rotate(mesh, rotate[0]*360, rotate[1]*360, rotate[2]*360);

    MAKE_VECTOR_PARSE(translate, 0);
    if (translate[0] != 0 || translate[1] != 0 || translate[2] != 0)
      MeshTools::Modification::translate(mesh, translate[0], translate[1], translate[2]);

    // turn the mesh into a subdivision surface
    MAKE_PARSE(bool, add_ghosts, true)
    MeshTools::Modification::all_tri(mesh); // convert all elements to triangles
    MeshTools::Subdiv::prepare_subdiv_mesh(mesh, !add_ghosts); // add ghosts, convert to subdivision elements

    // create the system
    EquationSystems equation_systems(mesh);
    ShellSystem& system = equation_systems.add_system<ShellSystem>("Shell");

    // initialize the data structures for the equation system
    equation_systems.init();
    // parse the config file and command line, load a saved state, and apply initial conditions
    system.set_parameters(config_file, command_line);
	
    // warn about unused parameters
    std::vector<std::string> ufos = config_file.unidentified_variables();
    if (ufos.size() > 0)
    {
      std::cout << INFO << "WARNING: Unrecognized config file parameter(s):";
      for (unsigned int i = 0; i < ufos.size(); ++i)
        std::cout << " " << ufos[i];
      std::cout << std::endl;
    }
    ufos = command_line.unidentified_arguments();
    if (ufos.size() > 0)
    {
      std::cout << INFO << "WARNING: Unrecognized command line argument(s):";
      for (unsigned int i = 0; i < ufos.size(); ++i)
        std::cout << " " << ufos[i];
      std::cout << std::endl;
    }
    // determine the element size distribution
    Real hmax = 0, havg = 0, hmin = std::numeric_limits<Real>::infinity();
    MeshBase::const_element_iterator el = mesh.elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.elements_end();
    for (; el != end_el; ++el)
    {
      const Elem* elem = *el;
      const Real this_hmin = elem->hmin();
      const Real this_hmax = elem->hmax();
      if (this_hmin < hmin) hmin = this_hmin;
      if (this_hmax > hmax) hmax = this_hmax;
      havg += (this_hmin + this_hmax) / 2;
    }
    havg /= mesh.n_elem();
    
    // print mesh info
    std::cout << "\nWITHOUT GHOSTS:\n"
              << "nodes   : " << (mesh.n_nodes() - system.get_n_ghost_nodes()) << "\n"
              << "elements: " << (mesh.n_elem() - system.get_n_ghost_elem()) << "\n"
              << "DOFs    : " << 3 * (mesh.n_nodes() - system.get_n_ghost_nodes()) << "\n";

    std::cout << "WITH GHOSTS:\n"
              << "nodes   : " << mesh.n_nodes() << "\n"
              << "elements: " << mesh.n_elem() << "\n"
              << "DOFs    : " << 3 * mesh.n_nodes() << "\n";

    std::cout << "ELEMENT SIZE:\n"
              << "min: " << hmin << "\n"
              << "avg: " << havg << "\n"
              << "max: " << hmax << "\n" << std::endl;

    std::cout.precision(std::numeric_limits<Real>::digits10 + 2);

    // load a saved state and/or apply initial conditions
    system.apply_initial_conditions(loadfile);

    if (data_time > 0)
      system.print_header();

    if (loadfile.is_open())
    {
      loadfile.close();
    }
    else
    {
      // measure and print the initial configuration
      if (data_time > 0)
      {
        system.measure_properties();
        system.print_properties();
      }

      // write an output file for the initial configuration
      if (vtk_time > 0)
      {
		sprintf(file_path_buffer, vtk_file_path_format.c_str(), time_step);
		system.write_vtk(file_path_buffer);
	  }
    }
    
    bool problem = false;
    unsigned int last_data_counter  = (data_time > 0 ? static_cast<unsigned int>(system.time / data_time) : 0);
    unsigned int last_vtk_counter   = (vtk_time  > 0 ? static_cast<unsigned int>(system.time / vtk_time ) : 0);
    unsigned int last_state_counter = 0;

    // integrate in time
    while (true)
    {
      if (system.check_stopping_criterion())
      {
        // exit normally if the stopping criterion is met
        exit_code = SUCCESS;
      }
      else if (walltime_limit > 0 && difftime(time(NULL), start_time) >= walltime_limit)
      {
        // exit semi-normally if the maximum wall time has been reached
        exit_code = TIME_LIMIT_REACHED;
      }
      else
      {
        // advance one step in time
        ++time_step;
        problem = system.advance(time_step);
    
    	// print data at regular time intervals
    	if (system.time > 2) vtk_time = 5e-4;
    
    	/*
    	if (system.time > 0.289)
		{
    		data_time = 3e-5;
    		vtk_time = 3e-5;
		}
        */

        const unsigned int data_counter = (data_time > 0 ? static_cast<unsigned int>(system.time / data_time) : 0);
        if (data_counter > last_data_counter)
        {
          last_data_counter = data_counter;
          system.print_properties();
        }

        // write a vtk file at regular time intervals
		const unsigned int vtk_counter = (vtk_time > 0 ? static_cast<unsigned int>(system.time / vtk_time) : 0);
        if (vtk_counter > last_vtk_counter)
        {
          last_vtk_counter = vtk_counter;
          sprintf(file_path_buffer, vtk_file_path_format.c_str(), time_step);
          system.write_vtk(file_path_buffer);
        }

        // exit normally if the maximum simulation time has been reached
        if (time_steps && time_step >= time_steps)
          exit_code = SUCCESS;

        // exit abnormally if a serious problem has occurred
        if (problem)
          exit_code = ERROR;
      }

      // write a state file at regular time intervals and at the end of the simulation
      const unsigned int state_counter = (state_time > 0 ? static_cast<unsigned int>(time_elapsed<double>(system.get_start_time()) / state_time) : 0);
      if (exit_code > UNFINISHED || state_counter > last_state_counter)
      {
        last_state_counter = state_counter;
        sprintf(file_path_buffer, state_filename_format.c_str(), time_step);
        std::ofstream savefile((data_dir + "/" + file_path_buffer).c_str());
        if (savefile.is_open())
        {
          savefile << std::setprecision(std::numeric_limits<Real>::digits10 + 2);
          savefile << time_step << "\n";
          system.save(savefile);
          savefile.close();
          std::system((std::string("ln -sf '") + file_path_buffer + "' '" + latest_state_filename + "'").c_str());
        }
        if (exit_code > UNFINISHED) break;
      }
    }
  }

  return exit_code;
}
