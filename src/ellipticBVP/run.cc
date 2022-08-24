//run method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <sys/stat.h>

//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void ellipticBVP<dim>::run(){

  const int dir_err = mkdir(userInputs.outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //initialization
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  computing_timer.enter_section("mesh and initialization");
  #else
  computing_timer.enter_subsection("mesh and initialization");
  #endif
  //read mesh;
  mesh();
  //initialize FE objects and global data structures
  init();
  initProjection();
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  computing_timer.exit_section("mesh and initialization");
  #else
  computing_timer.leave_subsection("mesh and initialization");
  #endif

  solve();
}
#include "../../include/ellipticBVP_template_instantiations.h"
