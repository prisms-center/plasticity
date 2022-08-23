//run method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <sys/stat.h>

//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void ellipticBVP<dim>::run(){

  const int dir_err = mkdir(userInputs.outputDirectory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //initialization
  computing_timer.enter_section("mesh and initialization");
  //read mesh;
  mesh();
  //initialize FE objects and global data structures
  init();
  initProjection();

  computing_timer.exit_section("mesh and initialization");
  pcout<<"sanity check in ellipticBVP<dim>::run\n";
  //solve();
  solve();
}
#include "../../include/ellipticBVP_template_instantiations.h"
