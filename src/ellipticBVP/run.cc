//run method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void ellipticBVP<dim>::run(){
  //initialization
  computing_timer.enter_section("mesh and initialization");
  //read mesh;
  mesh();
  //initialize FE objects and global data structures
  init();
  initProjection();
  //user model related variables and methods
#ifdef enableUserModel
  initQuadHistory();
#endif

  computing_timer.exit_section("mesh and initialization");

  //solve();
  solve();
}
