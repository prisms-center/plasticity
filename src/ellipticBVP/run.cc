//run method for ellipticBVP class

#ifndef RUN_ELLIPTICBVP_H
#define RUN_ELLIPTICBVP_H
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
  computing_timer.exit_section("mesh and initialization");

  //assembly
  computing_timer.enter_section("assembly");
  assemble();
  computing_timer.exit_section("assembly");

  //solve();
  solve();
}

#endif
