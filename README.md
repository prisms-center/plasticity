<B>What is PRISMS Plasticity?</B>

  It is a Finite Element Method (FEM) code for solving boundary value 
  problems arising in Continuum Plasticity and Crystal Plasticity. 
  It is build on top of the deal.II open source finite element 
  library [http://www.dealii.org]
  
  This code is developed by the PRedictive Integrated Structural
  Materials Science (PRISMS) Center [http://www.prisms-center.org/]
  at University of Michigan which is supported by the U.S. Department 
  of Energy (DOE), Office of Basic Energy Sciences, Division of Materials Sciences 
  and Engineering under Award #DE-SC0008637

<B>Installation and Use:</B>

  1) Install deal.II (version 8.2.1)<br>
     Download Binaries (OSX and Linux) or  Virtual Machine (VMI) from https://www.dealii.org/download.html <br>
     (OR) <br>
     Configure, compile and install the deal.II library with the 
     following configuration flags. Dependencies are MPI, p4est,
     PETSc libraries.<br>
     -DDEAL_II_WITH_MPI=ON, -DDEAL_II_WITH_LAPACK=ON, -DDEAL_II_WITH_P4EST=ON, -DDEAL_II_WITH_PETSC=ON

     Download: http://www.dealii.org/download.html <br>
     Installation instructions: http://www.dealii.org/8.2.1/readme.html <br>
     Installation instructions for external packages (P4EST, PETSC): https://www.dealii.org/developer/external-libs/

  2) Fork the repo https://github.com/prisms-center/plasticity on
  GitHub (How to fork: https://help.github.com/articles/fork-a-repo/)
  and clone the plasticity repository using your GitHub username:<br>
  $ git clone git@github.com:username/plasticity.git <br>
  (OR) <br>
  $ git clone https://github.com/username/plasticity.git <br>
  and <br>
  $ cd plasticity <br>
  $ git checkout master <br>
   
  3) Running plasticity applications, for example simple tension
  boundary value problem with continuum plasticity material model: 
  $ cd applications/continuumPlasticity/simpleTension
  For debug mode [default mode]:
  $ cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Debug
  For optimized mode:
  $ cmake CMakeLists.txt -DCMAKE_BUILD_TYPE=Release 
  and
  $ make  
  Execution (serial runs):
  $ make run
  Execution (parallel runs):
  $ mpiexec -np nprocs ./main
  [here nprocs denotes the number of processors]
  
  4) Updates: Since plasticity code is still under active development,
  regular code and documentation updates are pushed to the upstream
  repo (https://github.com/prisms-center/plasticity) and we strongly
  recommend users to synchronize their respective forks at regular
  intervals or when requested by the developers through the
  announcements on the mailing list. 
  (How to sync: https://help.github.com/articles/syncing-a-fork/)

  5) Visualization: Output of the primal fields is in standard vtk 
  format (*.pvtu, *.vtu files) which can be visualized with the 
  following open source applications:
  1. VisIt (https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
  2. Paraview (http://www.paraview.org/download/)

<B>Getting started:</B>

  Examples of various boundary value problems are located under the 
  applications/ folder. Easiest way to get started on the code is to 
  run the applications.

  Applications are intended to serve as (1) Demonstration of the
  capabilities of this library, (2) Provide a framework for
  further development of specialized/advanced applications by
  users. 

  Application or code under development/testing is preceded by an
  underscore, such as:
  _hcp, _bcc

  List of folders:
  src/: material models (continuum plasticity and crystal plasticity), 
  ellipticBVP (base class for parallel implementation of elliptic 
  boundary value problems), enrichment models (enhanced strain), 
  utility (static condensation class, crystal orientations and grain 
  information IO)
  applications/: continuum plasticity and crystal plasticity
  utils/: IntegrationTools (developed by the PRISMS center and available at
  https://github.com/prisms-center/IntegrationTools) and json headers 

<B>Documentation:</B>

  Under development. 
  Details of material model formulation are available under src/materialmodels
  e.g. src/materialModels/crystalPlasticity/CPFEM.pdf   
 	
<B>License:</B>

  GNU Lesser General Public License (LGPL). Please see the file
  LICENSE for details.

<B>Mailing List:</B>
  
  PRISMS.plasticity@umich.edu

<B>Further information, questions, issues and bugs:</B>

  Contact the plasticity mailing list at PRISMS.plasticity@umich.edu  



