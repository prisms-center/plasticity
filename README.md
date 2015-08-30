What is PRISMS-Plasticity?

  It is a Finite Element Method (FEM) code for solving boundary value 
  problems arising in Continuum Plasticity and Crystal Plasticity. 
  It is build on top of the deal.II open source finite element 
  library (http://www.dealii.org)
  
  This code is developed by the PRedictive Integrated Structural
  Materials Science (PRISMS) Center at University of Michigan which is
  supported by the U.S. Department of Energy, Office of Basic Energy
  Sciences, Division of Materials Sciences and Engineering under Award
  #DE-SC0008637

Installation and Use:

  1) Configure, compile and install the deal.II library with the 
     following configuration flags. Dependencies are MPI, p4est,
     PETSc and Trilinos libraries.
     -DDEAL_II_WITH_MPI=ON  
     -DDEAL_II_WITH_LAPACK=ON
     -DDEAL_II_WITH_P4EST=ON
     -DDEAL_II_WITH_PETSC=ON
     -DDEAL_II_WITH_TRILINOS=ON

     Download: http://www.dealii.org/download.html
     Installation instructions: http://www.dealii.org/8.1.0/readme.html
     Installation instructions for external packages (P4EST, PETSC, 
     TRILINOS): https://www.dealii.org/developer/external-libs/

  2) Fork the repo https://github.com/prisms-center/plasticity on
  GitHub (How to fork: https://help.github.com/articles/fork-a-repo/)
  and clone the plasticity repository using your GitHub username:
  $ git clone git@github.com:username/plasticity.git 
  (OR)
  $ git clone https://github.com/username/plasticity.git
  and
  $ cd plasticity  
  $ git checkout master
  [Note: plasticity is currently a private repository on GitHub and
  hence you need to be authorised to access the repository. Contact
  the developers/mailing-list to request access.]
   
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

Getting started:

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

Documentation:

  Under development. 
  Details of material model formulation are available under src/materialmodels
  e.g. src/materialModels/crystalPlasticity/CPFEM.pdf   
 	
License:

  GNU Lesser General Public License (LGPL). Please see the file
  LICENSE for details.

Mailing List:
  
  PRISMS.plasticity@umich.edu

Further information, questions, issues and bugs:

  Contact the plasticity mailing list at PRISMS.plasticity@umich.edu  


