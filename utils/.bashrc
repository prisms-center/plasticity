# .bashrc

# User specific aliases and functions

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
module load mkl
module load cmake
module load hdf5
alias qprisms="/scratch/prismsproject_flux/bpuchala/Public/scripts/prisms_usage.sh"
#export PETSC_DIR2=/nfs/mcfs_home/rudraa/Public/petsc/petsc-3.4.3
#export PETSC_DIR=/home/rudraa/software/numerics/petsc-3.3-p6
#export PETSC_DIR=/nfs/mcfs_home/rudraa/Public/petsc/petsc-3.4.3HypreSuperLU
export PETSC_DIR=/nfs/mcfs_home/rudraa/Public/petsc/petsc
export SLEPC_DIR=/nfs/mcfs_home/rudraa/Public/slepc/slepc-3.5.1
export PETSC_ARCH=shared-optimized
export TRILINOS_DIR=/nfs/mcfs_home/rudraa/Public/trilinos/ver10.12.2/install 
export TRILINOS_DIR2=/home/rudraa/software/numerics/trilinos-10.12.2/install
export PETIGA_DIR=/nfs/mcfs_home/rudraa/Public/petiga/PetIGA
export LD_LIBRARY_PATH=$PETIGA_DIR/$PETSC_ARCH/lib:$PETSC_DIR/$PETSC_ARCH/lib:$TRILINOS_DIR/lib:$LD_LIBRARY_PATH
#export DEAL_II_DIR=/nfs/mcfs_home/rudraa/Public/deal.II/ver8.1.0/installAvx2
export DEAL_II_DIR=/nfs/mcfs_home/rudraa/Public/deal.II/ver8.1.0/installAvxNew2


