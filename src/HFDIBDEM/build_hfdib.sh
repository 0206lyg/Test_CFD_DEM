#!/bin/bash

#SBATCH --job-name=hfdib_build
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G

module purge
module load gcc/11.4.1
module load openmpi/4.1.6
source "$HOME/OpenFOAM/OpenFOAM-8/etc/bashrc"

echo "=== Building HFDIBDEM library ==="
cd /home/lee6456/apps/lyjHFDIB-DEM/openHFDIB-DEM/src/HFDIBDEM

wclean libso || true
rm -rf Make/linux64GccDPInt32Opt
rm -f "$FOAM_USER_LIBBIN/liblyjHFDIBDEM.so"
wmake libso

echo "=== Building pimpleLYJHFDIBFoam solver ==="
cd /home/lee6456/apps/lyjHFDIB-DEM/openHFDIB-DEM/applications/solvers/incompressible/pimpleLYJHFDIBFoam
wclean || true
rm -rf Make/linux64GccDPInt32Opt
rm -f "$FOAM_USER_APPBIN/pimpleLYJHFDIBFoam"
wmake

echo "=== Build completed successfully ==="
which pimpleLYJHFDIBFoam || true
ls -lh "$FOAM_USER_LIBBIN/liblyjHFDIBDEM.so" || true
ls -lh "$FOAM_USER_APPBIN/pimpleLYJHFDIBFoam" || true
