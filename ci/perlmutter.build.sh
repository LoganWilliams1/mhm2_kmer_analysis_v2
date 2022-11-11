#!/bin/bash

# CI Build script for perlmutter

set -e

uname -a
pwd
date

if [ -z "${MHM2_SOURCE}" ] || [ -z "${CI_SCRATCH}" ]
then
  echo "please set the MHM2_SOURCE and CI_SCRATCH environmental variables"
  exit 1
fi

export INSTALL_PREFIX=${CI_SCRATCH}/install
export BUILD_PREFIX=${CI_SCRATCH}/build
export RUN_PREFIX=${CI_SCRATCH}/runs
export TMPDIR=${CI_SCRATCH}/tmp
mkdir -p ${CI_SCRATCH}/tmp

cd ${MHM2_SOURCE}
git describe
cd upcxx-utils
git describe
cd -

(
echo "Loading GPU environment"
envir=${MHM2_SOURCE}/contrib/environments/perlmutter/gnu.sh
source $envir
module list
env | grep '\(GASNET\|FI_\|UPC\)'
env | grep SLURM || true
env | grep TMP

df -h $TMPDIR

for t in Debug RelWithDebug RelWithDebInfo Release
do
  export MHM2_BUILD=${BUILD_PREFIX}-gpu-$t
  mkdir -p ${MHM2_BUILD}
  cd $MHM2_BUILD
  echo "Building gpu GNU ${t}"
  CXX=CC cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}-gpu-${t} -DCMAKE_BUILD_TYPE=${t} ${MHM2_CMAKE_EXTRAS} ${MHM2_SOURCE}
  make -j 1 all check install
  cp -p $envir ${INSTALL_PREFIX}-gpu-${t}/env.sh # store environment to support execution
done
)

(
echo "Loading CPU-only environment"
envir=${MHM2_SOURCE}/contrib/environments/perlmutter/gnu-cpuonly.sh
source $envir
module list
env | grep '\(GASNET\|FI_\|UPC\)'

for t in Debug RelWithDebug RelWithDebInfo Release
do
  export MHM2_BUILD=${BUILD_PREFIX}-cpu-$t
  mkdir -p ${MHM2_BUILD}
  cd $MHM2_BUILD
  echo "Building cpu GNU ${t}"
  CXX=CC cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}-cpu-${t} -DCMAKE_BUILD_TYPE=${t} ${MHM2_CMAKE_EXTRAS} ${MHM2_SOURCE}
  make -j 1 all check install
  cp -p $envir ${INSTALL_PREFIX}-cpu-${t}/env.sh # store environment to support execution
done
)

echo "Done building $(date) in ${SECONDS}"

