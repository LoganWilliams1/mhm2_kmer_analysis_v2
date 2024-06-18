#!/bin/bash

set -e

USAGE="$0 base_dir
Optionally set UPCXX_VER to download and install that version of UPCXX
Optionally set CI_UPCXX_CONFIGURE_OPTS to add extra options when building upcxx
Optionally set CI_CMAKE_OPTS to change the build options (such as -DENABLE_CUDA=off)

"


BASE=$1
if [ -z "$BASE" ]
then
	echo $USAGE
	exit 1
fi
UPCXX_VER=${UPCXX_VER:=2023.9.0}
GASNET_VER=${GASNET_VER:=2024.5.0}
echo "Using upcxx version $UPCXX_VER and gasnet ${GASNET_VER}"

CI_CMAKE_OPTS=${CI_CMAKE_OPTS}
echo "Using CI_CMAKE_OPTS=${CI_CMAKE_OPTS}"

git submodule init
git submodule sync
git submodule update
git describe --always
( cd upcxx-utils ; git describe --always )

export CI_INSTALL=${CI_INSTALL:=$BASE/ci-install-${CI_PROJECT_NAME}-upcxx-${UPCXX_VER}-${GASNET_VER}}
export HIPMER_DATA=${HIPMER_DATA:=${BASE}/scratch/}
export CI_SCRATCH=${CI_SCRATCH:=${BASE}/scratch/${CI_PROJECT_NAME}-${CI_COMMIT_SHORT_SHA}-${CI_COMMIT_REF_NAME}-${CI_COMMIT_TAG}-${CI_PIPELINE_ID}}
export RUN_PREFIX=${RUN_PREFIX:=${CI_SCRATCH}/runs}
export INSTALL_PREFIX=${INSTALL_PREFIX:=${CI_SCRATCH}}

mkdir -p ${HIPMER_DATA}
export GASNET_BACKTRACE=1

echo "Verifing and/or downloading test arctic data set"
export ARCTIC_URL=https://portal.nersc.gov/project/hipmer/MetaHipMer_datasets_12_2019/ArcticSynth/samples/
[ -f ${HIPMER_DATA}/arcticsynth-refs.fa ] || ( curl -LO ${ARCTIC_URL}/arcticsynth-refs.fa.gz && gunzip arcticsynth-refs.fa.gz && mv arcticsynth-refs.fa ${HIPMER_DATA}/arcticsynth-refs.fa )
for i in 0 1 2 3 4 5 6 7 8 9 10 11 ; do [ -f ${HIPMER_DATA}/arctic_sample_$i.fq ] || ( curl -LO ${ARCTIC_URL}/arctic_sample_$i.fq.gz && gunzip arctic_sample_$i.fq.gz && mv arctic_sample_$i.fq ${HIPMER_DATA}/arctic_sample_$i.fq ) ; done

echo "Establishing all tests under BASE=$BASE and CI_SCRATCH=$CI_SCRATCH"
set -x
mkdir -p ${CI_SCRATCH}
chmod a+rx ${CI_SCRATCH}
chmod g+s ${CI_SCRATCH}
mkdir -p ${RUN_PREFIX}

export MHM2_SOURCE=$(pwd)
uname -a
uptime
pwd
find * -type d -ls -maxdepth 3 || /bin/true
date

echo "Purging any old tests"
find ${BASE}/scratch -maxdepth 1  -name "${CI_PROJECT_NAME}-*-*-*"  -mtime +7 -type d -exec rm -rf '{}' ';' || /bin/true
df -h

echo "PATH=$PATH"
FAILED=""
echo "Checking for cmake, Berkeley UPC and UPC++"
which cmake && cmake --version || FAILED="${FAILED} cmake not found"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Checking or building upcxx"
which upcxx || UPCXXVER=${UPCXX_VER} ./upcxx-utils/contrib/install_upcxx.sh $CI_INSTALL --enable-gasnet-verbose ${CI_UPCXX_CONFIGURE_OPTS} || FAILED="${FAILED} could not install upcxx"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

which upcxx
upcxx --version

upcxx --version || FAILED="${FAILED} no upcxx was found"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Building all flavors with '${CI_CMAKE_OPTS}'"
mkdir -p ${RUN_PREFIX}
rm -rf $CI_SCRATCH/build
mkdir -p $CI_SCRATCH/build
cd $CI_SCRATCH/build

echo "Building Debug"
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}/mhm2-dbg -DCMAKE_BUILD_TYPE=Debug ${CI_CMAKE_OPTS} ${MHM2_SOURCE} || FAILED="${FAILED} Could not configure Debug"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make -j 16 all || make VERBOSE=1 all || FAILED="${FAILED} Could not build Debug"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make check install || FAILED="${FAILED} Could not build Debug"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make clean
export UPCXX_CODEMODE=debug
echo "Building RelWithDebInfo with UPCXX_CODEMODE=$UPCXX_CODEMODE"
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}/mhm2-rwdi -DCMAKE_BUILD_TYPE=RelWithDebInfo ${CI_CMAKE_OPTS} ${MHM2_SOURCE} || FAILED="${FAILED} could not configure RelWithDebInfo"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make -j 16 all check install || FAILED="${FAILED} Cuold not build RelWithDebInfo"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make clean
unset UPCXX_CODEMODE
echo "Building Release"
cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}/mhm2-rel -DCMAKE_BUILD_TYPE=Release ${CI_CMAKE_OPTS} ${MHM2_SOURCE} || FAILED="${FAILED} could not configure Release"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

make -j 16 all check install || FAILED="${FAILED} Could not build Release"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

if [ -z "$FAILED" ] ; then  true ; else echo "Something failed somehow - ${FAILED}"; false ; fi

cd -
rm -rf $CI_SCRATCH/build

echo "Completed Successfully: '$0 $@' at $(date) on $(uname) for ${SECONDS} s"
