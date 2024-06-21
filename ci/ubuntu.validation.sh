#!/bin/bash

set -e

USAGE="$0 base_dir
Optionally set UPCXX_VER to download and install that version of UPCXX

"


BASE=$1
if [ -z "$BASE" ]
then
	echo $USAGE
	exit 1
fi
UPCXX_VER=${UPCXX_VER:=2023.9.0}
echo "Using upcxx version $UPCXX_VER"

CI_CMAKE_OPTS=${CI_CMAKE_OPTS}
echo "Using CI_CMAKE_OPTS=${CI_CMAKE_OPTS}"

export CI_INSTALL=${CI_INSTALL:=$BASE/ci-install-${CI_PROJECT_NAME}-upcxx-${UPCXX_VER}}
export HIPMER_DATA=${HIPMER_DATA:=${BASE}/scratch/}
export CI_SCRATCH=${CI_SCRATCH:=${BASE}/scratch/${CI_PROJECT_NAME}-${CI_COMMIT_SHORT_SHA}-${CI_COMMIT_REF_NAME}-${CI_COMMIT_TAG}-${CI_PIPELINE_ID}}
export RUN_PREFIX=${RUN_PREFIX:=${CI_SCRATCH}/runs}
export INSTALL_PREFIX=${INSTALL_PREFIX:=${CI_SCRATCH}}

mkdir -p ${HIPMER_DATA}
export GASNET_BACKTRACE=1

echo "Validating all tests under BASE=$BASE and CI_SCRATCH=$CI_SCRATCH"

export MHM2_SOURCE=$(pwd)
df -h
uname -a
uptime
pwd
date

FAILED=""
which upcxx || FAILED="${FAILED} could not install upcxx"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]
upcxx --version || FAILED="${FAILED} no upcxx was found"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

set -x

cd ${CI_SCRATCH}
reads=${HIPMER_DATA}/arctic_sample_0.fq 
export DBG=${INSTALL_PREFIX}/mhm2-dbg/bin/
export RWDI=${INSTALL_PREFIX}/mhm2-rwdi/bin/
export REL=${INSTALL_PREFIX}/mhm2-rel/bin/

if [ ! -f arctic_sample_0.fq -a -f $reads ] ; then ln -s $reads ; fi
if [ ! -f arcticsynth-refs.fa -a -f ${HIPMER_DATA}/arcticsynth-refs.fa ] ; then ln -s ${HIPMER_DATA}/arcticsynth-refs.fa ; fi

echo "Starting Debug tiny CI using ${DBG}/mhm2"
${DBG}/ci_asm_qual_test-tiny.sh || FAILED="${FAILED} Could not run ci_asm_qual_test-tiny"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Starting RelWithDebInfo mhm2 on Arctic $reads using ${RWDI}/mhm2"
${RWDI}/ci_asm_qual_test.sh || FAILED="${FAILED} Could not run ci_asm_qual_test"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

if [ ! -f ./test-arctic-sample0/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta for rwdi in test-arctic-sample0" ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

mv test-arctic-sample0 ${RUN_PREFIX}/rwdi-test-arctic-sample0
echo "Starting debug run on with reduced workflow using ${DBG}/mhm2"
timeout -k 1m -s INT --foreground -v 25m ${DBG}/mhm2.py -r $reads -o ${RUN_PREFIX}/dbg-k3163 --kmer-lens 31,63 --scaff-kmer-lens 63 --checkpoint=no -v || FAILED="${FAILED} Could not run dbg"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "verify dbg-k3163 results $(ls -la ${RUN_PREFIX}/dbg-k3163)"
if [ ! -f ${RUN_PREFIX}/dbg-k3163/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta on dbg-k3163" ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

timeout -k 1m -s INT --foreground -v 10m ${DBG}/check_asm_quality.py --asm-dir ${RUN_PREFIX}/dbg-k3163 --expected-quals ${DBG}/../share/good-arctic-sample0-k31k63.txt --refs ${HIPMER_DATA}/arcticsynth-refs.fa || FAILED="${FAILED} did not pass check_asm_quality.py for dbg-k3163 with reduced workflow k31k63"
if [ -z "$FAILED" ] ; then  true ; else echo "Something failed somehow - ${FAILED}"; false ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Completed $0 Successfully at $(date)"
