#!/bin/bash

set -x
set -e

USAGE="$0 base_dir
Optionally set UPCXX_VER to download and install that version of UPCXX
Optionall set CI_EXTRA_PATH to include something build outside of this CI

"


BASE=$1
if [ -z "$BASE" ]
then
	echo $USAGE
	exit 1
fi
UPCXX_VER=${UPCXX_VER:=2021.3.0}
echo "Using upcxx version $UPCXX_VER"

CI_CMAKE_OPTS=${CI_CMAKE_OPTS}
echo "Using CI_CMAKE_OPTS=${CI_CMAKE_OPTS}"

CI_INSTALL=$BASE/ci-install-mhm2-upcxx-${UPCXX_VER}
PATH=$CI_INSTALL/bin:/bin:/usr/bin:/usr/local/bin:${CI_EXTRA_PATH}
export HIPMER_DATA=${BASE}/scratch/
mkdir -p ${HIPMER_DATA}
export CI_SCRATCH=${BASE}/scratch/mhm2-${CI_COMMIT_SHORT_SHA}-${CI_COMMIT_REF_NAME}-${CI_COMMIT_TAG}
export RUN_PREFIX=${CI_SCRATCH}/runs
export INSTALL_PREFIX=${CI_SCRATCH}
export GASNET_BACKTRACE=1

echo "Establishing all tests under BASE=$BASE and CI_SCRATCH=$CI_SCRATCH"
exec >  >(tee -ia ${CI_SCRATCH}/validation.log)
exec 2> >(tee -ia ${CI_SCRATCH}/validation.err >&2)
echo "Logging to ${CI_SCRATCH}/validation.log and .err at $(date) on $(uname -n) in $(pwd)"  
echo "Validating all tests under BASE=$BASE and CI_SCRATCH=$CI_SCRATCH"

export MHM2_SOURCE=$(pwd)
env
uname -a
pwd
find * -type d -ls
date

FAILED=""
which upcxx || FAILED="${FAILED} could not install upcxx"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]
upcxx --version || FAILED="${FAILED} no upcxx was found"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

env
df -h
uname -a
pwd
date
cd ${CI_SCRATCH}
reads=${HIPMER_DATA}/arctic_sample_0.fq 
export DBG=${INSTALL_PREFIX}/mhm2-dbg/bin/
export RWDI=${INSTALL_PREFIX}/mhm2-rwdi/bin/
export REL=${INSTALL_PREFIX}/mhm2-rel/bin/

echo "Starting RelWithDebInfo mhm2 on Arctic $reads using ${RWDI}/mhm2"
if [ ! -f arctic_sample_0.fq -a -f $reads ] ; then ln -s $reads ; fi
if [ ! -f arcticsynth-refs.fa -a -f ${HIPMER_DATA}/arcticsynth-refs.fa ] ; then ln -s ${HIPMER_DATA}/arcticsynth-refs.fa ; fi
${RWDI}/ci_asm_qual_test.sh || FAILED="${FAILED} Could not run ci_asm_qual_test"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

if [ ! -f ./test-arctic-sample0/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta for rwdi in test-arctic-sample0" ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

mv test-arctic-sample0 ${RUN_PREFIX}/rwdi-test-arctic-sample0
echo "Starting debug run on with reduced workflow using ${DBG}/mhm2"
${DBG}/mhm2.py -r $reads -o ${RUN_PREFIX}/dbg-k3163 --kmer-lens 31,63 --checkpoint=no -v || FAILED="${FAILED} Could not run dbg"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "verify dbg-k3163 results $(ls -la ${RUN_PREFIX}/dbg-k3163)"
if [ ! -f ${RUN_PREFIX}/dbg-k3163/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta on dbg-k3163" ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

${DBG}/check_asm_quality.py --asm-dir ${RUN_PREFIX}/dbg-k3163 --expected-quals ${DBG}/../share/good-arctic-sample0-k31k63.txt --refs ${HIPMER_DATA}/arcticsynth-refs.fa || echo "WARN did not pass check_asm_quality.py for dbg-k3163 with reduced workflow k31k63"
if [ -z "$FAILED" ] ; then  true ; else echo "Something failed somehow - ${FAILED}"; false ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Completed $0 Successfully at $(date)"
