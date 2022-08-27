#!/bin/bash

set -x
set -e

USAGE="$0 base_dir
Optionally set UPCXX_VER to download and install that version of UPCXX
Optionally set CI_EXTRA_PATH to include something built outside of this CI

"


BASE=$1
if [ -z "$BASE" ]
then
	echo $USAGE
	exit 1
fi
UPCXX_VER=${UPCXX_VER:=2021.3.0}
echo "Using upcxx version $UPCXX_VER"

CI_INSTALL=$BASE/ci-install-mhm2-upcxx-${UPCXX_VER}
PATH=$CI_INSTALL/bin:/bin:/usr/bin:/usr/local/bin:${CI_EXTRA_PATH}
export HIPMER_DATA=${BASE}/scratch/
export CI_SCRATCH=${BASE}/scratch/mhm2-${CI_COMMIT_SHORT_SHA}-${CI_COMMIT_REF_NAME}-${CI_COMMIT_TAG}
export RUN_PREFIX=${CI_SCRATCH}/runs
export GASNET_BACKTRACE=1
export INSTALL_PREFIX=${CI_SCRATCH}
export GASNET_BACKTRACE=1

exec >  >(tee -ia ${CI_SCRATCH}/accuracy.log)
exec 2> >(tee -ia ${CI_SCRATCH}/accuracy.err >&2)
echo "Logging to ${CI_SCRATCH}/accuracy.log and .err at $(date) on $(uname -n) in $(pwd)"  
echo "Running accuracy tests under BASE=$BASE and CI_SCRATCH=$CI_SCRATCH"

env
df -h
uname -a
pwd
date
upcxx --version
cd ${CI_SCRATCH}
FAILED=""

export REL=${INSTALL_PREFIX}/mhm2-rel/bin/
echo "Starting Release mhm2 on Arctic"
if [ ! -f arctic_sample_0.fq -a -f ${HIPMER_DATA}/arctic_sample_0.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_0.fq ; fi
if [ ! -f arctic_sample_1.fq -a -f ${HIPMER_DATA}/arctic_sample_1.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_1.fq ; fi
if [ ! -f arctic_sample_2.fq -a -f ${HIPMER_DATA}/arctic_sample_2.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_2.fq ; fi
if [ ! -f arctic_sample_3.fq -a -f ${HIPMER_DATA}/arctic_sample_3.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_3.fq ; fi
if [ ! -f arctic_sample_4.fq -a -f ${HIPMER_DATA}/arctic_sample_4.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_4.fq ; fi
if [ ! -f arctic_sample_5.fq -a -f ${HIPMER_DATA}/arctic_sample_5.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_5.fq ; fi
if [ ! -f arctic_sample_6.fq -a -f ${HIPMER_DATA}/arctic_sample_6.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_6.fq ; fi
if [ ! -f arctic_sample_7.fq -a -f ${HIPMER_DATA}/arctic_sample_7.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_7.fq ; fi
if [ ! -f arctic_sample_8.fq -a -f ${HIPMER_DATA}/arctic_sample_8.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_8.fq ; fi
if [ ! -f arctic_sample_9.fq -a -f ${HIPMER_DATA}/arctic_sample_9.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_9.fq ; fi
if [ ! -f arctic_sample_10.fq -a -f ${HIPMER_DATA}/arctic_sample_10.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_10.fq ; fi
if [ ! -f arctic_sample_11.fq -a -f ${HIPMER_DATA}/arctic_sample_11.fq ] ; then ln -s ${HIPMER_DATA}/arctic_sample_11.fq ; fi
if [ ! -f arcticsynth-refs.fa -a -f ${HIPMER_DATA}/arcticsynth-refs.fa ] ; then ln -s ${HIPMER_DATA}/arcticsynth-refs.fa ; fi

${REL}/ci_asm_qual_test.sh || FAILED="${FAILED} Could not run ci_asm_qual_test"
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

if [ ! -f ./test-arctic-sample0/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta on rel" ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]
mv test-arctic-sample0 ${RUN_PREFIX}/rel-test-arctic-sample0

#${REL}/ci_asm_qual_test-full.sh || FAILED="${FAILED} Could not run ci_asm_qual_test-full"
#echo "FAILED=${FAILED}" && [ -z "$FAILED" ]
#if [ ! -f ./test-arctic-samples/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta on rel" ; fi
#echo "FAILED=${FAILED}" && [ -z "$FAILED" ]
#mv test-arctic-samples ${RUN_PREFIX}/rel-test-arctic-samples

if [ -z "$FAILED" ] ; then  true ; else echo "Something failed somehow - ${FAILED}"; false ; fi
echo "FAILED=${FAILED}" && [ -z "$FAILED" ]

echo "Completed $0 Successfully at $(date)"

