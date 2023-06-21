#!/bin/bash

# CI Accuracy script for perlmutter

set -e
set -x

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
mkdir -p ${RUN_PREFIX}

cd ${MHM2_SOURCE}
git describe
cd upcxx-utils
git describe
cd -

export GASNET_BACKTRACE=1

cd $CI_SCRATCH
dt=$(date '+%Y%m%d_%H%M%S')

export ARCTIC=${SCRATCH}/GitlabCIData/
for f in arctic_sample_0.fq  arctic_samples.fq  arcticsynth-refs.fa
do
  if [ ! -f ${ARCTIC}/$f ]
  then
    echo "Missing $f in $ARCTIC"
    exit 1
  fi
done

export GASNET_BACKTRACE=1
slurm_jobs=
echo "Testing builds on perlmutter"

export GASNET_USE_HUGEPAGES=0

for arch in gpu cpu
do
  slurm_opts="--job-name='CIm${arch}a-${CI_COMMIT_SHORT_SHA} -C $arch --qos=debug --time=12:00 --account=m2865"
  nodes=8
  if [ "$arch" == "gpu" ] ; then slurm_opts="${slurm_opts}_g -G $((4*nodes))" ; fi
  inst=${INSTALL_PREFIX}-${arch}
  DBG=$inst-Debug/bin/mhm2.py
  REL=$inst-Release/bin/mhm2.py
  RWDI=$inst-RelWithDebInfo/bin/mhm2.py
  OPTS="-r ${ARCTIC}/arctic_samples.fq --checkpoint=no"
  echo "Submitting job on ${nodes} $arch nodes"
  old=${SLURM_MEM_PER_CPU}
  unset SLURM_MEM_PER_CPU
  job=$(sbatch --parsable ${slurm_opts} --nodes=$nodes --time=8:00 --wrap="set -x ; source ${inst}-Release/env.sh; module list; env|grep SLURM; env|grep UPC; env|grep FI; set -x ; ${REL} $OPTS -o ${RUN_PREFIX}/$arch-rel && echo Good")
  export SLURM_MEM_PER_CPU=${old}
  echo "${arch} JOB ${job}"
  slurm_jobs="$slurm_jobs $job"
done

for job in $slurm_jobs
do
  echo "Waiting for jobs to complete: $slurm_jobs at $(date)"
  while /bin/true
  do
    sleep 60 
    echo "Checking for ${job} at $(date)"
    sacct=$(sacct -j $job -o state -X -n 2>/dev/null || true)
    if [ -n "${sacct}" -a -z "$(echo "${sacct}" | grep ING)" ] ; then break ; fi
    squeue -u ${USER} -o "%.16i %.10q %.5P %.10j %.8u %.5T %.10M %.9l %.16S %.9p %.6D %.6C %R"
  done
  echo "sacct $sacct"
  sacct=$(sacct -j $job -X -n)
  echo "sacct $sacct"
  cat slurm-${job}.out
  wasgood=$(echo "${sacct}" | grep -v '0:0' || true)
  if [ -z "$wasgood" ] ; then  true ; else  echo "job ${job} failed somehow - ${wasgood}"; false ; fi
done

cmds="set -x ; env | grep SLURM ;"
waits="/bin/true"
for arch in gpu cpu
do
  for r in ${arch}-rel
  do
    d=${RUN_PREFIX}/$r
    if [ ! -f $d/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta in $d" ; fi
  done
  cmds="$cmds ${inst}-Release/bin/check_asm_quality.py --asm-dir ${RUN_PREFIX}/${arch}-rel --expected-quals ${inst}-Release/share/good-arcticsynth.txt --refs ${ARCTIC}/arcticsynth-refs.fa &"
  waits="$waits && wait -n"
done

echo "Submitting $cmds $waits"
OUT=perlmutter.accuracy.${CI_PROJECT_NAME}-${CI_COMMIT_SHORT_SHA}-${CI_COMMIT_REF_NAME}-${CI_COMMIT_TAG}-${CI_PIPELINE_ID}.${dt}.out
sbatch --output=$OUT --wait --qos=debug --time=30:00 --account=m342 -C cpu -c 8 --wrap="$cmds $waits"
cat $OUT

if [ -z "$FAILED" ] ; then echo "OK" ; else echo "Something failed somehow - ${FAILED}" ; false ; fi

echo "Done accuracy testing $(date) in ${SECONDS}"
