#!/bin/bash

# CI Validation script for perlmutter

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
mkdir -p ${RUN_PREFIX}

cd ${MHM2_SOURCE}
git describe
cd upcxx-utils
git describe
cd -

export ARCTIC=${SCRATCH}/GitlabCIData/
export GASNET_BACKTRACE=1

cd $CI_SCRATCH

for f in arctic_sample_0.fq  arctic_samples.fq  arcticsynth-refs.fa
do
  if [ ! -f ${ARCTIC}/$f ]
  then
    echo "Missing $f in $ARCTIC"
    exit 1
  fi
done

slurm_jobs=
echo "Testing builds on perlmutter"

for arch in gpu cpu
do
  slurm_opts="-C $arch --qos=debug --time=30:00 --account=m2865"
  nodes=5
  if [ "$arch" == "gpu" ]
  then
    nodes=2
    slurm_opts="${slurm_opts}_g -G $((4*nodes))"
  fi
  inst=${INSTALL_PREFIX}-${arch}
  DBG=$inst-Debug/bin/mhm2.py
  REL=$inst-Release/bin/mhm2.py
  RWD=$inst-RelWithDebug/bin/mhm2.py
  RWDI=$inst-RelWithDebInfo/bin/mhm2.py
  OPTS="-r ${ARCTIC}/arctic_sample_0.fq -v --checkpoint=yes"
  echo "Submitting job on ${nodes} $arch nodes"
  old=${SLURM_MEM_PER_CPU}
  unset SLURM_MEM_PER_CPU
  job=$(sbatch --parsable --job-name="CImvg-${CI_COMMIT_SHORT_SHA}" ${slurm_opts} --nodes=${nodes} --wrap="source $inst-Release/env.sh ; module list ; env|grep SLURM; env|grep GAS; env|grep UPC; env|grep FI; ${RWDI} $OPTS -o ${RUN_PREFIX}/$arch-rwdi-0 && echo 'sleeping for 90s so that srun can recover' && sleep 90 && ${DBG} ${OPTS} --kmer-lens 63 -o ${RUN_PREFIX}/$arch-dbg-0-k63 && echo Good")
  export SLURM_MEM_PER_CPU=${old}
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

for arch in gpu cpu
do
  d1=${arch}-rwdi-0 
  d2=${arch}-dbg-0-k63
  for r in $d1 $d2
  do
    d=${RUN_PREFIX}/$r
    if [ ! -f $d/final_assembly.fasta ] ; then FAILED="${FAILED} Did not find final_assembly.fasta in $d" ; fi
  done
  $inst-Release/bin/check_asm_quality.py --asm-dir ${RUN_PREFIX}/${d1} --expected-quals ${inst}-Release/share/good-arctic-sample0.txt --refs ${ARCTIC}/arcticsynth-refs.fa || FAILED="${FAILED} Did not pass check_asm_quality.py on ${d1}"
  $inst-Release/bin/check_asm_quality.py --asm-dir ${RUN_PREFIX}/${d2} --expected-quals ${inst}-Release/share/good-arctic-sample0-k63.txt --refs ${ARCTIC}/arcticsynth-refs.fa || FAILED="${FAILED} Did not pass check_asm_quality.py on ${d2}"
done

if [ -z "$FAILED" ] ; then echo "OK" ; else echo "Something failed somehow - ${FAILED}" ; false ; fi

echo "Done validating $(date) in ${SECONDS}"
