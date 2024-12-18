#!/bin/bash

usage() {
  echo "Usage: $0 [upc++ options] ./mhm2 [mhm2 options] -k <kmer length>"
  exit 1
}

# upc++ gdb debugging tool
export GASNET_BACKTRACE=1

# build command

args="$*"
k_found=false

while [[ "$#" -gt 0 ]]; do
  if [[ "$1" == "-k" ]]; then
    if [[ "$2" =~ ^[0-9]+$ ]]; then
      k_found=true
      break
    else
      usage
    fi
  else
    shift
  fi
done

if ! $k_found; then
  usage
fi

command="upcxx-run "$args" --scaff-kmer-lens 0 --use-kmer-filter=false"


proxy_path="$(pwd)/install/bin"

# run parent (edit path)
parent_path=~/mhm2/install/GPU/bin
cd "$parent_path"

echo
#echo -n "running parent..."
eval $command #> /dev/null
#echo -e "\tcomplete"

if [ $? -ne 0 ]; then
    echo
    echo "Running MHM2 parent failed"
    echo
    exit 1
fi

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

parent_count=$(grep 'Total kmers' "$file" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
parent_sups=$(grep 'Avg supermer' "$file" | sed -n 's/.*max \([0-9]*\).*/\1/p')
parent_time="$(grep 'Analyzing kmers' "$file" | sed -n 's/.*Analyzing kmers *\([0-9.]*\).*/\1/p') s"
parent_hc_kmers=$(grep 'High count kmer' "$file" | sed -E 's/.*count = ([0-9]+) kmer = ([A-Z]+)/\2: \1/' | sort -t: -k2,2n)

# run proxy
cd "$proxy_path"

command="upcxx-run "$args" "

#echo -n "running proxy..."
eval $command #> /dev/null
#echo -e "\tcomplete"

if [ $? -ne 0 ]; then
    echo
    echo "Running proxy failed"
    echo
    exit 1
fi

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

proxy_count=$(grep 'Total kmers' "$file" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
proxy_sups=$(grep 'Avg supermer' "$file" | sed -n 's/.*max \([0-9]*\).*/\1/p')
proxy_time=$(grep 'Analyzing kmers' "$file" | sed -n 's/.*Analyzing kmers *\(.*\)/\1/p')
proxy_hc_kmers=$(grep 'High count kmer' "$file" | sed -E 's/.*count = ([0-9]+) kmer = ([A-Z]+)/\2: \1/' | sort -t: -k2,2n)



echo -e "\n----------------------------\n"
echo "MHM2 Time: $parent_time"
echo "Proxy Time: $proxy_time"
echo
echo "MHM2 kmers: $parent_count"
echo "Proxy kmers: $proxy_count"
echo
echo "MHM2 high count kmers:"
echo "$parent_hc_kmers"
echo
echo "Proxy high count kmers:"
echo "$proxy_hc_kmers"