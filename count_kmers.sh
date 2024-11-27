#!/bin/bash

usage() {
  echo "Usage: $0 [mhm2 options] -k <kmer length>"
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

command="./mhm2.py "$args" -s 0 --use-kmer-filter=false"


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

line=$(grep 'Total kmers' "$file")
parent_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
line=$(grep 'Avg supermer' "$file")
parent_sups=$(echo "$line" | sed -n 's/.*max \([0-9]*\).*/\1/p')
line=$(grep 'Analyzing kmers' "$file")
parent_time="$(echo "$line" | sed -n 's/.*Analyzing kmers *\([0-9.]*\).*/\1/p') s"

# run proxy
cd "$proxy_path"

command="./mhm2.py "$args" "

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

line=$(grep 'Total kmers' "$file")
proxy_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
line=$(grep 'Avg supermer' "$file")
proxy_sups=$(echo "$line" | sed -n 's/.*max \([0-9]*\).*/\1/p')
line=$(grep 'Analyzing kmers' "$file")
proxy_time=$(echo "$line" | sed -n 's/.*Analyzing kmers *\(.*\)/\1/p')


echo -e "\n----------------------------\n"
echo "MHM2 Time: $parent_time"
echo "Proxy Time: $proxy_time"
echo
echo "MHM2 max supermer inserts: $parent_sups"
echo "Proxy max supermer inserts: $proxy_sups"
echo
echo "MHM2 kmers: $parent_count"
echo "Proxy kmers: $proxy_count"
