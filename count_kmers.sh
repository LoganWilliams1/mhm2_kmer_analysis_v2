#!/bin/bash

usage() {
  echo "Usage: $0 [mhm2 options] -k <kmer length>"
  exit 1
}

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

command="./mhm2.py "$args" -s 0"


proxy_path="$(pwd)/install/bin"

# run parent (edit path)
parent_path=~/mhm2/install/bin
cd "$parent_path"

echo
#echo -n "running parent..."
eval $command #> /dev/null
#echo -e "\tcomplete"

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

line=$(grep 'Total kmers' "$file")
parent_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
line=$(grep 'Analyzing kmers' "$file")
parent_time=$(echo "$line" | sed -n 's/.*Analyzing kmers *\(.*\)/\1/p')

# run proxy
cd "$proxy_path"

#echo -n "running proxy..."
eval $command #> /dev/null
#echo -e "\tcomplete"

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

line=$(grep 'Total kmers' "$file")
proxy_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
line=$(grep 'Analyzing kmers' "$file")
proxy_time=$(echo "$line" | sed -n 's/.*Analyzing kmers *\(.*\)/\1/p')


echo -e "\n----------------------------\n"
echo "MHM2 Time: $parent_time"
echo "Proxy Time: $proxy_time"
echo
echo "MHM2 kmers: $parent_count"
echo "Proxy kmers: $proxy_count"

