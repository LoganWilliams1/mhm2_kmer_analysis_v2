#!/bin/bash

proxy_path=$(pwd)

# run parent (edit path)
parent_path=~/mhm2/install/bin

cd "$parent_path"

./mhm2.py "$@"

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

line=$(grep 'Total kmers' "$file")

parent_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')

# run proxy
cd "$proxy_path/install/bin"

./mhm2.py "$@"

log_dir=$(ls -t | head -n 1)

file="$log_dir/mhm2.log"

line=$(grep 'Total kmers' "$file")

proxy_count=$(echo "$line" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')

echo "----------------------------------------"

echo

echo "Parent kmers: $parent_count"
echo "Proxy kmers: $proxy_count"
