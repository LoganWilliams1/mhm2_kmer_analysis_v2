#!/bin/bash

# upc++ gdb debugging tool
export GASNET_BACKTRACE=1

if [ ! -d "install/bin" ]; then
    echo
    echo "Install the app to run the sample input"
    echo
    exit 1
else
    cd "install/bin"
fi


if [ -e "sample_output" ]; then
  rm -rf "sample_output"
fi

use_qf=FALSE
if [[ "${1,,}" =~ --use-qf=true ]]; then
    use_qf=TRUE
fi

upcxx-run -n 4 ./mhm2 -p ../../sample_1.fastq ../../sample_2.fastq \
        --adapter-refs ./../share/all_adapters.fa -o sample_output --use-qf="$use_qf"

#echo -n "running proxy..."
# run > /dev/null
#echo -e "\tcomplete"

if [ $? -ne 0 ]; then
    echo
    echo "Running proxy failed"
    echo
    exit 1
fi


file="sample_output/mhm2.log"

kmer_count=$(grep 'Total kmers' "$file" | sed -n 's/.*Total kmers: *\([0-9]\+\).*/\1/p')
sups=$(grep 'Avg supermer' "$file" | sed -n 's/.*max \([0-9]*\).*/\1/p')
time=$(grep 'Analyzing kmers' "$file" | sed -n 's/.*Analyzing kmers *\(.*\)/\1/p')


if [ -d "../../.build/src/kcount/CMakeFiles/kcount_cpu" ]; then
    expected_count=34479151
elif [ -d "../../.build/src/kcount/kcount-kokkos" ] || { [ -d "../../.build/src/kcount/kcount-gpu" ] && [ "$use_qf" = FALSE ]; }; then
    expected_count=34477443
elif [ -d "../../.build/src/kcount/kcount-gpu" ] && [ "$use_qf" = TRUE ]; then
    expected_count=~34283000
else
  echo "Unexpected error: build failed"
  exit 1
fi

echo -e "\n----------------------------\n"
echo "Kmer Analysis Elapsed Time: $time"
echo
echo "Sample Kmer Count:    $kmer_count"
echo "Expected:             $expected_count"              
echo
if [ "$kmer_count" == "$expected_count" ]; then
    echo "Sample run executed successfully"
elif [ "$use_qf" = TRUE ] && (( $kmer_count > 34281000 && $kmer_count < 34284000 )); then
    echo "Sample run executed successfully"
    
else
    echo "Sample run failed"
fi
echo