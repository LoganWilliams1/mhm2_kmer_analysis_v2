#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse
import re



def fail(msg):
    print("\nError: "+msg)
    sys.exit(1)


def main():

    argparser = argparse.ArgumentParser()

    argparser.add_argument("-n", default=0, type=int)

    argparser.add_argument("-N", default=1, type=int)

    upcxx_args, mhm2_args = argparser.parse_known_args()

    if upcxx_args.n == 0: fail("Number of processes (-n) must be specified")

    input_reads = False

    for m_arg in mhm2_args:
        if m_arg[0] == '-':
            if len(m_arg) == 2 and (m_arg[1] == 'r' or m_arg[1] == 'p' or m_arg[1] == 'u'):
                input_reads = True
                break
    if not input_reads:
        fail("FASTQ input must be specified")

    cmd = ["upcxx-run", "-n", str(upcxx_args.n), "-N", str(upcxx_args.N), "--"]

    mhm2_binary_path = os.path.split(sys.argv[0])[0] + "/mhm2"
    adapters_path = os.path.split(sys.argv[0])[0] + "/../share/all_adapters.fa"

    cmd.extend([mhm2_binary_path, *mhm2_args, "--adapter-refs", adapters_path])


    print("Executing:\n", " ".join(cmd))


    # avoid slow gasnet memory probe
    os.environ["GASNET_PHYSMEM_PROBE"] = "0"

    # proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.chdir("install/bin")

    result = subprocess.run(cmd)

    if result.returncode: fail("Run failed")

    log_dir = max(os.listdir(), key=lambda d: os.path.getmtime(d) if os.path.isdir(d) else 0)

    file = os.path.join(log_dir, "mhm2.log")

    def extract_value(pattern, line):
        match = re.search(pattern, line)
        return match.group(1) if match else None

    kmer_count = None
    sups = None
    time = None
    expected_count = 34477443

    with open(file, 'r') as f:
        for line in f:
            if 'Total kmers' in line:
                kmer_count = extract_value(r'.*Total kmers: *([0-9]+).*', line)
            elif 'Avg supermer' in line:
                sups = extract_value(r'.*max ([0-9]*).*', line)
            elif 'Analyzing kmers' in line:
                time = extract_value(r'.*Analyzing kmers *(.*)', line)

    print("\n----------------------------\n")
    print("Kmer Analysis Elapsed Time:\t", time, "\n")
    print("Unique Kmer Count:\t\t", kmer_count, "\n")


    
    

if __name__ == "__main__":
    main()