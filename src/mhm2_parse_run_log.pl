#!/usr/bin/env perl

use strict;
use warnings;

use File::stat;
use File::Basename;

our %stats;
# translated modules from hipmer's parse_run_log.pl as best as possible with same number of columns
our @modules = qw {ReadFastq MergeReads kcount TraverseDeBruijn klign LocalAssm klign-kernel NA2 NA3 Localize Checkpoint WriteOutputs PostAlignment TotalMinusIO NA4 CGraph};

# allow many names to map to a standard module
our %module_map;
foreach my $m (@modules) {
    $module_map{$m} = $m;
}
$module_map{'UFX'} = 'kcount';
$module_map{'LoadFQ'} = 'ReadFastq';
$module_map{'Meraculous'} = 'TraverseDeBruijn';
$module_map{'merAligner'} = 'klign';
$module_map{'FindNewMers'} = 'Checkpoint';
$module_map{'PrefixMerDepth'} = 'WriteOutputs';
$module_map{'ProgressiveRelativeDepth'} = 'PostAlignment';
$module_map{'ContigMerDepth'} = 'TotalMinusIO';
$module_map{'oNo'} = 'klign-kernel';
$module_map{'ParCC'} = 'NA2';
$module_map{'GapClosing'} = 'NA3';
$module_map{'ContigEndAnalyzer'} = 'NA4';

our @metrics = qw {Date Operator DataSetName GBofFASTQ NumReads MinDepth ContigKmers DistinctKmersWithFP 
                   MinDepthKmers Assembled>5kbp Version HipmerWorkflow Nodes Threads CoresPerNode NumRestarts 
                   ManualRestarts TotalTime NodeHours CalculatedTotalTime NumStages StageTime };
our @fields = (@metrics, @modules, "RunDir", "RunOptions", "NumRawReads", "RawGbp", "NumMergedReads", 
                "NumRawPairs", "MergedGbp", "TotalAssembled", "InitialFreeMemGB", "PeakMemGB", 
                "InitialFreeGPUGB", "PeakGPUGB", "PeakKmerCountLoadFactor");
foreach my $module (@modules) {
    $stats{$module} = 0;
}
$stats{'CachedIO'} = 'Assembled>5kbp';
$stats{'NumStages'} = 0;
$stats{'StageTime'} = 0;
$stats{'tmpMergedReads'} = 0; # 2021-08-06 15:06:19   merged 14255494 (15.44%) pairs
$stats{'tmpLoadedReads'} = 0; # Loaded 126665730 tot_bases=18044017886 names=0B
$stats{'tmpLoadedBases'} = 0; # Loaded 126665730 tot_bases=18044017886 names=0B



sub printStats {
   foreach my $field (@fields) {
     if (not defined $stats{$field}) { print STDERR "No $field found\n"; }
   }
   print join(";", @fields) . "\n";
   print join(";", @stats{@fields}) . "\n";
}

$stats{"Operator"} = $ENV{"USER"};
my ($diploid, $zygo, $cgraph, $meta);
$stats{"NumRestarts"} = 0;
$stats{"ManualRestarts"} = 0;
$stats{"NumRawReads"} = 0;
$stats{"NumRawPairs"} = 0;
$stats{"RawGbp"} = 0;
$stats{'MergedGbp'} = 0;
$stats{'NumMergedReads'} = 0;

our %h_units = ( 
'B' => 1./1024./1024./1024., ' B' => 1./1024./1024./1024., 
'K' => 1./1024./1024., 'KB' => 1./1024./1024., 
'M' => 1./1024., 'MB' => 1./1024.,
'G' => 1., 'GB' => 1.,
'T' => 1024., 'TB' => 1024.);
my $stage_pat = '\w+\.[hc]pp';

my $completed = 0;
my $num_gpus = 0;
my $gpu_peak_is_more_accurate = 0;
my $since_read_kmers = 0;
my $post_processing = 0;
my $knowGB = 0;
my $firstUFX = 1;
$stats{"MinDepth"} = 2.0;
my $restarted = 0;
while (<>) {
    s/[\000-\037]\[(\d|;)+m//g; # remove any control characters from the log
    if (/Total size of \d+ input .* is (\d+\.\d\d)(.)/) {
      my $unit = $h_units{$2};
      $stats{"GBofFASTQ"} = ($1+0.0) * $unit;
    }

    if (/Executed as: (.+)/) {
        if (not defined $stats{"RunOptions"}) {
            $stats{"RunOptions"} = $1;
        }
        if (/--restart/) {
            $stats{"NumRestarts"}++;
            $restarted = 1;
        }
    }
    if (/MHM2 version (\S+) with upcxx-utils /) {
        $stats{"Version"} = $1;
    }
    if (!$restarted && /Starting run with (\d+) processes on (\d+) node.? at/) {
        $stats{"Threads"} = $1;
        $stats{"Nodes"} = $2;
    }
    if (/#  config_file:\s+(\S+)/) {
        #fixme
        $stats{"Config"} = $1;
    }
    if (/ kmer-lens = \s+\[(\d+.*?)\]/ || / kmer-lens = \s+(\d+.*)/) {
        $stats{"ContigKmers"} = $1;
        $stats{"ContigKmers"} =~ s/ *$//;
        $stats{"ContigKmers"} =~ s/[ ,]/-/g;
    }
    if (/output = \s+(\S+)/) {
        $stats{"RunDir"} = $1;
        if ( -d $stats{"RunDir"} ) { # get the user too
           my $uid = stat($stats{"RunDir"})->uid;
           $stats{"Operator"} = getpwuid($uid);
        }
    }
    if (/ scaff-kmer-lens =\s+(\S+)/) {
        $cgraph = ($1 ne "" and $1 ne "0") ? 1 : 0;
    }

    if (/ min-depth-thres =\s+ (\d+\.?\d*)$/) {
        $stats{"MinDepth"} = $1;
    }
    if (/ checkpoint = \s+(\S+)/) {
        $stats{"Checkpoint"} = $1;
    }
    
    if (/ ${stage_pat}:Merge reads: ([\d\.]+) s/) {
        $stats{"MergeReads"} = $1;
    }
    if (/ ${stage_pat}:Analyze kmers: ([\d\.]+) s/) {
        $stats{"kcount"} = $1;
    }
    if (/ ${stage_pat}:Traverse deBruijn graph: ([\d\.]+) s/) {
        $stats{"TraverseDeBruijn"} = $1;
    }
    if (/ ${stage_pat}:Alignments: ([\d\.]+) s/) {
        $stats{"klign"} = $1;
    }
    if (/ ${stage_pat}:Kernel alignments: ([\d\.]+) s/) {
        $stats{"klign-kernel"} = $1;
    }
    if (/ ${stage_pat}:Local assembly: ([\d\.]+) s/) {
        $stats{"LocalAssm"} = $1;
    }
    if (/ ${stage_pat}:Traverse contig graph: ([\d\.]+) s/) {
        $stats{"CGraph"} = $1;
    }
    if (/ FASTQ total read time: ([\d\.]+)/) {
        $stats{'ReadFastq'} += $1;
    }
    if (/Loading reads into cache \s+([\d\.]+) s/) {
        $stats{'ReadFastq'} += $1;
    }
    if (/ merged FASTQ write time: ([\d\.]+)/ || / Contigs write time: ([\d\.]+)/) {
        $stats{'WriteOutputs'} += $1;
    }
    if (/ > 5kbp: \s+(\d+) /) {
        $stats{'Assembled>5kbp'} = $1;
    }
    if (/Total assembled length:\s+(\d+)/) {
        $stats{'TotalAssembled'} = $1;
    }

    if (/Finished in ([\d\.]+) s at (\d+)\/(\d+)\/(\d+) /) {
        $stats{"Date"} = "20" . $4 . "-" . $2 . "-" . $3;
        $stats{"CalculatedTotalTime"} = $1;
    }
    if (/Initial free memory across all .*, ([\d\.]+)(.?B) max/) {
        $stats{"InitialFreeMemGB"} = $1 * $h_units{$2};
    }
    if (/Peak memory used across all .*, ([\d\.]+)(.?B) max/) {
        $stats{"PeakMemGB"} = $1 * $h_units{$2};
    }
    if (/read kmers in hash/ || /stats for read kmers pass/) {
        $since_read_kmers = 0;
    }
    if ($since_read_kmers++ < 4 && / load factor.* avg,? ([\d\.]+) max/) {
        $stats{"PeakKmerCountLoadFactor"} = $1;
    }

    if (/Available GPU memory per rank for kmers hash table is .* accounting for PnP of ([\d\.]+)(.?B)/) {
        $gpu_peak_is_more_accurate = $1 * $h_units{$2};
    }
    if (/GPU read kmers hash table used ([\d\.]+)(.?B) memory on GPU out of ([\d\.]+)(.?B)/) {
        my $gb = $1 * $h_units{$2};
        my $max = $3 * $h_units{$4};
        if ($gpu_peak_is_more_accurate > 0) {
            $gb += $gpu_peak_is_more_accurate;
        }
        
        if ((not defined $stats{'PeakGPUGB'}) || $stats{'PeakGPUGB'} < $gb) {
            $stats{'PeakGPUGB'} = $gb
        }
        
        if ((not defined $stats{'InitialFreeGPUGB'}) || $stats{'InitialFreeGPUGB'} < $max) {
            $stats{'InitialFreeGPUGB'} = $max;
        }
    }
    
    if (/ Finished in ([\d\.]+) s/) {
        $stats{"TotalTime"} += $1;
    }

    if (/Post processing/) {
        $post_processing = 1;
    }
    if ($post_processing && /Aligning reads to contigs \s+([\d\.]+) s/) {
        $stats{"PostAlignment"} = $1;
    }

    if (/Loaded .*: reads=(\d+) tot_bases=(\d+) names=/) {
        $stats{'NumMergedReads'} += $1;
        $stats{'MergedGbp'} += $2 / 1000000000;
    }
    if (/Finished reading .* tot_num_reads=(\d+) tot_num_pairs=(\d+) tot_num_bases=(\d+)/) {
        $stats{"NumRawReads"} += $1;
        $stats{"NumRawPairs"} += $2;
        $stats{"RawGbp"} += $3 / 1000000000;
    }

    if ($firstUFX) {
        if (/Processed a total of (\d+) reads/) {
            $stats{"NumReads"} = $1;
        }
        if (/Loaded (\d+) tot_bases=(\d+) names=/) {
            $stats{'tmpLoadedReads'} += $1;
            $stats{'tmpLoadedBases'} += $2;
        }
        if (/merged (\d+) .* pairs/) {
            $stats{'tmpMergedReads'} += $1;
        }
        if (/[pP]urged (\d+) .* singleton kmers out of (\d+)/) { # kcount_cpu & kcount gpu
            $stats{'Purged'} = $1;
            $stats{'DistinctKmersWithFP'} = $2;
            $stats{'MinDepthKmers'} = $2 - $1;
        }
        if (/For (\d+) kmers, average kmer count /) {
            if ((not defined $stats{'MinDepthKmers'}) || $1 > $stats{'MinDepthKmers'}) { # fix for pre Issue190
                $stats{'MinDepthKmers'} = $1;
            }
        }
        if ((not defined $stats{"DistinctKmersWithFP"}) && (/Found (\d+) .* unique kmers/ || /Number of elements in hash table: (\d+)/)) { # legacy
            $stats{"DistinctKmersWithFP"} = $1;
        }

        if (not defined $stats{"MinDepthKmers"}) { # legacy
            if (/After purge of kmers < .*, there are (\d+) unique kmers/) {
                if (defined $stats{"MinDepthKmers"} && (not defined $stats{"DistinctKmersWithFP"})) {
                    $stats{"DistinctKmersWithFP"} = $stats{"MinDepthKmers"}; # the previous one
                }
                $stats{"MinDepthKmers"} = $1; # the last one
            }
            if (not defined $stats{'MinDepthKmers'}) {
                if (/ hash table final size is (\d+) entries and final load factor/) {
                    $stats{'MinDepthKmers'} = $1;
                }
            }
        }
        if (defined $stats{"DistinctKmersWithFP"} && defined $stats{'Purged'}) {
            if (not defined $stats{'MinDepthKmers'}) {
                $stats{'MinDepthKmers'} = $stats{"DistinctKmersWithFP"} - $stats{'Purged'};
            }
        }
        if (not defined $stats{'DistinctKmersWithFP'}) { # legacy
          if (/: purged \d+ .* singleton kmers out of (\d+)/) {
            $stats{'DistinctKmersWithFP'} = $1;
          }
        }
        if (/Completed contig round k =/) {
            $firstUFX = 0;
        }
    }
    if ((not defined $stats{'MinDepthKmers'}) && (defined $stats{'DistinctKmersWithFP'} && defined $stats{'Purged'})) {
        $stats{'MinDepthKmers'} = $stats{'DistinctKmersWithFP'} - $stats{'Purged'};
    }
    if (/Available number of GPUs on this node (\d+)/) {
        $num_gpus = $1;
    }
    if (/Total time before close and finalize/) {
        $completed = 1;
    }
}
$stats{"CoresPerNode"} = $stats{"Threads"} / $stats{"Nodes"};

if ((not defined $stats{"TotalTime"}) || $stats{"CalculatedTotalTime"} > $stats{"TotalTime"}) {
    $stats{"TotalTime"} = $stats{"CalculatedTotalTime"};
}

$stats{"TotalMinusIO"} = $stats{"TotalTime"} - $stats{"ReadFastq"} - $stats{"WriteOutputs"};

$stats{"TotalTime"} =~ s/\..*//;
$stats{"NodeHours"} = $stats{"TotalTime"} * $stats{"Nodes"} / 3600.0;
$stats{"NodeHours"} =~ s/(\.\d)\d*/$1/;

$stats{"HipmerWorkflow"} = "Normal";
if ($post_processing) {
    $stats{"HipmerWorkflow"} .= " with Post-Alignment";
}

$stats{"DataSetName"} = $stats{"RunDir"};
if ($stats{"DataSetName"} =~ /\//) {
    $stats{"DataSetName"} = basename($stats{"DataSetName"});
}
foreach my $module (@modules) {
    $stats{$module} =~ s/\..*//;
}

$stats{"GBofFASTQ"} =~ s/(\.\d)\d*/$1/;
if (not defined $stats{'NumReads'}) {
  if ((not defined $stats{'tmpLoadedReads'}) or (not defined $stats{'tmpMergedReads'}) or $stats{'tmpLoadedReads'} < $stats{'tmpMergedReads'}) { print STDERR "Wrong number of loaded vs merged reads  $stats{'tmpLoadedReads'} vs  $stats{'tmpMergedReads'}\n"; }
  $stats{'NumReads'} = $stats{'tmpMergedReads'} + $stats{'tmpLoadedReads'} - $stats{'tmpMergedReads'}
}

if ($num_gpus == 0) {
    $stats{'InitialFreeGPUGB'} = '';
    $stats{'PeakGPUGB'} = '';
}

printStats();

if ($completed) {
    print "MHM2, version " . $stats{"Version"} . ", was executed on " . $stats{"NumReads"} . " reads" . " and " 
        . $stats{"GBofFASTQ"} . " GB of fastq " . "for " . $stats{"TotalTime"} . " seconds in a job over " 
        . $stats{"Nodes"} . " nodes (" . $stats{"Threads"} . " threads and " . $num_gpus . " gpus) using the " 
        . $stats{"HipmerWorkflow"} . " workflow.\n";
    if (defined $stats{'PeakKmerCountLoadFactor'} && $stats{'PeakKmerCountLoadFactor'} < 0.667) {
        my $mem_utility =($stats{'PeakKmerCountLoadFactor'} / 0.667);
        if ($num_gpus == 0 && defined $stats{'PeakMemGB'} && defined $stats{'InitialFreeMemGB'}
            && $stats{'PeakMemGB'} < $stats{'InitialFreeMemGB'}) {
            $mem_utility *= ($stats{'PeakMemGB'} / $stats{'InitialFreeMemGB'});
        } elsif ($num_gpus > 0 && $gpu_peak_is_more_accurate && defined $stats{'PeakGPUGB'} && defined $stats{'InitialFreeGPUGB'}) {
            $mem_utility *= ($stats{'PeakGPUGB'} / $stats{'InitialFreeGPUGB'});
        }
        my $tgt_nodes = ($stats{'Nodes'} * (1.+2.*$mem_utility)/3.); # 2/3 of the way to mem_utility target
        if ($tgt_nodes < $stats{'Nodes'}) {
            print "This " . $stats{'Nodes'} . " node job had enough memory and might run more efficiently on fewer nodes (Mem utility = " 
                . int($mem_utility * 100. + 0.5) . "\%). Suggestion: " . int($tgt_nodes*10)/10. 
                . " nodes for max " . int($stats{TotalTime} / $mem_utility+0.5) . " s.\n";
        }
    }
} else {
    print "MHM2 failed after running for " . $stats{"TotalTime"} . " seconds\n";
}

