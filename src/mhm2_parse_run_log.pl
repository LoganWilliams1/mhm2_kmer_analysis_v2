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

our @metrics = qw {Date Operator DataSetName GBofFASTQ NumReads MinDepth ContigKmers DistinctKmersWithFP MinDepthKmers Assembled>5kbp Version HipmerWorkflow Nodes Threads CoresPerNode NumRestarts ManualRestarts TotalTime NodeHours CalculatedTotalTime NumStages StageTime };
our @fields = (@metrics, @modules, "RunDir", "RunOptions");
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
   print join("\t", @fields) . "\n";
   print join("\t", @stats{@fields}) . "\n";
}

$stats{"Operator"} = $ENV{"USER"};
my ($diploid, $zygo, $cgraph, $meta);
$stats{"NumRestarts"} = 0;
$stats{"ManualRestarts"} = 0;

our %h_units = ( 'B' => 1./1024./1024./1024., 'K' => 1./1024./1024., 'M' => 1./1024., 'G' => 1., 'T' => 1024.);
my $stage_pat = '\w+\.[hc]pp';

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
    if (/ kmer-lens = \s+(\d+.*)/) {
        $stats{"ContigKmers"} = $1;
        $stats{"ContigKmers"} =~ s/ *$//;
        $stats{"ContigKmers"} =~ s/ /-/g;
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

    if (/Finished in ([\d\.]+) s at (\d+)\/(\d+)\/(\d+) /) {
        $stats{"Date"} = "20" . $4 . "-" . $2 . "-" . $3;
        $stats{"CalculatedTotalTime"} = $1;
    }
    
    if (/Overall time taken (including any restarts): ([\d\.]+) s/) { # fixme is not printed in mhm2.log!
        $stats{"TotalTime"} += $1;
    }

    if (/Post processing/) {
        $post_processing = 1;
    }
    if ($post_processing && /Aligning reads to contigs \s+([\d\.]+) s/) {
        $stats{"PostAlignment"} = $1;
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

        if (/Found (\d+) .* unique kmers/ || /Number of elements in hash table: (\d+)/) {
            $stats{"DistinctKmersWithFP"} = $1;
        }
        if (defined $stats{"DistinctKmersWithFP"} && /Purged (\d+) kmers \( /) {
            $stats{'MinDepthKmers'} = $stats{"DistinctKmersWithFP"} - $1;
        }
        if (/After purge of kmers < .*, there are (\d+) unique kmers/) {
            if (defined $stats{"MinDepthKmers"}) {
                $stats{"DistinctKmersWithFP"} = $stats{"MinDepthKmers"}; # the previous one
            }
            $stats{"MinDepthKmers"} = $1; # the last one
        }
        if (not defined $stats{'MinDepthKmers'}) {
           if (/ hash table final size is (\d+) entries and final load factor/) {
            $stats{'MinDepthKmers'} = $1;
           }
        }
        if (not defined $stats{'DistinctKmersWithFP'}) {
          if (/read kmers hash table: purged \d+ .* singleton kmers out of (\d+)/) {
            $stats{'DistinctKmersWithFP'} = $1;
          }
        }
        if (/Completed contig round k =/) {
            $firstUFX = 0;
        }
    }

}
$stats{"CoresPerNode"} = $stats{"Threads"} / $stats{"Nodes"};

if (not defined $stats{"TotalTime"}) {
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
  if (not defined $stats{'tmpLoadedReads'} or not defined $stats{'tmpMergedReads'} or $stats{'tmpLoadedReads'} < $stats{'tmpMergedReads'}) { print STDERR "Wrong number of loaded vs merged reads  $stats{'tmpLoadedReads'} vs  $stats{'tmpMergedReads'}\n"; }
  $stats{'NumReads'} = $stats{'tmpMergedReads'} + $stats{'tmpLoadedReads'} - $stats{'tmpMergedReads'}
}

printStats();

print "MHM2, version " . $stats{"Version"} . ", was executed on " . $stats{"NumReads"} . " reads" . " and " 
       . $stats{"GBofFASTQ"} . " GB of fastq " . "for " . $stats{"TotalTime"} . " seconds in a job over " 
       . $stats{"Nodes"} . " nodes (" . $stats{"Threads"} . " threads) using the " . $stats{"HipmerWorkflow"} . " workflow.\n";

