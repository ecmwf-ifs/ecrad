#!/usr/bin/perl

#
# drhook_merge_walltime.pl
#
# For merging wall clock time results from different MPI-tasks
# i.e. DR_HOOK_OPT=prof 
#
# Original script by Eckhard Tschirschnitz, Cray, 2006 (Mflop/s)
#
# Usage: cat drhook.* | perl -w drhook_merge_walltime.pl
#

use strict;
use warnings;

# this expects concatenated dr_hook listings (wall clock time listings)

my $bignum = 999999999;

my $skip = 1;
my $str = "THRESHOLD_PERCENT";
my $threshold_percent = exists($ENV{$str}) ? $ENV{$str} : 0.1;
my $tottim = 0;

my $maxwall = 0;
my $minwall = $bignum;
my $avgwall = 0;
my $stdevwall = 0;

my $nproc = 0; # no of MPI-tasks
my $maxomp = 0; # max. no of OpenMP-threads encountered
my $exe = ""; # the name of the executable

my %nself = ();
my %sumself = ();
my %sum2self = ();
my %maxself = ();
my %minself = ();
my %ompself = ();

my %numcalls = ();

my %taskhits = ();
my %omphits = ();

RECORD: for (<>) {
  chomp; # get rid of the newline character

  next RECORD if (m/^\s*$/); # a blank line

  if (m/^\s*Profiling\b/) { 
    if ($nproc == 0) {
      $exe = $1 if (m/program='([^\'].*)'/);
    }
    $nproc++;
    $skip = 1;
  }
  elsif (m/^\s+Wall-time\s+is\s+(\S+)\s+/) {
    my $value = $1;
    $maxwall = $value if ($value > $maxwall);
    $minwall = $value if ($value < $minwall);
    $avgwall += $value;
    $stdevwall += $value * $value;
    $skip = 1;
  }
  elsif (m/^\s+1\s+/) { # the first record of (any) file encountered
    &CumulativeUpdate();  
    %ompself = ();
    $skip = 0;
  }

  next RECORD if ($skip);

  #         rank  %time cumul self    total  #_of_calls self:ms/call tot:ms/call routine_name
  if (m/^\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+\S+\s+(.*)\s*/) {
    my $self = $1;
    my $ncalls = $2;
    my $name = $3;
    my $tid = 0;
    $name =~ s/\s+//g;
    $name =~ s/^[*]//;
    #print "$self $name\n";
    if ($name =~ m/^(.*)[@](\d+)/) {
      $tid = $2;
      $maxomp = $tid if ($tid > $maxomp);
      $name = $1;
    }
    unless (exists($sumself{$name})) {
      $nself{$name} = 0;
      $sumself{$name} = 0;
      $sum2self{$name} = 0;
      $maxself{$name} = 0;
      $minself{$name} = $bignum;
      $numcalls{$name} = 0;
      $taskhits{$name} = 0;
      $omphits{$name} = 0;
    }
    $numcalls{$name} += $ncalls;
    $taskhits{$name} += 1 if ($tid == 1);
    $omphits{$name} += 1;
    
    # Account the most expensive OpenMP thread only per $name label
    $ompself{$name} = 0 unless (exists($ompself{$name}));
    $ompself{$name} = $self if ($self > $ompself{$name});
  }
}

if ($nproc > 0) {
  # One final time ...
  &CumulativeUpdate();  
    
  print  STDOUT "The name of the executable : $exe\n";
  print  STDOUT "Number of MPI-tasks        : $nproc\n";
  print  STDOUT "Number of OpenMP-threads   : $maxomp\n";
  print  STDOUT "export $str   : $threshold_percent%\n";
  $avgwall /= $nproc;
  if ($nproc > 1) {
    $stdevwall = ($stdevwall - $nproc * $avgwall * $avgwall)/($nproc - 1);
    $stdevwall = ($stdevwall > 0) ? sqrt($stdevwall) : 0; # be careful with rounding of errors
  }
  else {
    $stdevwall = 0;
  }
  #printf STDOUT ("Total time (average)       : %.3f secs \n", $avgwall);
  printf STDOUT ("Wall-times over %d MPI-tasks (secs) : Min=%.3f, Max=%.3f, Avg=%.3f, StDev=%.3f\n", 
		 $nproc, $minwall, $maxwall, $avgwall, $stdevwall);
  my $avgpercent = $threshold_percent * $avgwall / 100.0;

  printf STDOUT ("Routines whose average time > %.2f%% (%.3f secs) of the total average will be included in the listing\n",
		 $threshold_percent,$avgpercent);

  printf STDOUT ("%7s %10s %10s %10s %8s %8s %12s %13s %8s %8s : %s\n",
		 "Avg-%", "Avg.secs", 
		 "Min.secs", "Max.secs", "St.dev", 
		 "Imbal-%", 
		 "# of calls",
		 "Avg.msec/call",
		 "Tasks",
		 "OpenMP",
		 "Name-of-the-routine");

  open(PIPE,"|sort -nr|sed 's/ %/%/'");

  # Values accounted for
  my $acc_avgpercent = 0;
  my $acc_avgtime = 0;
  my $acc_maxtime = 0;
  my $acc_mintime = $bignum;

  foreach my $name (keys %sumself) {
    my $value = $sumself{$name};
    my $navg = $nself{$name};
    my $avgself = $value/$navg;
    if ($avgself > $avgpercent) {
      my $percent = $value/$tottim*100;
      my $stdev = 0;
      if ($navg > 1) {
	$stdev = $sum2self{$name};
	$stdev = ($stdev - $navg * $avgself * $avgself)/($navg - 1);
	$stdev = ($stdev > 0) ? sqrt($stdev) : 0; # be careful with rounding of errors
      }
      printf PIPE ("%6.2f %% %10.3f %10.3f %10.3f %8.3f %7.2f%% %12.0f %13.3f %8d %8d : %s\n",
		   $percent, $avgself,
		   $minself{$name},$maxself{$name},$stdev,
		   ($maxself{$name} - $minself{$name})/$maxself{$name}*100,
		   $numcalls{$name},
		   $avgself / $numcalls{$name} * 1000.0,
		   $taskhits{$name}, 
		   $omphits{$name} / $taskhits{$name},
		   $name);
      # Update values accounted for
      $acc_avgpercent += $percent;
      $acc_avgtime += $avgself;
      $acc_mintime = $minself{$name} if ($acc_mintime > $minself{$name});
      $acc_maxtime = $maxself{$name} if ($acc_maxtime < $maxself{$name});
    }
  }
  close(PIPE);
  printf STDOUT ("%6.2f%% %10.3f %10.3f %10.3f\n",
		 $acc_avgpercent,$acc_avgtime,
		 $acc_mintime, $acc_maxtime);
}

exit 0;

sub CumulativeUpdate {
    foreach my $name (keys %ompself) {
	my $self = $ompself{$name};
	$nself{$name} += 1;
	$sumself{$name} += $self;
	$sum2self{$name} += $self * $self;
	$maxself{$name} = $self if ($self > $maxself{$name});
	$minself{$name} = $self if ($self < $minself{$name});
	$tottim += $self;
    }
}

# A sample of typical output
__DATA__
The name of the executable : /lustre/tmp/saarines/rundir/LUMO/19/fmiopt/teho.gnu/h1/N7+1.T140+20xt1xh1.dcmp_10x14.ppn_20+20.nproma_60.hgp_undef.fmiopt-534667/MASTERODB
Number of MPI-tasks        : 140
Number of OpenMP-threads   : 1
export THRESHOLD_PERCENT   : 0.2%
Wall-times over 140 MPI-tasks (secs) : Min=399.410, Max=400.320, Avg=399.876, StDev=0.220
Routines whose average time > 0.20% (0.800 secs) of the total average will be included in the listing
  Avg-%   Avg.secs   Min.secs   Max.secs   St.dev  Imbal-%   # of calls Avg.msec/call    Tasks   OpenMP : Name-of-the-routine
  5.32%     21.270     19.341     24.105    0.920   19.76%     10598700         0.002      140        1 : LAITRI
  5.18%     20.730     16.166     23.143    1.239   30.15%       698740         0.030      140        1 : VDFHGHTNHL
  5.06%     20.215      0.469     70.180   14.105   99.33%         8540         2.367      140        1 : SLCOMM
  4.12%     16.474     12.309     17.213    0.952   28.49%       933656         0.018      140        1 : FPCINCAPE
  3.61%     14.429     11.903     16.617    0.774   28.37%       698740         0.021      140        1 : INITAPLPAR
  3.02%     12.055      0.671     45.288   10.498   98.52%         6860         1.757      140        1 : SLCOMM2A
  2.72%     10.875      5.275     16.432    2.938   67.90%       698740         0.016      140        1 : RAIN_ICE
  2.70%     10.793      4.965     18.614    3.433   73.33%       698740         0.015      140        1 : RAIN_ICE:RAIN_ICE_SEDIMENTATION_STAT
  2.55%     10.182      9.763     10.552    0.212    7.48%        25340         0.402      140        1 : TRMTOL
  2.54%     10.169      8.510     10.509    0.444   19.02%       698740         0.015      140        1 : APL_AROME
  2.51%     10.054      8.519     10.192    0.421   16.41%      1397480         0.007      140        1 : COMPUTE_FUNCTION_THERMO_MF
  2.51%     10.038      8.181     11.849    1.156   30.96%        25340         0.396      140        1 : TRLTOG
  2.31%      9.252      5.219     21.083    4.795   75.25%        21000         0.441      140        1 : TRLTOM
  1.93%      7.711      1.976     24.881    6.385   92.06%        22540         0.342      140        1 : TRGTOL
  1.91%      7.653      3.985     12.186    2.001   67.30%       698740         0.011      140        1 : RAIN_ICE:RAIN_ICE_SLOW
  1.84%      7.364      6.510      7.976    0.403   18.38%      1413160         0.005      140        1 : ELARCHE
  1.46%      5.825      3.269      8.145    1.265   59.86%        71300         0.082      140        1 : RRTM_RTRN1A_140GP
  1.33%      5.301      4.389      5.637    0.229   22.14%      5589920         0.001      140        1 : TRIDIAG_MASSFLUX
  1.30%      5.208      4.416      5.282    0.218   16.40%      1397480         0.004      140        1 : TURB:COMPUTE_FUNCTION_THERMO
  1.25%      4.997      4.445      5.741    0.406   22.57%      9260880         0.001      140        1 : RPASSF
  1.12%      4.489      3.801      4.562    0.189   16.68%       698740         0.006      140        1 : ACTQSAT
  0.93%      3.699      2.879      4.401    0.341   34.58%       698740         0.005      140        1 : CONDENSATION
  0.90%      3.588      0.507      9.260    1.812   94.52%       698740         0.005      140        1 : RAIN_ICE:RAIN_ICE_FAST_RG

etc.

  0.23%      0.922      0.891      0.978    0.025    8.90%       718900         0.001      140        1 : SIPTP
  0.21%      0.840      0.583      1.262    0.160   53.80%        13440         0.063      140        1 : ESEIMPLS
  0.21%      0.834      0.691      0.878    0.037   21.30%       698740         0.001      140        1 : LATTEX
  0.21%      0.825      0.698      0.845    0.035   17.40%    112297500         0.000      140        1 : SWTT1
  0.18%      0.893      0.000      1.548    0.567  100.00%       477554         0.002      113        1 : ECUME_FLUX
 85.41%    341.669      0.000     70.180
