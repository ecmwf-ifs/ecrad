#!/usr/bin/perl
#!/usr/local/apps/perl/current/bin/perl
#
# (C) Copyright 2010- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Modified to cope with rttov files which have built in interface
# blocks that just need to be parsed out.
# Paul Burton 2010

use strict;
use warnings;
#use lib "/home/rd/rdx/bin/prepifs/perl";
use Fortran90_stuff;
use Data::Dumper;
$Data::Dumper::Indent = 1;

my $rttov_intf=0;

{
  my (@files);
  (@files) = @ARGV;
  &setup_parse();
  my $locintfbldir=$ENV{LOC_INTFBDIR} or die "LOC_INTFBDIR not defined ";
  my $intfbldir=$ENV{INTFBDIR} or die "INTFBDIR not defined ";
  our $study_called;
  FILE:for (@files) {
    my (@interface_block);
    my (@line_hash);
    chomp;
 # Read in lines from file
    my $fname = $_;

    my $int_block_fname=$fname;
    $int_block_fname=~s/\.F90/.intfb.h/;
    $int_block_fname=~s#.*/(.+)$#$1#;
    my $ofname=$intfbldir.'/'.$int_block_fname;
    my $nfname=$locintfbldir.'/'.$int_block_fname;
    
# Do nothing if file hasn't changed since intfb already created

    if ( (-f $nfname) &&
         ( (stat($nfname))[9] > (stat($fname))[9] )) {
      print "INTERFACE BLOCK $int_block_fname UP TO DATE\n";
      next FILE;
    }

# skip .h files, Nils
    my $base = $_;
    ( $base = $fname ) =~ s/\.(\w+)\s*$//;
    my $suffix = $1;
    next if ( $suffix eq "h" );
# end Nils
    my @statements=();
    my %prog_info=();
    my @lines = &readfile($fname);

    $rttov_intf=0;
    if ($fname=~/^satrad\//) {
      if (grep(/^!INTF_END/,@lines)) {
        print "Working on rttov file $fname\n";
        $rttov_intf=1;
        &create_rttov_interface_block(\@lines,\@interface_block,\%prog_info);
      } else {
        # satrad file without INTF_END marker - ignore
        print "Ignoring satrad file $fname (no INTF_END marker)\n";
        next;
      }
    } else {
      print "Working on file $fname\n";
      &expcont(\@lines,\@statements);
      $study_called=0;
      &study(\@statements,\%prog_info);
#     print Dumper(\%prog_info);
    }
    unless($prog_info{is_module}) {
      if ($rttov_intf) {
        @lines=@interface_block;
      } else {
        &create_interface_block(\@statements,\@interface_block);
        &cont_lines(\@interface_block,\@lines,\@line_hash);
      }
      if ( -f $nfname ) {
	my @nldlines=&readfile($nfname);
	if(&eq_array(\@nldlines, \@lines)){
          print "INTERFACE BLOCK $int_block_fname NOT UPDATED \n" ;
          next FILE;
        }
      }
      if ( -f $ofname ) {
	my @oldlines=&readfile($ofname);
	if(&eq_array(\@oldlines, \@lines)){
          print "INTERFACE BLOCK $int_block_fname UNCHANGED \n" ;
          next FILE;
        }
      }
      print "WRITE INTERFACE BLOCK $int_block_fname \n";
      print "$nfname \n";
      &writefile($nfname,\@lines);
    }
  }
}
sub eq_array {
    my ($ra, $rb) = @_;
    return 0 unless $#$ra == $#$rb;
    for my $i (0..$#$ra) {
	return 0 unless $ra->[$i] eq $rb->[$i];
    }
    return 1;
}

sub create_rttov_interface_block {
  my ($lines,$intfblk,$prog_info)=(@_);
  my ($line,$on,$what,@intf);
  $on=1;
  $what="";
  for $line (@{$lines}) {
    $on=1 if ($line=~/^\s*!\s*INTF_ON\s*$/);
    $on=0 if ($line=~/^\s*!\s*INTF_OFF\s*$/);
    last if ($line=~/^\s*!INTF_END\s*$/);
    if ($what eq "") {
      $what="SUBROUTINE" if ($line=~/^\s*SUBROUTINE/i);
      $what="FUNCTION" if ($line=~/\s*FUNCTION/i);
      $what="MODULE" if ($line=~/^\s*MODULE/i);
      $what="PROGRAM" if ($line=~/^\s*PROGRAM/i);
    }
    next if ($line=~/^\s*!/ || $line=~/^\s*$/ || !$on );
    push(@intf,$line);
  }
  die "Cannot guess if file contains a SUBROUTINE or FUNCTION.\n" unless $what;
  $prog_info->{is_module}=1 if ($what eq "MODULE");
  @{$intfblk}=("INTERFACE\n",@intf,"END $what\n","END INTERFACE\n");
}
    
