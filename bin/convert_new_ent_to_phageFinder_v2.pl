#!/usr/bin/env perl

# perl script to convert the Wu GenbankGenomeRetriever .ent file to phage_finder_info.txt

# .ent file: feat_name\tend5 end3 GC%\tannotation

# phage_finder_info.txt: asmbl_id\tlength_of_assembly\tfeat_name\tend5\tend3\tannotation
# modified 06/06/2019 to populate .ent file data into a hash so that the data can be sorted by asmbl_id and end5

#   Usage: convert_ent_to_phageFinder.pl <.ent file>

use TIGR::FASTAreader;
use TIGR::Foundation;
use TIGR::FASTArecord;

die <<ENDARGS
 Bad argument list.  Call like this:
    $0 <.ent file>
ENDARGS
  unless ($#ARGV == 0);

unless (open(ENT,$ARGV[0]))  {
  die ("can't open file $ARGV[0].\n");
}
my $file = "$ARGV[0]";
$file =~ s/\.ent.*$//;
my $fasta_file = "$file.con";
my %ContigHash = ();
my %EntHash = ();
my $key = ""; #hash key = asmbl_id
my $yek = ""; #hash key = end5

$tf = new TIGR::Foundation;

if (!defined $tf){
  die ("Bad foundation\n");
}

$fr = new TIGR::FASTAreader($tf, \@errors, $fasta_file);

if (!defined $fr){
  die ("Bad reader\n");
}

while ($fr->hasNext){
  $rec = $fr->next();
  $id = $rec->getIdentifier();
  $ContigHash{$id}->{'size'} = length($rec->getData());
}

open (OUTFILE, ">phage_finder_info.txt") || die ("can't create file phage_finder_info.txt\n");
while (<ENT>)  {
  chomp;
  next if (/^\s*$/); # skip blank lines in .ent file
  @a = split(/\t+/);
  if ($a[0] =~ /:(\w.*)-/)  { # need asmbl_id for draft genomes, but not in the locus name
    $asmbl_id = $1;
    $locus = $a[0];
    $locus =~ s/:$asmbl_id//;
  }
  elsif ($a[0] =~ /-(\w+)/) {
    $asmbl_id = $1;
    $locus = $a[0];
  }

  if ($a[4] > 0) {
    $end5 = $a[2];
    $end3 = $a[3];
  }
  else {
    $end5 = $a[3];
    $end3 = $a[2];
  }

  $a[5] =~ s/similar to.*$//;
  $a[5] =~ s/identified by.*//;
  $a[5] =~ s/^ *//;

  $EntHash{$asmbl_id}{$end5}->{'end3'} = $end3;
  $EntHash{$asmbl_id}{$end5}->{'annotation'} = $a[5] . " (" . $a[1] . ")";
  $EntHash{$asmbl_id}{$end5}->{'locus'} = $locus;
  $EntHash{$asmbl_id}{$end5}->{'length'} = $ContigHash{$asmbl_id}->{'size'};
}

foreach $key (sort {$a cmp $b} keys %EntHash)  { # sort by asmbl_id
  foreach my $yek (sort {$a <=> $b} keys %{$EntHash{$key}}) { # sort by end5
    if ($EntHash{$key}{$yek}->{'locus'} !~ /RNA/) {
      print OUTFILE "$key\t$EntHash{$key}{$yek}->{'length'}\t$EntHash{$key}{$yek}->{'locus'}\t$yek\t$EntHash{$key}{$yek}->{'end3'}\t$EntHash{$key}{$yek}->{'annotation'}\n";
    }
  }
}
close (ENT);
close (OUTFILE);
exit(0);
