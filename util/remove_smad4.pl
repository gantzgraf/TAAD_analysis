#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use Text::CSV;
use Data::Dumper;

my $idir = "$RealBin/../input_files";
#my $uk   = "$idir/UK_Most_Damaging_Data.csv";
#my $yale = "$idir/Yale_Most_Damaging_Data.csv";
if (@ARGV < 1){
    die <<EOT
Usage: $0 variant_list.csv 

EOT
    ;
}

my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: "
                                                    .Text::CSV->error_diag ();
my @keep = ();
while (my $var = shift){
    open my $FH, "<:encoding(utf8)", $var or die "Could not open $var: $!";
    my $header = $csv->getline($FH);
    my $sym = 0;
    $sym++ until $sym > $#{$header} or $header->[$sym] eq 'Symbol';
    die "Column 'Symbol' not found in $var\n" if $sym > $#{$header};
    while ( my $row = $csv->getline($FH) ) {
        next if $row->[$sym] eq 'SMAD4';
        push @keep, $row;
    }
    close $FH;
    $csv->eol ("\n");
    $csv->print (\*STDOUT, $header);
    $csv->print (\*STDOUT, $_) for @keep;
}

