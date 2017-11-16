#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use Text::CSV;
use Data::Dumper;
use lib "$RealBin/lib/dapPerlGenomicLib";
use VcfReader 0.3;

my $idir = "$RealBin/../input_files";
#my $uk   = "$idir/UK_Most_Damaging_Data.csv";
#my $yale = "$idir/Yale_Most_Damaging_Data.csv";
my $uk   = "$idir/UK_All_Variants_Data.csv";
my $yale = "$idir/Yale_All_Variants_Data.csv";
if (2 < @ARGV or @ARGV < 1){
    die <<EOT
Usage: $0 variant_list.txt [ExAC VCF]

Your variant list should contain two columns - one for gene name and the other
for the HGVSc notation. If an ExAC VCF is provided, allele frequencies will 
also be reported. This script is designed to work with the following ExAC VCF:
ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz

EOT
    ;
}

my $var_list = shift;
my $exac = shift;
my %search_args = ();
if ($exac){
    %search_args = VcfReader::getSearchArguments($exac);
}
my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: "
                                                    .Text::CSV->error_diag ();

# store in silico predictions in hash 
# keys are gene symbols to hashes of HGVSc ID to hashes of in silico tool to 
# prediction value
my %var_to_insilico = ();
parse_variants($uk);
parse_variants($yale);

open (my $FH, $var_list) or die "Could not read $var_list: $!\n";
print join("\t", qw/ Symbol HGVSc Polyphen SIFT Condel CADD/);
print "\tExAC_AF\tExAC_Max_AF" if %search_args;
print "\n";
while (my $line = <$FH>){
    chomp $line;
    my @split = split(/\s+/, $line);
    die "Need at least two columns - could not parse: $line\n" if @split < 2;
    my @row = @split[0..1];
    if (exists $var_to_insilico{$split[0]}->{$split[1]}){
        foreach my $f (qw/Polyphen SIFT Condel CADD/){
            push @row, $var_to_insilico{$split[0]}->{$split[1]}->{$f};
        }
        if (%search_args){
            if ($var_to_insilico{$split[0]}->{$split[1]}->{ExAC_AF}){
                push @row, $var_to_insilico{$split[0]}->{$split[1]}->{ExAC_AF};
            }else{
                push @row, 'N/A';
            }
            if ($var_to_insilico{$split[0]}->{$split[1]}->{ExAC_Max_AF}){
                push @row, $var_to_insilico{$split[0]}->{$split[1]}->{ExAC_Max_AF};
            }else{
                push @row, 'N/A';
            }
        }
    }else{
        warn "Could not find variant info for $split[0] - $split[1]\n";
        my $x = %search_args ? 6 : 4;
        push @row, ('N/A') x $x;
    }
    print join("\t", @row) . "\n";
}

#################################################
sub parse_variants{
    my $var = shift;
    open my $FH, "<:encoding(utf8)", $var or die "Could not open $var: $!";
    my $header = $csv->getline($FH);
    my @req = qw/Symbol HGVSc Polyphen SIFT Condel CADD Chrom Pos Ref Alt/;
    my %cols = ();
    foreach my $r (@req){
        my $i = 0;
        $i++ until $i > $#{$header} or $header->[$i] eq $r;
        die "Column '$r' not found in $var\n" if $i > $#{$header};
        $cols{$r} = $i;
    }
    while ( my $row = $csv->getline($FH) ) {
        my $symbol = $row->[$cols{Symbol}];
        my $hgvsc = $row->[$cols{HGVSc}];
        $hgvsc =~ s/^ENST\d+\.\d+://;
        for (my $i = 0; $i < @$header; $i++){
            my $col = $header->[$i];
            my $value = (defined $row->[$i] and $row->[$i] =~ /\S/) ? $row->[$i] : "N/A";
            $var_to_insilico{$symbol}->{$hgvsc}->{$col} = $value;
        }
        if (%search_args){
            my ($af, $popmax) = get_exac_freq(  $row->[$cols{Chrom}],
                                                $row->[$cols{Pos}],
                                                $row->[$cols{Ref}],
                                                $row->[$cols{Alt}],
                                             );
            $var_to_insilico{$symbol}->{$hgvsc}->{ExAC_AF} = $af; 
            $var_to_insilico{$symbol}->{$hgvsc}->{ExAC_Max_AF} = $popmax; 
        }
    }
    close $FH;
}

#################################################
sub get_exac_freq{
    my ($chrom, $pos, $ref, $alt) = @_;
    print STDERR "Searching ExAC VCF for $chrom:$pos-$ref/$alt...\n";
    ($pos, $ref, $alt) = VcfReader::reduceRefAlt($pos, $ref, $alt);
    my @exac_hits = VcfReader::searchForPosition(
        %search_args,
        chrom => $chrom,
        pos   => $pos,
    );
    my $af = 'N/A';
    my $popmax = 'N/A';
HIT: foreach my $exac_line(@exac_hits){
        my @exac_split = split( "\t", $exac_line );
        my %exac_min = VcfReader::minimizeAlleles( \@exac_split, );
        foreach my $i (keys(%exac_min)){
            if ($exac_min{$i}->{POS} == $pos and 
                $exac_min{$i}->{REF} eq $ref and 
                $exac_min{$i}->{ALT} eq $alt ){
                my @afs = split(",", 
                          VcfReader::getVariantInfoField(\@exac_split, "AF"));
                my @popmax_acs = split(",", 
                    VcfReader::getVariantInfoField(\@exac_split, "AC_POPMAX"));
                my @popmax_ans = split(",", 
                    VcfReader::getVariantInfoField(\@exac_split, "AN_POPMAX"));
                $af = $afs[$i-1];
                if ($popmax_ans[$i-1] and $popmax_ans[$i-1] ne 'NA'){
                    $popmax = sprintf("%.3g", $popmax_acs[$i-1]/$popmax_ans[$i-1]);
                }
                print STDERR "Found matching variant with AF = $af (POPMAX = $popmax)\n";
                last HIT;
            }
        }
    }
    return ($af, $popmax)
}
