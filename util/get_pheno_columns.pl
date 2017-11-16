#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use Text::CSV;
use Data::Dumper;

my $idir = "$RealBin/../input_files";
my $uk   = "$idir/UK_Phenotype_Data.csv";
my $yale = "$idir/Yale_Phenotype_Data.csv";
die "Usage: $0 sample_list.txt\n" if @ARGV != 1;
my $samples = shift;
open (my $SAMP, $samples) or die "Could not open $samples: $!\n";
my $csv = Text::CSV->new ( { binary => 1 } ) or die "Cannot use CSV: "
                                                    .Text::CSV->error_diag ();
my %uk_pheno;
my %yale_pheno;
open my $UK, "<:encoding(utf8)", $uk or die "Could not open $uk: $!";
my $header = $csv->getline($UK);
while ( my $row = $csv->getline($UK) ) {
    $row->[0] =~ s/\s//g;
    if ($row->[0] =~ /^(\d{2})\w+(\d{4})/){
        my $id = "UK_$1_$2";
        if (exists $uk_pheno{$id}){
            warn "Duplicate ID '$id' ($row->[0])\n";
            next;
        }
        for (my $i = 0; $i < @$header; $i++){
            my $col = $header->[$i];
            $uk_pheno{$id}->{$col} = $row->[$i];
        }
    }else{
        warn "Could not parse UK ID '$row->[0]'\n";
    }
}
close $UK;
open my $YALE, "<:encoding(utf8)", $yale or die "Could not open $yale: $!";
while ( my $row = $csv->getline($YALE) ) {
    $row->[0] =~ s/\s//g;
    if ($row->[0] =~ /^(\w+_)?(\d+_\d+)/){
        my $id = "Y_$2";
        if (exists $yale_pheno{$id}){
            warn "Duplicate ID '$id' ($row->[0])\n";
            next;
        }
        for (my $i = 0; $i < @$header; $i++){
            my $col = $header->[$i];
            my $value = (defined $row->[$i] and $row->[$i] =~ /\S/) ? $row->[$i] : "N/A";
            $yale_pheno{$id}->{$col} = $value;
        }
    }else{
        warn "Could not parse Yale ID '$row->[0]'\n";
    }
}
close $YALE;

my @pheno_fields = ("Primary Diagnosis  indication for surgery: Aneurysm / Dissection / Rupture / Transection / IMH / PAU",
                    "known syndrome - Marfan / LDS / EDS",
                    "Fhx- known or probable",
                    "Age at Diagnosis (any aortic disease)", 
                    "Location of Primary Diagnosis  Ascending, Arch, Descending, Thoracoabdominal, Infrarenal",
                    "Maximal Aortic Size (cm) ",
                    );
print join("\t", (  "Primary Diagnosis", 
                    "Clinical Diagnosis",
                    "Family History",
                    "Age at Diagnosis",
                    "Primary Anatomical Presentation",
                    "Maximum Aortic Diameter (cm)",
)) . "\n";
while(my $l = <$SAMP>){
    my ($id) = split(/\s+/, $l);
    #check sample type and lookup phenotype hash
    my $hash_ref;
    if ($id =~ /^UK/){
        $hash_ref = \%uk_pheno;
    }elsif ($id =~ /^Y/){
        $hash_ref = \%yale_pheno;
    }else{
        warn "Could not determine cohort for sample $id\n";
        next;
    }
    my @row = ($id);
    if (not exists $hash_ref->{$id}){
        warn "Could not find phenotype data for sample '$id'\n";
        foreach my $p (@pheno_fields){
            push @row, "NOT_FOUND";
        }
    }else{
        foreach my $p (@pheno_fields){
            push @row, $hash_ref->{$id}->{$p};
        }
    }
    print join("\t", @row) . "\n";
}
close $SAMP;
