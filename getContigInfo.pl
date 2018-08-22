#! /usr/bin/perl -w
use strict;
use warnings;
use Fatal;
use IO::File;
use IO::Handle;
sub usage {
print <<EOF;
Usage: $0 <contigs.fa> <minimum contig length>
EOF
exit;
}
!@ARGV && &usage();

my $minLength;
$minLength = $ARGV[1] || 0;

my %contig_length = ();
my %contig_n = ();
my $contig_id = "";

my $fh = &open_maybe_compressed($ARGV[0]);
while (my $line = $fh->getline) {
    chomp $line;
    if($line =~ m/^>(.+)/) {
        $contig_id = $1;
        $contig_length{$contig_id} = 0;
        $contig_n{$contig_id} = 0;
    } 
        else {
        $contig_length{$contig_id} += length $line;
        $contig_n{$contig_id} += $line =~ tr/[Nn]/N/;
    }
}

my $contig_count = 0;
my $contig_filtered = 0;
my $n_count = 0;
my $total_read_length = 0;

my $hundred = 0;
my $thousand = 0;
my $five_thousand = 0;
my $ten_thousand = 0;

while (my ($k,$v) = each %contig_length) {
    if ($v > $minLength) {
        $total_read_length += $v;
        $contig_count++;
        $n_count += $contig_n{$k};
        $v > 100 and $hundred++;
        $v > 1000 and $thousand++;
        $v > 5000 and $five_thousand++;
        $v > 10000 and $ten_thousand++;
    }
    else {
        $contig_filtered++;
    }
}

#my @length_order = sort {$a<=>$b} values(%contig_length);
my @keys = sort { $contig_length{$b} <=> $contig_length{$a} } keys(%contig_length);
my @length_order = @contig_length{@keys};
my $N90 = 0;
my $N90_ctg = 0;
my $N50 = 0;
my $N50_ctg = 0;
my $last_ctg;
for (my $i = 0; $i <= $#length_order; $i++) {
    $length_order[$i] < $minLength and $last_ctg = $i-1 and last;
    $N90 += $length_order[$i];
    if ($N90 > ($total_read_length * 0.9) and $N90_ctg == 0) {
        $N90_ctg = $i;
    }
    $N50 += $length_order[$i];
    if ($N50 > ($total_read_length * 0.5) and $N50_ctg == 0) {
        $N50_ctg = $i;
    }
}
$last_ctg = $last_ctg || -1;
my $n_content=sprintf("%.3f",$n_count/$total_read_length*100);

print "input fasta file: $ARGV[0]\n";
print "length cutoff: $minLength bp ($contig_filtered records are filtered)\n" if $minLength > 0;
print "\n";
print "minimum length: $length_order[$last_ctg] ($keys[$last_ctg])\n";
print "maximum length: $length_order[0] ($keys[0])\n";
print "2nd long contig: $length_order[1] ($keys[1])\n" if $#length_order > 1;
print "3rd long contig: $length_order[2] ($keys[2])\n" if $#length_order > 2;
print "total length: $total_read_length\n";
print "number of Ns: $n_count ($n_content%)\n";
#print "N content: $n_content%\n";
print "avg. length: ".int ($total_read_length/$contig_count)."\n";
print "L90: ",$N90_ctg+1," ($keys[$N90_ctg])\n";
print "N90: $length_order[$N90_ctg]\n";
print "L50: ",$N50_ctg+1," ($keys[$N50_ctg])\n";
print "N50: $length_order[$N50_ctg]\n\n";

print "contig counts:\n";
print " > $minLength bp: $contig_count\n";
print "( >   100 bp: $hundred)\n";
print "( >  1000 bp: $thousand)\n";
print "( >  5000 bp: $five_thousand)\n";
print "( > 10000 bp: $ten_thousand)\n\n";

# --------------------------------------------------------------------------------------------------
sub open_maybe_compressed {
    my($fname) = @_;
    my @try = ("pbzip2 -dc", "pigz -dc", "bzip2 -dc", "gzip -dc", "cat");
    my $fh;
    for my $exe (@try) {
        my $io = IO::File->new("$exe \Q$fname\E 2> /dev/null |");
        my $c = $io->getc;
        if (defined $c) {
            $io->ungetc(ord $c);
            #print STDERR "Using $exe for $fname\n";
            return $io;
        }
    }
    die "could not open $fname";
}
