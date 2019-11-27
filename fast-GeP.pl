#!/usr/bin/perl
#Copy (C) 2018-2019  Massey University.
#Written by Ji Zhang, MD PhD

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY. See the GNU General Public License for 
#more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Revision notes:
#version: 1.0.1


# Switch:
#   -h: Show help
#   -v: Do not show progress on the screen.
#   -l: Produce extended results.
#   -n: Do not produce pairwise comparison files for all the isolates. 
#   -b: Using DIAMOND instead of BLAST+ as aligner in the primary search.
#   -V: Print version.

# Option:
#   -g: Name of the text file listing the genomic sequences need to be analyzed. Each file name should occupy a line.
#   -r: Name of the reference genome sequence (GenBank format).
#   -c: Minimum coverage of alignment to define an allele. [Default = 100].
#   -i: Minimum identity percentage to define an allele. [Default = 80].
#   -t: Minimum coverage of alignment to define a truncated allele. [Default = 50].
#   -d: Number of N filled in concatinating the assemblies. [Default = 200].
#   -f: Number of threads. [Default = 4].

use strict;
use Getopt::Std;
use Benchmark; 
my $T0 = time();
###################
getopts('g:r:d:c:i:t:e:f:v:V b a h l o m n ');
our ($opt_g, $opt_r, $opt_d, $opt_c, $opt_i, $opt_t, $opt_e, $opt_f, $opt_v, $opt_b, $opt_a, $opt_h, $opt_l, $opt_o, $opt_n, $opt_V);
my $distance_d = $opt_d || "200";
my $coverage_c = $opt_c || "100";
my $identity_i = $opt_i || "80";
my $coverage_t = $opt_t || "50";
my $threads_f = $opt_f || "4";
my $evalue_e = $opt_e || "0.001";
###################
my @date = readpipe("date");
if($opt_g && $opt_r){
	print "Fast-GeP analysis initiated at @date" if ($opt_v == 0);
}else{
	$opt_h = 1;
}
###################
if($opt_V == 1){
	print "fast-GeP 1.0.1\n";
	exit;
}
###################
if($opt_h == 1 && $opt_V == 0){
print ("
Fast Genome Profiler (Fast-GeP) - extraction of allele profiles from genomic sequences with genome-by-genome approach.
version: 1.0.1
date: 27.11.2019

Example command: 
fast-GeP.pl -g list_genomes.txt -r reference.gbk

or,
perl fast-GeP.pl -g list_genomes.txt -r reference.gbk

Switch:
   -h: Show help
   -v: Do not show progress on the screen.
   -l: Produce extended results.
   -n: Do not produce pairwise comparison files for all the isolates. 
   -b: Using DIAMOND instead of BLAST+ as aligner in the primary search.
   -V: Print version.

Option:
   -g: Name of the text file listing the genomic sequences need to be analyzed. Each file name should occupy a line.
   -r: Name of the reference genome sequence (GenBank format).
   -c: Minimum coverage of alignment to define an allele. [Default = 100].
   -i: Minimum identity percentage to define an allele. [Default = 80].
   -t: Minimum coverage of alignment to define a truncated allele. [Default = 50].
   -d: Number of N filled in concatinating the assemblies. [Default = 200]. 
   
");
exit;
}
###################
open(GENOME, "<$opt_g") or die "Cannot open genome list file!";
open(REPORT, ">>output.report.tmp");
open(OUT, ">genome_list.tmp");
open(OUT2, ">genome_list2.tmp");
open(OUT3, ">genome_list3.tmp");
print OUT2 "\n";
print OUT3 "genome\n";
my $counter_isolates = 0;
my(@isolate_names);
while(<GENOME>){
	if($_ =~ m/^\s/){
		next;
	}else{
		chomp;
		$counter_isolates++;
		my $isolate_name = $_;
		$isolate_name =~ s/\..*//g;
		push @isolate_names, $isolate_name;
		print OUT "$_\n";
		print OUT2 "$_\n";
		print OUT3 "$_\n";
	}
}
print REPORT "Fast-GeP analysis started at @date.<br> $opt_r was used as reference genome.<br>";
print REPORT "Allele sequences were searched with DIAMOND (identity threshold >= $identity_i).<br>" if ($opt_b == 1);
print "Allele sequences were searched with DIAMOND (identity threshold >= $identity_i).\n" if ($opt_b == 1 && $opt_v == 0);
print REPORT "Allele sequences were searched with BLAST+ (identity threshold >= $identity_i).<br>" if ($opt_b == 0);
print "Allele sequences were searched with BLAST+ (identity threshold >= $identity_i).\n" if ($opt_b == 0 && $opt_v == 0);
close GENOME;
close OUT;
###################
my ($genome, $seq, $counter, $spacer);
open(GENOME, "<genome_list.tmp");
while(<GENOME>){
	chomp;
	$genome = $_;
	open(SEQIN, "<$genome") or die "Cannot open genome sequence $genome!";
	print("concatinating $genome ...\n") if ($opt_v == 0);
	open(OUT, ">>$genome.combined.fas");
	$/=">";
	print OUT ">", "$genome\n";
	while(<SEQIN>){
		$spacer = ();
		$counter = 0;
		while($counter<$distance_d){
			$spacer = "$spacer"."N";
			$counter++;
		}
		if($_ =~ /^>/){
			next;
		}else{
			my $seq = $_;
			$seq =~ s/.*\n//;
			$seq =~ s/>//;
			$seq =~ s/\s//g;
			$seq =~ s/N/n/g;
			print OUT "$seq", "\n$spacer\n";
		}
	}
	close OUT;
	close SEQIN;
	$/="\n";
}
close GENOME;
################### parsing the reference GenBank file
print("parsing the reference GKB file $opt_r...\n") if ($opt_v == 0);
my $GBK_name = $opt_r;
open (IN, "<$GBK_name") or die "Cannot open the GBK file!";
$/="\n//";
my @infile = <IN>;
$/="\n";
close IN;
open (INF, ">>$GBK_name.xls.tmp");
open (FAA, ">>$GBK_name.faa.tmp");
open (DIC, ">>output.gene_list.txt");
my $gene_sn = 0;
print INF "aa_id\tcontig_id\taa_length\tdistance\tstart\tstop\tlocus_tag\tproduct\tnt_seq\taa_seq\n";
my(@extraction_seqs, %extraction_heads_hash);
my $aa_sn = 0;
my $count_GBKs = 0;
foreach(@infile){
	$count_GBKs++;
	$_ =~ s/\n//gs;
	$_ =~ s/  +gene  +/\n$&/gs;
	$_ =~ s/  +CDS  +/\n$&/gs;
	$_ =~ s/  +rRNA  +/\n$&/gs;
	$_ =~ s/  +tRNA  +/\n$&/gs;
	$_ =~ s/ORIGIN  +1/\n$&/;
	open (GBKs, ">$GBK_name.$count_GBKs.tmp");
	print GBKs "$_";
	close GBKs;
	my (@strands, @starts, @stops, @locus_tags, @products, @aa_seqs, @gene_types);
	my ($contig_id, $poz, $strand, $start, $stop, $gene_type, $locus_tag, $product, $nt_total);
	my ($count_genes, $nt_len, $nt_seq, $aa_seq, $aa_len, $distance, $stop_prev);
	open(GBKs, "<$GBK_name.$count_GBKs.tmp");
	while(<GBKs>){
		chomp;
		$contig_id = "$GBK_name.$count_GBKs";
		if($_ =~ /^  +CDS  +/){
			chomp;
			$aa_seq = $_;
			if($aa_seq =~ m/  +\/translation\=\"/){
				$aa_seq =~ s/.*  +\/translation\=\"//;
				$aa_seq =~ s/\".*//;
				$aa_seq =~ s/ //g;
				$poz = $_;
				$poz =~ s/\>//g;
				$poz =~ s/\<//g;
				if($poz =~ m/  +complement\([0-9]+\.\.[0-9]+\)/){
					$strand = "-";
				}else{
					$strand = "+";
				}
				$poz =~ /[0-9]+\.\.[0-9]+/;
				$poz = $&;
				$start = $poz;
				$stop = $poz;
				$stop =~ s/[0-9]+\.\.//;
				$start =~ s/\.\.[0-9]+//;
				$gene_type = "CDS";
				if($_ =~ m/  \/locus_tag\=\"/){
					$locus_tag = $_;
					$locus_tag =~ s/.*  \/locus_tag\=\"//;
					$locus_tag =~ s/\".*//;
				}elsif($_ =~ m/[0-9]+\.peg\.[0-9]+/){
					$locus_tag = $_;
					$locus_tag =~ m/[0-9]+\.peg\.[0-9]+/;
					$locus_tag = $&;
					$locus_tag =~ s/.*peg/peg/;
				}else{
					$locus_tag = "NA";
				}
				$product = $_;
				$product =~ s/.*  +\/product\=\"//;
				$product =~ s/\".*//;
				$product =~ s/ +/ /g;
				push @strands, $strand;
				push @gene_types, $gene_type;
				push @starts, $start;
				push @stops, $stop;
				push @locus_tags, $locus_tag;
				push @aa_seqs, $aa_seq;
				push @products, $product;
			}else{
				next;
			}
		}elsif($_ =~ /^ORIGIN  +1/){
			$nt_total = $_;
			$nt_total =~ s/$&//;
			$nt_total =~ s/[0-9]+//gs;
			$nt_total =~ s/ //gs;
			$nt_total =~ s/\/\///;
		}else{
			next;
		}
	}
	close GBKs;
	system("rm -f $GBK_name.$count_GBKs.tmp");
	$count_genes = -1;
	foreach(@starts){
		$count_genes++;
		$start = $_;
		$stop = $stops["$count_genes"];
		$gene_type = $gene_types["$count_genes"];
		$locus_tag = $locus_tags["$count_genes"];
		$product = $products["$count_genes"];
		$aa_seq = $aa_seqs["$count_genes"];
		$aa_len = length($aa_seq);
		$strand = $strands["$count_genes"];
		$start = $start-1;
		$nt_len = $stop - $start;
		$nt_seq = substr $nt_total, $start, $nt_len;
		if($strand eq "-"){
		$nt_seq = reverse($nt_seq);
		$nt_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
		$start = $start+1;
		}else{
			$start = $start+1;
		}
		if($gene_type eq "CDS"){
			$aa_sn++;
			if($count_genes == 0){
				$distance = $start - 1;
				$stop_prev = $stop;
				# $distance = 999999999999;
			}else{
				$distance = $start - $stop_prev - 1;
				$stop_prev = $stop;
			}
			print FAA ">gene$aa_sn", "_", "$aa_len\n$aa_seq\n";
			print INF "gene$aa_sn\t$contig_id\t$aa_len\t$distance\t$start\t$stop\t$locus_tag\t$product\t$nt_seq\t$aa_seq\n";
			print DIC "gene$aa_sn\t$product; $locus_tag\n";
			$gene_sn++;
			open(NT, ">extracted.gene$gene_sn.fas");
			print NT ">ref_gene$gene_sn\n$nt_seq\n";
			close NT;
			my $ffn_head = "gene$aa_sn"."_"."$nt_len";
			$extraction_heads_hash{$nt_seq} = $ffn_head;
			push @extraction_seqs, $nt_seq;
		}else{
			next;
		}
	}
}
close IN;
close INF;
close FAA;
my $aa_total = $aa_sn;
###################
system("makeblastdb -in $opt_r.faa.tmp -dbtype prot -logfile makeblastdb.log") if ($opt_b == 0);
system("diamond makedb --in $opt_r.faa.tmp -d db --quiet") if ($opt_b == 1);
open(GENOME, "<genome_list.tmp");
open(OUT, ">>output.txt");
open(MATRIX, ">>pa_matrix.tmp");
print MATRIX "genome\t";
print OUT "genome\taaid\tstrand\ts_len\taa_len\tcoverage\tidentity\tq_start\tq_end\ts_start\ts_end\n";
my $sn = 0;
my (@aa_presences);
while($sn<$aa_total){
	$aa_presences["$sn"] = "0";#0 code for no valid hit
	$sn++;
	print MATRIX "gene$sn\t";
}
print MATRIX "\n";
while(<GENOME>){
	chomp;
	$genome = $_;
	print("primary search for $genome ...\n") if ($opt_v == 0);
	my(@blastx, @presences);
	my($qs_prev, $qe_prev, $coverage_prev, $aaid_prev, $strand_prev, $identity_prev, $key, $copy);
	@presences = @aa_presences;
	@blastx = readpipe("diamond blastx -d db -q $genome.combined.fas --outfmt 6 --quiet --evalue $evalue_e --query-gencode 11 --max-target-seqs 0 --max-hsps 9 --more-sensitive --threads $threads_f") if ($opt_b == 1);
	@blastx = readpipe("blastx -query $genome.combined.fas -db $opt_r.faa.tmp -seg no -num_threads $threads_f -evalue $evalue_e -max_target_seqs 9999999 -query_gencode 11 -max_hsps 9 -outfmt 6") if ($opt_b == 0);
	print MATRIX "$genome";
	foreach(@blastx){
		my($aaid, $aa_len, $identity, $q_start, $q_end, $s_start, $s_end, $strand, $temp);
		my($s_len, $q_len, $coverage);
		chomp;
		my @m8_line = split(/\t/,$_);
		$aaid = $m8_line[1];
		$aaid =~ s/_.*//g;
		$aa_len = $m8_line[1];
		$aa_len =~ s/.*_//g;
		$identity = $m8_line[2];
		$q_start = $m8_line[6];
		$q_end = $m8_line[7];
		$s_start = $m8_line[8];
		$s_end = $m8_line[9];
		if($q_end>$q_start){
			$strand = "PLUS";
		}else{
			$strand = "MINUS";
			$temp = $q_end;
			$q_end = $q_start;
			$q_start = $temp;
		}
		$s_len = abs($s_end-$s_start)+1;
		$q_len = abs($q_end-$q_start)+1;
		$coverage = $s_len/$aa_len;
		$coverage = $coverage*100;
		if($identity>=$identity_i && $coverage>=$coverage_t){
			print OUT "$genome\t$aaid\t$strand\t$s_len\t$aa_len\t$coverage\t$identity\t$q_start\t$q_end\t$s_start\ts_end\n";
			if($aaid eq $aaid_prev){
				$copy++;
				$key = $aaid;
				$key =~ s/gene//;
				$key = $key - 1;
				$presences["$key"] = "2"; #2 code for multiple copy	
			}else{
				if($coverage_prev >= $coverage_c && $copy == 1){
					my($seqin, $extract_start, $extract_len, $extraction);
					open(SEQIN, "<$genome.combined.fas");
					open(ETC, ">>extracted.$aaid_prev.fas");
					$/ = undef;
					$seqin = <SEQIN>;
					$seqin =~ s/>.*\n//;
					$seqin =~ s/\s//g;
					$/ = "\n";
					if($strand_prev eq "MINUS"){
						$extract_start = $qs_prev - 4;
						$extract_len = $qe_prev - $qs_prev + 4;
						$extraction = substr $seqin, $extract_start, $extract_len;
						$extraction = reverse($extraction);
						$extraction =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
					}else{
						$extract_start = $qs_prev - 1;
						$extract_len = $qe_prev - $qs_prev + 4;
						$extraction = substr $seqin, $extract_start, $extract_len;
					}
					unless($extraction =~ m/[^atgcATGC]/){
						print ETC ">$genome.$strand_prev$aaid_prev\n$extraction\n";						
					}
					close SEQIN;
					close ETC;
					$key = $aaid_prev;
					$key =~ s/gene//;
					$key = $key - 1;
					if($extraction =~ m/[^atgcATGCN]/){
						$presences["$key"] = "9"; #9 code for containing ambiguity codes
					}elsif($extraction =~ m/N/){
						$presences["$key"] = "3"; #3 code for a truncated gene
					}else{
						$presences["$key"] = "1"; #1 code for a complete and single-copy gene
					}
					push @extraction_seqs, $extraction;
					my $extraction_head = "$aaid_prev"."_"."$extract_len";
					$extraction_heads_hash{$extraction} = $extraction_head;
				}elsif($coverage_prev < $coverage_c && $copy == 1){
					$key = $aaid_prev;
					$key =~ s/gene//;
					$key = $key - 1;
					$presences["$key"] = "3"; #3 code for a truncated gene
				}
				$qs_prev = $q_start;
				$qe_prev = $q_end;
				$aaid_prev = $aaid;
				$coverage_prev = $coverage;
				$strand_prev = $strand;
				$identity_prev = $identity;
				$copy = 1;
			}
		}else{
			next;
		}
	}
	foreach(@presences){
		chomp;
		print MATRIX "\t$_";
	}
	print MATRIX "\n";
}
close OUT;
close GENOME;
close MATRIX;
################### remove redundant allele sequences from extraction, build complementary search NT db
print "generating database for complementary search...\n" if ($opt_v == 0);
open(IN, "<pa_matrix.tmp");
open(OUT, ">>pa_matrix.transposed.tmp");
my($inline, $col_sum, $col);
my(@column);
while (<IN>){
	if($_ =~ m/^\s/){
		next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;	
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' pa_matrix.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
my %seen = ();
my @unique_allele = grep {!$seen{$_}++} @extraction_seqs;
open(IN, "<pa_matrix.transposed.tmp");
open(NTDB, ">db.ffn.tmp");
foreach(<IN>){
	unless($_ =~ m/^genome|^\s/){
		my $locus_id = $_;
		$locus_id =~ s/\t.*//g;
		chomp $locus_id;
		my $pa_string = $_;
		$pa_string =~ s/gene[0-9]+\t/\t/;
		if($pa_string =~ m/\t0/ || $pa_string =~ m/\t3/){
			foreach(@unique_allele){
				chomp;
				my $extraction_head = $extraction_heads_hash{$_};
				my $gene_id = $extraction_head;
				$gene_id =~ s/_.*//g;
				if($gene_id eq $locus_id){
					print NTDB ">$extraction_head\n$_\n";
				}
			}
		}
	}
}
close IN;
close NTDB;
system("makeblastdb -in db.ffn.tmp -dbtype nucl -logfile makeblastdb.log");
###################
open(IN, "<pa_matrix.tmp");
my %pa_string_hash = ();
while(<IN>){
	chomp;
	unless($_ =~ m/^genome|^\s/){
		my $genome = $_;
		$genome =~ s/\t.*//g;
		my $pa_string = $_;
		$pa_string =~ s/$genome\t//;
		$pa_string_hash{$genome} = $pa_string;
	}
}
close IN;
################### complementary search
open(OUT2, ">>output.blastn.txt");
print OUT2 "genome\tgnen_id\tstrand\tq_start\tq_end\tlength\textraction\n";
open(GENOME, "<genome_list.tmp");
while(<GENOME>){
	chomp;
	$genome = $_;
	print("complementary search for $genome ...\n") if ($opt_v == 0);
	my @blastn = readpipe("blastn -query $genome.combined.fas -db db.ffn.tmp -num_threads $threads_f -dust no -penalty -3 -reward 2 -word_size 11 -outfmt 6 -max_target_seqs 9999999");
	my($strand, @gene_ids, %q_start_hash, %q_end_hash, %strands_hash);
	foreach(@blastn){
		chomp;
		my @inline = split (/\t/,$_);
		my $gene_id = $inline[1];
		$gene_id =~ s/_.*//g;
		my $allele_len = $inline[1];
		$allele_len=~ s/.*_//g;
		my $identity = $inline[2];
		my $q_start = $inline[6];
		my $q_end = $inline[7];
		my $s_start = $inline[8];
		my $s_end = $inline[9];
		my $s_len = abs($s_end - $s_start) + 1;
		my $coverage = $s_len/$allele_len*100;
		if($s_end < $s_start){
			$strand = "MINUS";
		}else{
			$strand = "PLUS";
		}
		if($coverage >= $coverage_c && $identity ne 100 && $identity >= $identity_i){
			push @gene_ids, $gene_id;
			$q_start_hash{$gene_id} = $q_start;
			$q_end_hash{$gene_id} = $q_end;
			$strands_hash{$gene_id} = $strand;
		}
	}
	my %seen = ();
	my @unique_ids = grep {!$seen{$_}++} @gene_ids;
	my $pa_string = $pa_string_hash{$genome};
	my @pa = split("\t", $pa_string);
	my $sn = 0;
	my @gene_missing = ();
	foreach(@pa){
		chomp;
		$sn++;
		my $gene_id = "gene$sn";
		if($_ eq 0 || $_ eq 3){
			push @gene_missing, $gene_id;
		}
	}
	open(SEQIN, "<$genome.combined.fas");
	$/ = undef;
	my $seqin = <SEQIN>;
	$seqin =~ s/>.*\n//;
	$seqin =~ s/\s//g;
	$/ = "\n";
	close SEQIN;
	foreach(@gene_missing){
		chomp;
		my $gene_id = $_;
		my $q_start = $q_start_hash{$gene_id};
		if($q_start){
			my $q_end = $q_end_hash{$gene_id};
			my $strand = $strands_hash{$gene_id};
			open(SEQOUT, ">>extracted.$gene_id.fas");
			my($extract_start, $extract_len, $extraction);
			if($strand eq "MINUS"){
				$extract_start = $q_start - 1;
				$extract_len = abs($q_end - $q_start) + 1;
				$extraction = substr $seqin, $extract_start, $extract_len;
				$extraction = reverse($extraction);
				$extraction =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
			}else{
				$extract_start = $q_start - 1;
				$extract_len = abs($q_end - $q_start) + 1;
				$extraction = substr $seqin, $extract_start, $extract_len;
			}
			unless($extraction =~ m/[^atgcATGC]/){
				print OUT2 "$genome\t$gene_id\t$strand\t$extract_start\t$q_end\t$extract_len\t$extraction\n";
				print SEQOUT ">$genome.$strand$gene_id\n$extraction\n";
			}
			close SEQOUT;
		}
	}
}
close OUT2;
close GENOME;
################### transforming presence/absence matrix to allele matrix
print("Allele calling completed! \nSummarizing results ...\n") if ($opt_v == 0);
open(LIST, "<pa_matrix.transposed.tmp");
open(MATRIX, ">pre_allele_matrix.tmp");
open(GENES, ">gene_list.tmp");
my $gene_sn = 0;
while(<LIST>){
	chomp;
	unless($_ =~ m/^genome|^\s/){
		my $gene = $_;
		$gene =~ s/\t.*//g;
		print GENES "$gene\n";
		my $pa_string = $_;
		$pa_string =~ s/gene[0-9]+\t//;
		my @pa = split("\t", $pa_string);
		$/ = ">";
		open(FFN, "<extracted.$gene.fas");
		$gene_sn++;
		my(%seq_hash, @seqs);
		foreach(<FFN>){
			chomp;
			unless($_ =~ m/^>/){
				my($head, $isolate_name, $seq, $seqlen);
				$head = $_;
				$head =~ s/\n.*//g;
				$isolate_name = $head;
				$isolate_name =~ s/\..*//g;
				$seq = $_;
				$seq =~ s/.*\n//;
				$seq =~ s/>//;
				$seq =~ s/\s//g;
				$seq = lc $seq;
				$seq_hash{$isolate_name} = $seq;
				push @seqs, $seq;
			}
		}
		close FFN;
		system("rm -f extracted.$gene.fas");
		$/ = "\n";
		my %seen = ();
		my @unique_allele = grep {!$seen{$_}++} @seqs;
		my @allele_output = @unique_allele;
		open(NT, ">nt.$gene.fas");
		my(%allele_hash);
		my $allele_number = 0;
		foreach(@allele_output){
			my $allele_length = length($_);
			unless($allele_length eq 0){
				$allele_number++;
				print NT ">$gene", "_$allele_length", "_$allele_number\n$_\n";
				$allele_hash{$_} = $allele_number;
			}
		}
		close NT;
		my $allele_string = $gene;
		my $counter_poz = -1;
		foreach(@pa){
			$counter_poz++;
			my $allele_assigned = ();
			if($_ eq 1){
				my $isolate_name = $isolate_names["$counter_poz"];
				my $seq = $seq_hash{$isolate_name};
				$allele_assigned = $allele_hash{$seq};
			}elsif($_ eq 0){
				my $isolate_name = $isolate_names["$counter_poz"];
				my $seq = $seq_hash{$isolate_name};
				$allele_assigned = $allele_hash{$seq};
				unless($allele_assigned){
					$allele_assigned = "Missing";
				}	
			}elsif($_ eq 3){
				my $isolate_name = $isolate_names["$counter_poz"];
				my $seq = $seq_hash{$isolate_name};
				$allele_assigned = $allele_hash{$seq};
				unless($allele_assigned){
					$allele_assigned = "Truncated";
				}
			}elsif($_ eq 2){
				$allele_assigned = "Paralogue";
			}elsif($_ eq 9){
				$allele_assigned = "Ambiguity_Codes";
			}
			$allele_string = "$allele_string\t$allele_assigned";
		}
		print MATRIX "$allele_string\n";
	}
}
close LIST;
close MATRIX;
close GENES;
open(IN, "<pre_allele_matrix.tmp");
open(OUT, ">>pre2_allele_matrix.tmp");
my($inline, $col_sum, $col);
my(@column);
while (<IN>){
	if($_ =~ m/^\s/){
		next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' pre_allele_matrix.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
system("paste genome_list3.tmp pre2_allele_matrix.tmp >allele_matrix.tmp");
###################
my ($inline, $col_sum, $col);
my (@column);
open(IN, "<allele_matrix.tmp");
open(OUT, ">>allele_matrix.transposed.tmp");
while (<IN>){
	if($_ =~ m/^\s/){
	next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
###################
open(IN, "<allele_matrix.transposed.tmp");
open(OUT, ">allele_matrix.sorted.tmp");
open(OUT2, ">allele_matrix.cleaned.tmp");
my ($inline, $in, $genome, $gene, $gene_previous, $allele, $coordinates, $copy, $counter,
		 $first_col, $second_col, $number_of_column, $allele_sum1, $allele_sum2, $duplicates, $shared);
my (@in, @in2);
my $line_counter = 0;
my $na = 0;
my $dup = 0;
my $same = 0;
my $dif = 0;
while (<IN>){
	$line_counter++;
	chomp;
	if($line_counter == 1){
		next;
	}
	else{
		$inline = $_;
		$inline =~ s/^\s//;
		$genome = $inline;
		$genome =~ s/\t.*//g;
		$inline =~ s/$genome\t//;
		$inline =~ s/\t$//;
		@in = split (/\t/,$inline);
		$number_of_column = @in;
		$allele_sum1 = $in[0] * $number_of_column;
		if($inline =~ /^\s/){
			next;
		}
		else{
			if($inline =~ /Missing/ || $inline =~ /Ambiguity/ || $inline =~ /Truncated/){
				$na++;
				print OUT "$inline\t", "incomplete", "\n";
			}
			elsif($inline =~ /Paralogue/){
				$dup++;
				print OUT "$inline\t", "paralogue", "\n";
			}
			else{
				my($element1, $element2, $judgement);
				$judgement = "identical";
				foreach(@in){
					$element1 = $_;
					foreach(@in){
						$element2 = $_;
						if($element1 != $element2){
							$judgement = "polymorphic";
						}
						else{
							next;
						}
					}
				}
				if($judgement eq "identical"){
					$same++;
					print OUT "$inline\t", "identical", "\n";
					print OUT2 "$inline\n";
				}
				else{
					$dif++;
					print OUT "$inline\t", "polymorphic", "\n";
					print OUT2 "$inline\n";
				}
			}
		}
	}
}
$shared = $same + $dif;
system("paste gene_list.tmp allele_matrix.sorted.tmp >allele_matrix.sorted2.tmp");
my $allele_called = $shared + $dup;
print REPORT "Number of the CDSs in the reference genome: $aa_total<br>",
				"Number of the genomic sequence scanned: $counter_isolates<br>",
				"Number of loci shared by the $counter_isolates genomic sequences: $allele_called<br>",
				"Number of the shared-loci that was found identical: $same<br>",
				"Number of the shared-loci that was excluded because of hypothetical gene duplication: $dup<br>",
				"Number of the shared-loci that was used to construct distance/differnece matrix: $shared ($dif were polymorphic)<br>",
				"Number of loci that was excluded because of incomplete information (missing, truncation or containing nucleotide ambiguity): ", "$na<br>";
close IN;
close OUT;
close OUT2;
###################
open(IN, "<allele_matrix.sorted2.tmp");
open(OUT, ">>allele_matrix.sorted.transposed.tmp");
open(OUT2, ">allele_calling.tmp");
open(LIST, "<genome_list.tmp");
open(LIST2, ">genome_list2.tmp");
print LIST2 "\n";
while(<LIST>){
	chomp;
	print LIST2 "$_\n";
	print OUT2 "\t$_";
}
print OUT2 "\tsummary\n";
system("cat allele_calling.tmp allele_matrix.sorted2.tmp>allele_calling.txt");
my ($inline, $in, $genome, $gene, $gene_previous, $allele, $coordinates, $copy, $counter,
		 $first_col, $second_col, $number_of_column, $allele_sum1, $allele_sum2, $duplicates, $shared);
my (@column, @genome);
while (<IN>){
	if($_ =~ m/^\s/){
		next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.sorted2.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
system("paste genome_list2.tmp allele_matrix.sorted.transposed.tmp >output.allele_matrix.txt");
close IN;
close OUT;
close OUT2;
close LIST;
close LIST2;
###################
open(IN, "<output.allele_matrix.txt");
my ($isolate_key, $profile_value, $gene_list, $type_list, $gene_key, $type_value, $key, $counter,
		$gene, $gene_counter, $type, $allele1, $allele2, $isolate1, $isolate2, $profiles1, $profiles2);
my (@infile, @genes, @isolates, @types, @profile1, @profile2);
my (%gene_hash, %isolate_hash);
while(<IN>){
	chomp;
	if($_ =~ /^\t/){
		next;
	}
	else{
		$_ =~ s/\t$//;
		$isolate_key = $_;
		$isolate_key =~ s/\t.*//g;
		push @isolates, $isolate_key;
		$profile_value =$_;
		$profile_value =~ s/$isolate_key\t//;
				$profile_value =~ s/\n//;
		$isolate_hash{$isolate_key} = $profile_value;
	}
}
open(IN, "<output.allele_matrix.txt");
@infile = <IN>;
$gene_list = $infile[0];
$gene_list =~ s/^\t//;
$gene_list =~ s/\t\n//;
$type_list = $infile[-1];
$type_list =~ s/^\t//;
$type_list =~ s/\t\n//;
@genes = split('\t', $gene_list);
@types = split('\t', $type_list);
$counter = 0;
foreach(@genes){
	$gene_key = $_;
	$type_value = @types["$counter"];
	$gene_hash{$gene_key} = $type_value;
	$counter++;
}
foreach(@isolates){
	$isolate1 = $_;
	$profiles1 = $isolate_hash{"$isolate1"};
	$gene_counter = 0;
	foreach(@genes){
		$gene = $_;
		$type = $types["$gene_counter"];
		$gene_counter++;
		if($type eq "polymorphic"){
			$gene_counter = $gene_counter - 1;
			@profile1 = split ('\t', $profiles1);
			$allele1 = $profile1["$gene_counter"];
			foreach(@isolates){
				$isolate2 = $_;
				$profiles2 = $isolate_hash{"$isolate2"};
				@profile2 = split ('\t', $profiles2);
				$allele2 = $profile2["$gene_counter"];
				if($allele1 == $allele2){
					next;
				}
				else{
					open (PROFILE, ">>pairwise.$isolate1.$isolate2.txt");
					print PROFILE "$gene\t", "$allele1\t", "$allele2\n";
					close PROFILE;
				}
			}
			$gene_counter++;
		}
		else{
			next;
		}
	}
}
###################
my ($inline, $col_sum, $col);
my (@column);
open(IN, "<allele_matrix.cleaned.tmp");
open(OUT, ">>allele_matrix.cleaned.transposed.tmp");
while (<IN>){
	if($_ =~ m/^\s/){
	next;
	}
	else{
		chomp;
		$inline = $_;
		$col_sum = ($inline =~ s/\t/\t/g) + 1;
		$col = 1;
		while($col<=$col_sum){
			@column = readpipe("awk '{print \$$col}' allele_matrix.cleaned.tmp");
			foreach(@column){
				chomp;
				print OUT "$_", "\t";
			}
			$col++;
			print OUT "\n";
		}
		last;
	}
}
close IN;
close OUT;
###################
open(IN_1, "<allele_matrix.cleaned.transposed.tmp");
open(GENOME, "<genome_list.tmp");
open(OUT, ">>dif.tmp");
open(OUT2, ">>dis.tmp");
open(BAPS, ">output.BAPS.txt");
open(STRU, ">Structure.tmp");
my ($isolate_1, $isolate_2, $allele_counter, $dif_counter, $allele1, $allele2, $distance,
		$inline, $outline, $genome_name, $line_counter, $genome_counter, $line_counter, $in);
my (@allele_1, @allele_2, @inline);
$line_counter = 0;
foreach (<IN_1>){
	chomp;
	$line_counter++;
	$in = $_;
	$in =~ s/\t$//;
	print BAPS "$in\t", "$line_counter\n";
	print STRU "$in\n";
	$isolate_1 = $_;
	@allele_1 = split ("\t",$isolate_1);
	open(IN_2, "<allele_matrix.cleaned.transposed.tmp");
	foreach(<IN_2>){
		chomp;
		$allele_counter = 0;
		$dif_counter = 0;
		$isolate_2 = $_;
		@allele_2 = split ("\t",$isolate_2);
		foreach(@allele_2){
			$allele2 = $_;
			$allele1 = @allele_1[$allele_counter];
			$allele_counter++;
			if ($allele1!=$allele2){
			$dif_counter++;
			}
			else{next;
			}
		}
		$distance = $dif_counter/$allele_counter;
		print OUT "$dif_counter\t";
		print OUT2 "$distance\t";
		close IN_2;
	}
		print OUT "\n";
		print OUT2 "\n";
}
close IN_1;
close OUT;
close BAPS;
close STRU;
system("paste genome_list.tmp Structure.tmp >output.Structure.txt");
system("paste genome_list.tmp dis.tmp >dis.txt.tmp");
system("paste genome_list.tmp dif.tmp >dif.txt.tmp");
open DIF, "<dif.txt.tmp";
open DIS, "<dis.txt.tmp";
open DIF2, ">output.difference_matrix.txt";
open DIS2, ">output.Splitstree.nex";
open (GENOME, "<genome_list.tmp");
$genome_counter = 0;
foreach(<GENOME>){
	$genome_counter++;
}
close GENOME;
print DIS2 
"#NEXUS\n\n",
"BEGIN taxa;", "\n",
"   DIMENSIONS ntax=$genome_counter;\n",
"TAXLABELS\n";
print DIF2 "\t";
open GENOME, "<genome_list.tmp";
foreach(<GENOME>){
	if($_ =~ /^\s/){
		next;
	}
	else{
		chomp;
		print DIS2 "   $_\n";
		print DIF2 "$_\t";
	}
}
print DIF2 "\n";
print DIS2
";\n",
"END;\n\n",
"BEGIN distances;\n",
"   DIMENSIONS ntax=", "$genome_counter;\n",
"   FORMAT\n",
"      triangle=LOWER\n",
"      diagonal\n",
"      labels\n",
"      missing=?\n",
"   ;\n",
"MATRIX\n";
$line_counter = 1;
foreach (<DIS>){
	chomp;
	$line_counter++;
	$inline = $_;
	@inline = split ("\t",$inline);
	$genome_name = $inline[0];
	splice @inline, "$line_counter";
	$outline = join("\t", @inline);
	print DIS2 "$outline\n";
}
print DIS2
";\n",
"END;\n";
$line_counter = 1;
foreach (<DIF>){
	chomp;
	$line_counter++;
	$inline = $_;
	@inline = split ("\t",$inline);
	$genome_name = $inline[0];
	splice @inline, "$line_counter";
	$outline = join("\t", @inline);
	print DIF2 "$outline\n";
}
close DIF;
close DIS;
close DIF2;
close DIS2;
###################
open (IN, "<output.difference_matrix.txt");
open (HTML,'>>difference_matrix.tmp');

my ($isolate1, $isolate2, $list_isolate, $line_counter, $inline, $x_counter, $value, $x, $y, $isolates_number, $counter, $start);
my (@infile, @isolates, @values);
@infile = <IN>;
$list_isolate = $infile[0];
$list_isolate =~ s/^\t//;
$list_isolate =~ s/\n//;
@isolates = split ("\t", $list_isolate);
print HTML <<END_OF_HTML;
<HTML><HEAD><TITLE>Difference matrix</TITLE></HEAD><BODY><div><H2><a name="top">Difference matrix</a></H2></div><TABLE BORDER=1><TR>
END_OF_HTML
print HTML "<td><br></td>";
foreach(@isolates){
	print HTML "<td>$_</td>";
}
print HTML "</TR>";
$line_counter = 0;
foreach(@infile){
	chomp;
	if($_ =~ /^\t/){
		next;
	}
	else{
		$inline = $_;
		@values = split ("\t", $inline);
		$y = shift @values;
		$line_counter++;
		print HTML "<TR>", "<td>$y</td>";
		$x_counter = 0;
		foreach(@values){
			$value = $_;
			$x = $isolates["$x_counter"];
			$x_counter++;
			print HTML "<TH><a href=\"../html_$T0/$x.vs.$y.html\">$value</a></TH>"
		}
		print HTML "</TR>";
	}
}
print HTML "</TABLE>";

if($opt_n == 0){
	system("mkdir html_$T0");
	open (GENE, "<output.gene_list.txt");
	my ($alias, $gene, $allele1, $allele2);
	my (@inline);
	my (%gene_hash);
	while(<GENE>){
		chomp;
		$alias = $_;
		$alias =~ s/\t.*//;
		$gene = $_;
		$gene =~ s/.*\t//;
		$gene_hash{$alias} = $gene;
	}
	close GENE;
	$isolates_number = @isolates;
	$counter = 0;
	foreach(@isolates){
		$isolate1 = $_;
		$counter++;
		$start = $counter;
		while($start<$isolates_number){
			$isolate2 = $isolates["$start"];
			$start++;
			open (PAIR, "<pairwise.$isolate1.$isolate2.txt");
			open(HTML2, ">$isolate1.vs.$isolate2.html");
print HTML2 <<END_OF_HTML;
<H3><a name="$isolate1.$isolate2">$isolate1 VS. $isolate2</a></H3><TABLE BORDER=1><TR><TH>alias</TH><TH>$isolate1</TH><TH>$isolate2</TH><TH>product</TH></TR>
END_OF_HTML
			foreach(<PAIR>){
				chomp;
				@inline = split ("\t", $_);
				$alias = $inline[0];
				$gene = $gene_hash{$alias};
				$allele1 = $inline[1];
				$allele2 = $inline[2];
print HTML2 <<END_OF_HTML;
<TR><Td><a href="../aln_$T0/$alias.aln.txt">$alias</a></Td><Th>$allele1</Th><Th>$allele2</Th><Td>$gene</Td></TR>
END_OF_HTML
		}
print HTML2 <<END_OF_HTML;
</TABLE><div><a href="../output_$T0/report.html">index</a></div>
END_OF_HTML

print HTML2 <<END_OF_HTML;
</BODY></HTML>
END_OF_HTML
			system("mv $isolate1.vs.$isolate2.html html_$T0");
			close PAIR;
			close HTML2;
		}
	}
}
print HTML <<END_OF_HTML;
</BODY></HTML>
END_OF_HTML
close HTML;
################### align all the allele sequences
system ("mkdir tmp_$T0");
system ("mkdir output_$T0");
if($opt_l == 1){
	print "producing a 'long' result as requested ..." if ($opt_v == 0);
	open(GENE, "<gene_list.tmp");
	my ($gene, $seqfile, $allele_counter);
	while(<GENE>){
		chomp;
		$gene = $_;
		$seqfile = "nt."."$gene".".fas";
		open (SEQ, "<$seqfile");
		open (OUT, ">$gene.allele.fas");
		$allele_counter = 0;
		while(<SEQ>){
			chomp;
			if($_ =~ />/){
				$_ =~ s/_[0-9]+_/_/;
				$allele_counter ++;
				print OUT "$_\n";
			}
			else{
				print OUT "$_\n";
			}
		}
		close OUT;
		close SEQ;
		if($allele_counter>1){
			system ("mafft --clustalout --quiet $gene.allele.fas >$gene.aln.txt");
			system ("mafft --quiet $gene.allele.fas >$gene.aln.fas");
			system ("rm -f $gene.allele.fas");
		}
		else{
			system ("mv $gene.allele.fas $gene.aln.fas");
		}
	}
	close GENE;
}
################### output core genome files (in FASTA format and eXtended Multi-FASTA format).
if($opt_l == 1){
	my($inline, $line_counter, $gene, $seqname, $allele, $allele_db, $seq, $isolate);
	my(@isolates, @isolates2, @alleles);
	open(IN, "<allele_matrix.transposed.tmp");
	open(OUT, ">allele_profiles.tmp");
	$line_counter = 0;
	while(<IN>){
		chomp;
		$inline = $_;
		$line_counter++;
		if($line_counter == 1){
			print OUT "$inline\n";
		}
		elsif($inline=~ /\t[A-Za-z]/){
			next;
		}
		else{
			print OUT "$inline\n";
		}
	}
	close IN;
	close OUT;
	open(IN, "<allele_profiles.tmp");
	open(OUT, ">clonalframe.dat");
	$line_counter = 0;
	while(<IN>){
		$inline = $_;
		$inline =~ s/\t\n//;
		$line_counter++;
		if($line_counter == 1){
			@isolates = split ("\t", $inline);
			shift @isolates;
			next;
		}
		elsif($inline=~/\t\t/){
			next;
		}
		else{
			@isolates2 = @isolates;
			@alleles = split ("\t", $inline);
			$gene = shift @alleles;
			print OUT "#", "$gene\n";
			foreach(@alleles){
				$isolate = shift @isolates2;
				open(ISOLATE, ">>$isolate.core.fas.tmp");
				$allele = $_;
				open (DB, "<$gene.aln.fas");
				$/=">";
				foreach(<DB>){
					$seqname = $_;
					$seqname =~ s/\n.*//gs;
					$allele_db = $seqname;
					$allele_db =~ s/.*_//g;
					$allele_db =~ s/new//;
					$seq = $_;
					$seq =~ s/$seqname//;
					$seq =~ s/\n//g;
					$seq =~ s/>//;
					if($allele_db == $allele){
						print OUT ">", "$isolate\n", "$seq\n";
						print ISOLATE "$seq\n";
					}
					else{
						next;
					}
				}
				close DB;
				close ISOLATE;
				$/="\n";
			}
		}
		print OUT "=", "\n";
	}
	close OUT;
	open(OUT, ">core_genomes.fas");
	foreach(@isolates){
		$isolate = $_;
		$/="";
		open(IN, "<$isolate.core.fas.tmp");
		foreach(<IN>){
			print OUT ">$isolate.core\n$_";
			system("rm -f $isolate.core.fas.tmp");
		}
		close IN;
	}
	close OUT;
	system ("mv core_genomes.fas ./output_$T0/");
	system ("mv clonalframe.dat ./output_$T0/");
	system ("rm -f gene*.aln.fas");
	system ("mkdir aln_$T0");
	system ("mv *.aln.* ./aln_$T0/");
}
###################
my @files = readpipe("ls");
foreach(<@files>){
	chomp;
	my $filename = $_;
	if($filename =~ m/pairwise\./){
		system "rm -f $filename";
	}
}
system ("rm -f *.allele.fas");
my $T1 = time();
my $time = $T1-$T0;
my $minute = $time/60;
my $minute_rounded = sprintf("%.1f", $minute);
my $hour = $time/3600;
my $hour_rounded = sprintf("%.1f", $hour);
print REPORT "Total running time: $hour_rounded hours (i.e. $minute_rounded minutes, i.e. $time seconds)"; 
close REPORT;
print "All done! Open the folder output_$T0 to view the results.\nThank you for using Fast-GeP! Come again :)\n\n" if ($opt_v == 0);
system ("cat output.report.tmp difference_matrix.tmp >difference_matrix.html");
system ("mv difference_matrix.html ./output_$T0/");
system ("mv output.allele_matrix.txt ./tmp_$T0/allele_matrix.txt");
system ("mv output.difference_matrix.txt ./tmp_$T0/difference_matrix.txt");
system ("mv output.Splitstree.nex ./output_$T0/Splitstree.nex");
system ("mv output.Structure.txt ./output_$T0/allele_profiles.txt");
system ("mv output.gene_list.txt ./output_$T0/gene_list.txt");
system ("mv output.BAPS.txt ./output_$T0/BAPS.txt");
system ("mv allele_calling.txt ./output_$T0/allele_calling.txt");
#system ("mv *.combined.fas ./tmp_$T0/");
system ("rm -f *.combined.fas");
system ("mv *.tmp ./tmp_$T0/");
system ("mv output.* ./output_$T0/");
system ("mv ./output_$T0/output.txt ./tmp_$T0/primary_search.txt");
system ("mv ./output_$T0/output.blastn.txt ./tmp_$T0/complementary_search.txt");
system ("rm -f *.*hr");
system ("rm -f *.*in");
system ("rm -f *.*sq");
system ("rm -f *.dmnd");
system ("rm -f makeblastdb.*");
system ("mkdir scheme_$T0");
system ("mv nt.* ./scheme_$T0");
