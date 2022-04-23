#!/usrbin/env perl
use strict;
use warnings;

my @samples = qw/JGALT1 JGALT2 JGWT1 JGWT2/;
my %info;
open  FIN,"refer/tRNA_mature.len";
while(<FIN>){
	chomp;
	my ($tRNA,$len) = split /\t/,$_;
	$info{$tRNA}{length} = $len;
}
close FIN;
open FIN,"refer/tRNA.info";
while(<FIN>){
	chomp;
	my ($tRNA,$type) = split /\t/,$_;
	my ($x,$y) = split /_/,$type;
	$info{$tRNA}{tRNA_type} = $x;
	$info{$tRNA}{Anticodon} = $y;
}
close FIN;

my %data_a;
my %data_b;
my %data_z;
foreach my $sample(@samples){
	open FIN,"out/tRFs/align/$sample.mature.sam";
	while(<FIN>){
		chomp;
		my ($read_name,$flag,$tRNA,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$other) = split /\t/,$_;
		if($cigar =~ /^(\d+)M$/){
			my $read_len = $1;
			my $start = $pos;
			my $end   = $pos + $read_len - 1;

			#print $read_name,"\n";
			#print $read_len,"\t",$start,"\t",$end,"\n";
			#sleep(1);
			
			$data_z{$tRNA}{$sample}++;
			if($start == 1){
				if(    $end <= 20){
					$data_a{$tRNA}{$sample}{"tRF-5a"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "tRF-5a";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}elsif($end <= 26){
					$data_a{$tRNA}{$sample}{"tRF-5b"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "tRF-5b";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}elsif($end <= 31){
					$data_a{$tRNA}{$sample}{"tRF-5c"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "tRF-5c";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}elsif($end < 38){
					$data_a{$tRNA}{$sample}{"5tiR"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "5tiR";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}else{
					$data_a{$tRNA}{$sample}{"5other"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "5other";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}
			}elsif($end == $info{$tRNA}{length}){
				if(    $start >= $end - 19){
					$data_a{$tRNA}{$sample}{"tRF-3a"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "tRF-3a";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}elsif($start >= $end - 26){
					$data_a{$tRNA}{$sample}{"tRF-3b"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "tRF-3b";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}elsif($start > 31 && $start < 38){
					$data_a{$tRNA}{$sample}{"3tiR"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "3tiR";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}else{
					$data_a{$tRNA}{$sample}{"3other"}++;
					$data_b{$seq}{$tRNA}{expr}{$sample}++;
					$data_b{$seq}{$tRNA}{type} = "3other";
					$data_b{$seq}{$tRNA}{site} = $start."-".$end;
				}
			}elsif(( $start > 26 && $start <= 31 ) && ($end > 38 &&  $end <= 43)){
				$data_a{$tRNA}{$sample}{"tRF-2"}++;
				$data_b{$seq}{$tRNA}{expr}{$sample}++;
				$data_b{$seq}{$tRNA}{type} = "tRF-2";
				$data_b{$seq}{$tRNA}{site} = $start."-".$end;
			}
		}
	}
	close FIN;
	open FIN,"out/tRFs/align/$sample.down100.sam";
	while(<FIN>){
		chomp;
		my ($read_name,$flag,$tRNA,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$seq,$qual,$other) = split /\t/,$_;
		if($cigar =~ /^(\d+)M$/){
			my $read_len = $1;
			my $start = $pos;
			my $end   = $pos + $read_len - 1;

			if($start == 1){
				$data_a{$tRNA}{$sample}{"tRF-1"}++;
				$data_b{$seq}{$tRNA}{expr}{$sample}++;
				$data_b{$seq}{$tRNA}{type} = "tRF-1";
				$data_b{$seq}{$tRNA}{site} = "downstrem:".$start."-".$end;
			}
		}
	}
	close FIN;
}

open  FOUT,">out/tRFs/tRNA.expression.xls";
print FOUT "tRNA\t".join("\t",@samples)."\n";
foreach my $tRNA (%data_z){
	if($tRNA =~ /^HASH/){
		next;
	}
	print FOUT $tRNA;
	foreach my $sample(@samples){
		if(exists $data_z{$tRNA}{$sample}){
			print FOUT "\t",$data_z{$tRNA}{$sample};
		}else{
			print FOUT "\t",0;
		}
	}
	print FOUT "\n";

}
close FOUT;

open  FOUT,">out/tRFs/tRFs.expression.detail.xls";
print FOUT "tRF/tiRNA_Sequence\ttRF/tiRNA_Type\ttRF/tiRNA_Region\tSource_tRNA\ttRNA_Length\ttRNA_Type\tAnticodon\t".join("\t",@samples)."\n";
foreach my $seq(sort keys %data_b){
	foreach my $tRNA(sort keys %{$data_b{$seq}}){
		print FOUT $seq,"\t",$data_b{$seq}{$tRNA}{type},"\t",$data_b{$seq}{$tRNA}{site},"\t",$tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
		foreach my $sample(@samples){
			if(exists $data_b{$seq}{$tRNA}{expr}{$sample}){
				print FOUT "\t",$data_b{$seq}{$tRNA}{expr}{$sample};
			}else{
				print FOUT "\t",0;
			}
		}
		print FOUT "\n";
	}
}
close FOUT;

__END__
open  F5, ">out/tRFs/5-tRFs.expression.detail.xls";
open  F3, ">out/tRFs/3-tRFs.expression.detail.xls";
open  FM, ">out/tRFs/mid-tRFs.expression.detail.xls";
open  F55,">out/tRFs/5-tiRNA.expression.detail.xls";
open  F33,">out/tRFs/3-tiRNA.expression.detail.xls";
print F5 "tRNA_ID\ttRNA_length\ttRNA_type\tAnticodon\t".join("\t",qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/)."\n";
print F3 "tRNA_ID\ttRNA_length\ttRNA_type\tAnticodon\t".join("\t",qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/)."\n";
print FM "tRNA_ID\ttRNA_length\ttRNA_type\tAnticodon\t".join("\t",qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/)."\n";
print F55 "tRNA_ID\ttRNA_length\ttRNA_type\tAnticodon\t".join("\t",qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/)."\n";
print F33 "tRNA_ID\ttRNA_length\ttRNA_type\tAnticodon\t".join("\t",qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/)."\n";
foreach my $tRNA(sort keys %info){
	print F5 $tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
	print F3 $tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
	print FM $tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
	print F55 $tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
	print F33 $tRNA,"\t",$info{$tRNA}{length},"\t",$info{$tRNA}{tRNA_type},"\t",$info{$tRNA}{Anticodon};
	foreach my $sample(qw/WT_CK_1_b WT_CK_2_b WT_CK_3_b WT_HS_1_b WT_HS_2_b WT_HS_3_b mop1_CK_1_b mop1_CK_2_b mop1_CK_3_b mop1_HS_1_b mop1_HS_2_b mop1_HS_3_b/){
		my $num_5a = ( defined($data{$tRNA}{$sample}{"5a"}) ? $data{$tRNA}{$sample}{"5a"} : 0 );
		my $num_5b = ( defined($data{$tRNA}{$sample}{"5b"}) ? $data{$tRNA}{$sample}{"5b"} : 0 );
		my $num_5c = ( defined($data{$tRNA}{$sample}{"5c"}) ? $data{$tRNA}{$sample}{"5c"} : 0 );
		my $num_5d = ( defined($data{$tRNA}{$sample}{"5d"}) ? $data{$tRNA}{$sample}{"5d"} : 0 );

		my $num_3a = ( defined($data{$tRNA}{$sample}{"3a"}) ? $data{$tRNA}{$sample}{"3a"} : 0 );
		my $num_3b = ( defined($data{$tRNA}{$sample}{"3b"}) ? $data{$tRNA}{$sample}{"3b"} : 0 );
		my $num_3d = ( defined($data{$tRNA}{$sample}{"3d"}) ? $data{$tRNA}{$sample}{"3d"} : 0 );

		my $num_md = ( defined($data{$tRNA}{$sample}{"md"}) ? $data{$tRNA}{$sample}{"md"} : 0 );

		print F5 "\t",$num_5a+$num_5b+$num_5c;
		print F3 "\t",$num_3a+$num_3b;
		print FM "\t",$num_md;
		print F55 "\t",$num_5d;
		print F33 "\t",$num_3d;
	}
	print F5 "\n";
	print F3 "\n";
	print FM "\n";
	print F55 "\n";
	print F33 "\n";
}
close F5;
close F3;
close FM;
close F55;
close F33;
