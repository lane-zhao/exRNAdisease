#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use PBS::Queue;

my $version  = "v2.0";
my $script   = "$Bin/script";
my $software = "$Bin/software";
my $database = "$Bin/database"; 

my $queue;
my $workdir  = getcwd;

my ($help, $file_input, $file_meta, $file_comp, $outdir, $strand);
my ($refer, $gtf, $go_asso, $kegg_asso, $novel, $de_pval, $de_fdr, $de_fc, $de_log2fc, $gene_type);
my ($species);
my ($noqc, $circrna);
my ($max_cpus, $max_jobs, $queue_name, $job_prefix, $job_local, $queue_rerun);
my %bkg_info;
my (%meta_grp2sn,%meta_sn2grp);
my (%comp_grp2sn,%comp_sn2grp, %comp_param);
my (%samples,%samples_type);
my %all_dependencies;
my %results;
my %meta;
my %group_info;
my @sample_names;

&GetOptions(
	"help!"          => \$help,
	"input:s"        => \$file_input,
	"meta:s"         => \$file_meta,
	"comp:s"         => \$file_comp,
	"outdir:s"       => \$outdir,
	
	"refer:s"        => \$refer,
	"gtf:s"          => \$gtf,
	"go_asso:s"      => \$go_asso,
	"kegg_asso:s"    => \$kegg_asso,
	"species:s"      => \$species,
	
	"strand:s"       => \$strand,
	
	"novel!"         => \$novel,
	"noqc!"          => \$noqc,
	"circrna!"       => \$circrna,
	
	"gene_type:s"    => \$gene_type,
	"de_pval:f"      => \$de_pval,
	"de_fdr:f"       => \$de_fdr,
	"de_fc:f"        => \$de_fc,
	
	"max_cpus:i"     => \$max_cpus,
	"max_jobs:i"     => \$max_jobs,
	"queue_name:s"   => \$queue_name,
	"job_prefix:s"   => \$job_prefix,
	"local!"         => \$job_local,
	"rerun!"         => \$queue_rerun,
);

sub usage{
	my $program = basename($0);
	die "
Version: $version
Aurthor: Guantao Zheng
Contact: guantao.zheng\@origin-gene.com
Create : 20190418 v1.0
Update : 20190516 v1.1
         20200702 v2.0
Note   : Before this workflow, you must do:
         1. check the consistency of contig or chr id between genome and gtf file!!!
         2. sort the gtf file!!!
Usage  : $program [options]
Options:
  Input Options:
    --input          <string>    fastq file list.[require]
    --meta           <string>    meta information for grouping comparison. [ optional ]
    --comp           <string>    comparison configure file. [ optional ]
    --refer          <string>    reference genome sequence. [ required ]
    --gtf            <string>    gene annotation file in gtf format. [ required ]
    --strand         <string>    non-strand-specific(no,default) or strand-specific(RF)

  Output Options:
    --outdir         <string>    output dir, default is 'out'.

  Analysis Options:
    --go_asso        <string>    association file between gene and go term. [optional]
    --kegg_asso      <string>    association file between gene and go term. [optional]
    --species        <string>    species name in kegg format, default is ko
    --gene_type      <string>    gene to type conversion file.
    --de_pval        <float>     the threshold of de pval, default is 0
                                 if it set 0, this option will be missed
    --de_fdr         <string>    the threshold of de pdr,  default is 0.05
                                 if it set 0, this option will be missed
    --de_fc          <number>    the threshold of fold change, default is 2
                                 if it set 0, this option will be missed

  Step Control Options:
    -noqc           <NA>        skip fastq quality control
    -circrna        <NA>        run circRNA analysis based on lncRNA data
    
  Task Options:
    -max_cpus       <number>    max cpu number limitation, default is 200 [ optional ]
    -max_jobs       <number>    max job number limitation, default is 20  [ optional ]
    -queue_name     <string>    queue name, default is yxsw
    -local          <NA>        local model
    -rerun          <NA>        rerun model
###################
#P#PS1:<input> contain 3 colums, sample_id, read1.fq(gz), read2.fq(gz)
#		sample_1	/path/sample_1_R1.fq	/path/sample_1_R2.fq
#		sample_2	/path/sample_2_R1.fq	/path/sample_2_R2.fq
#		sample_3	/path/sample_3_R1.fq	/path/sample_3_R2.fq
#		......
#		For different lines, samples name must be unique !!!
#PS2:[meta] contain at leat 2 colums, 1st is [Sample], 2nd, 3rd, 4th ... are [Group Standard]
#		SampleID	Group_STD_1	Group_STD_2	Group_STD_3
#		sample_1	g1	w	A
#		sample_2	g1	d	A
#		sample_3	g2	d	
#		sample_4	g2	w	B
#		sample_5	g3	d	B
#		sample_6	g3	w	
#		For different lines, samples name must be unique !!!
#PS3:[comp] 1st is [Sample] or [Group Standard], 2nd, 3rd, 4th ... are group name
#		Sample	sample_1,sample_3
#		Group_STD_1	g1,g3
#		Group_STD_1	g1,g2,g3
#		Group_STD_2	w,d
#		For different lines, samples name must be unique !!!
#PS4:[gene_type] gene type file
#		gene_1	mRNA
#		gene_2	mRNA
#		gene_3	sRNA
###################		
\n";
}

&parsing_parameters();
&parsing_sample_info;
&parsing_meta_and_comp_info;
&fastq_qc();
&run_rockhopper();
&gene_annotation();
&sRNA_analysis();
&refer_preparation();
&alignment();
&explevel();
&sample_cluster();
&expdiff_analysis();
&promoter();
&snp_indel();
&transcriptome_assesement();
&result_collection();
$queue->jointhreads();

print "all jobs done!\n";

###########################################################################

sub parsing_parameters{

	&usage if ($help);
	
	if(!defined( $file_input )){
		print STDERR "Error: -input must be specified!\n";
		&usage;
	}elsif(!(-e $file_input)){
		print STDERR "Error: input [$file_input] must be exist!\n";
		&usage;
	}

	if(!defined( $refer )){
		print STDERR "Error: -refer must be specified!\n";
		&usage;
	}elsif(!(-e $refer)){
		print STDERR "Error: refer [$refer] must be exist!\n";
		&usage;
	}else{
		$refer = abs_path($refer);
	}
	
	if( !defined( $gtf ) ){
		print STDERR "Error:-gtf must be specified!\n";
		&usage;
	}elsif(!(-e $gtf)){
		print STDERR "Error: gtf [$gtf] must be exist!\n";
		&usage;
	}else{
		$gtf = abs_path($gtf);
	}
	
	# if( !defined( $go_asso ) ){
		# print STDERR "Error:-go must be specified!\n";
		# &usage;
	# }elsif(!(-e $go_asso)){
		# print STDERR "Error: go association file  [$go_asso] must be exist!\n";
		# &usage;
	# }else{
		# $go_asso = abs_path($go_asso);
	# }
	
	# if( !defined( $kegg_asso ) ){
		# print STDERR "Error:-kegg must be specified!\n";
		# &usage;
	# }elsif(!(-e $kegg_asso)){
		# print STDERR "Error: kegg association file  [$kegg_asso] must be exist!\n";
		# &usage;
	# }else{
		# $kegg_asso = abs_path($kegg_asso);
	# }
	
	$species ||= "ko";
	$strand  ||= "no";
	$de_pval   = 0     unless(defined $de_pval);
	$de_fc     = 2     unless(defined $de_fc);
	$de_fdr    = 0.05  unless(defined $de_fdr);
	if($de_fc > 0){
		$de_log2fc = log($de_fc)/log(2);
	}else{
		$de_log2fc = 0;
	}
	$queue_name  ||= "yxsw";
	$job_prefix  ||= "prna";
	$queue_rerun ||= 0;
	$job_local   ||= 1;
	$max_cpus    ||= 200;
	$max_jobs    ||= 40;
	$outdir      ||= "out";
	$outdir = abs_path( $outdir );
	
	if($strand ne "no" && $strand ne "yes"){
		print STDERR "Error:-strand must be yes or no!\n";
		&usage;
	}
	
	if($de_pval == 0 && $de_fdr == 0 && $de_fc == 0){
		print STDERR "Error: at least one of -de_pval, -de_fdr, and -de_fc must be none-zero !\n";
		&usage;
	}

	system("mkdir -p $outdir");
	
	$queue = PBS::Queue->new(
		{
			'queue_name'  => $queue_name, 
			'job_prefix'  => $job_prefix,
			'max_cpus'    => $max_cpus,
			'max_jobs'    => $max_jobs,
			'job_local'   => $job_local,
			'queue_rerun' => $queue_rerun,
		}
	);
}

sub parsing_sample_info{
	open  FIN,$file_input or die "can not open raw fq list file!";
	while(<FIN>){
		chomp;
		next if /^#/;
		my ($sample,$fq1,$fq2) = split /\t/,$_;
		if(exists $samples{$sample}){
			die "sample_name [$sample] is duplicate! please check input file[$file_input]!";
		}
		$samples{$sample}{R1} = $fq1;
		$samples{$sample}{R2} = $fq2;
		push @sample_names,$sample;
	}
	close FIN;
}

sub parsing_meta_and_comp_info{
	system("mkdir -p $outdir/group_info");
	my %mat_info;
	my %all_std;
	if(defined $file_meta && defined $file_comp){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		my $head = <FIN>; chomp($head);
		unless($head =~ "^SampleID\t"){
			die "Meta file format error[header], please check!";
		}
		my ($a,@gs) = split /\t/,$head;
		while(<FIN>){
			chomp;
			my ($sample,@G) = split /\t/,$_;
			unless(exists $samples{$sample}){
				die "sample_name[$sample] in meta file does not existe in fasta file!";
			}
			if(scalar @G != scalar @gs){
				die "meta file format error[content], please check!";
			}
			for(my $i = 0; $i < scalar(@gs); $i++){
				if($G[$i] ne ""){
					$meta{$gs[$i]}{$G[$i]}{$sample}++;
				}
			}
			$meta{"Sample"}{$sample}{$sample}++;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			next if /^#/;
			my ($name,$class,$tmp) = split /\t/,$_;
			unless(exists $meta{$class}){
				die "group std name [$class] does not exists in meta info file [$file_meta]!";
			}
			my @G =split /,/,$tmp;
			if( scalar(@G) == 1 ){
				die "the number of elements in comparison should not be 1!";
			}elsif( scalar(@G) == 0 ){
				@G = sort keys %{$meta{$class}};
				$name = $class.".all" unless $name;
			}else{
				$name = $class.".".join("_vs_",@G) unless $name;
			}
			if(exists $all_std{$name}){
				die "std[$name] was exists!!!";
			}
			$all_std{$name}++;
			foreach my $ele(@G){
				die "group name [$ele] does not exists in group std [$class] !!" unless exists $meta{$class}{$ele};
				foreach my $sam(sort keys %{$meta{$class}{$ele}}){
					$group_info{$name}{$ele}{$sam}++;
					$mat_info{$sam}{$name} = $ele;
				}
			}
		}
		close FIN;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "sample\tgroup\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
	}elsif(defined $file_meta && !(defined $file_comp)){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		# SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription
		my $head = <FIN>; chomp($head);
		unless($head =~ "^SampleID\t"){
			die "Meta file format error, please check!";
		}
		my ($a,@gs) = split /\t/,$head;
		while(<FIN>){
			chomp;
			my ($sample,@G) = split /\t/,$_;
			unless(exists $samples{$sample}){
				die "sample_name[$sample] in meta file does not existe in fasta file!";
			}
			if(scalar @G != scalar @gs){
				die "meta file format error, please check!";
			}
			for(my $i = 0; $i < scalar(@gs); $i++){
				if($G[$i] ne ""){
					$group_info{$gs[$i]}{$G[$i]}{$sample}++;
					$mat_info{$sample}{$gs[$i]} = $G[$i];
					$all_std{$gs[$i]}++;
				}
			}
			$group_info{"Sample"}{$sample}{$sample}++;
			$mat_info{$sample}{"Sample"} = $sample;
			$all_std{"Sample"}++;
		}
		close FIN;
		close FOUT;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "#id\t$class\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
	}elsif(defined $file_comp && !(defined $file_meta)){
		open FIN,"$file_input";
		while(<FIN>){
			chomp;
			my ($sample) = split /\t/,$_;
			$meta{"Sample"}{$sample}{$sample}++;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			my ($name,$class,$tmp) = split /\t/,$_;
			if($class ne "Sample"){
				die "class must be \"Sample\"!";
			}
			my @G =split /,/,$tmp;
			if( scalar(@G) == 1 ){
				die "the number of elements in comparison should not be 1!";
			}elsif( scalar(@G) == 0 ){
				@G = sort keys %{$meta{$class}};
				$name = $class.".all" unless $name;
			}else{
				$name = $class.".".join("_vs_",@G) unless $name;
			}
			$all_std{$name}++;
			foreach my $ele(@G){
				die "group name [$ele] does not exists in group std [$class] !!" unless exists $meta{$class}{$ele};
				foreach my $sam(sort keys %{$meta{$class}{$ele}}){
					$group_info{$name}{$ele}{$sam}++;
					$mat_info{$sam}{$name} = $ele;
				}
			}
		}
		close FIN;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "#id\t$class\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
		
	}else{
		open FIN,"$file_input";
		open F1,">$outdir/group_info/Sample_all.info";
		open F2,">$outdir/group_info/Sample_all.txt";
		print F1 "#id\tSample\n";
		while(<FIN>){
			chomp;
			my ($sample) = split /\t/,$_;
			$group_info{"Sample_all"}{$sample}{$sample}++;
			$mat_info{$sample}{"Sample_all"} = $sample;
			print F1 $sample,"\t",$sample,"\n";
			print F2 $sample,"\t",$sample,"\n";
		}
		close F1;
		close F2;
		close FIN;
		
		open FOUT,">$outdir/group_info/group_info.tsv";
		print FOUT "SampleID\tSample_all\n";
		foreach my $sample(sort keys %{$group_info{"Sample_all"}}){
			print FOUT "$sample\t$sample\n";
		}
		close FOUT;
	}
	$file_meta = "$outdir/group_info/group_info.tsv";
}

sub rawdata_quality_control{
	foreach my $sn(sort keys %samples){
		$all_dependencies{"rawdata_qc.sample_$sn"}++;
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("rawdata_qc.sample_$sn");
		$queue->set_job_disc("样本[$sn]原始数据质量控制分析");
		$queue->set_job_depend("");
		$queue->set_file_check("$outdir/rawdata_qc/$sn/$sn.fastp.html,$outdir/rawdata_qc/$sn/$sn.fastp.json");
		$queue->set_work_dir($workdir);
		
		$queue->addcommonds("module load fastp/0.21.0");
		$queue->addcommonds("module load seqtk/1.3");
		$queue->addcommonds("module load ucsc/1.0");
		$queue->addcommonds("module load fastqc/0.11.8");
		$queue->addcommonds("module load blast/2.2.26");
		
		$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/");
		
		if($samples_type{$sn} eq "SE"){
			my @fqs = split /,/,$samples{$sn}{S};
			if(scalar(@fqs) == 1){
				$queue->addcommonds("fastp --in1 $fqs[0] --out1 $outdir/rawdata_qc/$sn/$sn.cleandata.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
				$queue->addcommonds("$script/parse_fastp.pl --json $outdir/rawdata_qc/$sn/$sn.fastp.json --outdir $outdir/rawdata_qc/$sn --prefix $sn");
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.cleandata.fq 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
				foreach my $x(qw/rawdata cleandata/){
					$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/$sn.$x.fastqc");
					$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.$x.fq --outdir $outdir/rawdata_qc/$sn/$sn.$x.fastqc &> $outdir/rawdata_qc/$sn/$sn.$x.fastqc.log");
				}
			}else{
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.rawdata.fq");
				foreach my $fq(@fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.fq");
					}
				}
				$queue->addcommonds("fastp --in1 $outdir/rawdata_qc/$sn/$sn.rawdata.fq --out1 $outdir/rawdata_qc/$sn/$sn.cleandata.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
				$queue->addcommonds("$script/parse_fastp.pl --json $outdir/rawdata_qc/$sn/$sn.fastp.json --outdir $outdir/rawdata_qc/$sn --prefix $sn");
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.cleandata.fq 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
				foreach my $x(qw/rawdata cleandata/){
					$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/$sn.$x.fastqc");
					$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.$x.fq --outdir $outdir/rawdata_qc/$sn/$sn.$x.fastqc &> $outdir/rawdata_qc/$sn/$sn.$x.fastqc.log");
				}
			}
		}else{
			my @R1_fqs = split /,/,$samples{$sn}{R1};
			my @R2_fqs = split /,/,$samples{$sn}{R2};
			
			if(scalar(@R1_fqs) != scalar(@R2_fqs)){
				die "[Error]: sample[$sn] has diff number of files in R1 and R2!\n";
			}elsif(scalar(@R1_fqs) == 1){
				$queue->addcommonds("fastp --in1 $R1_fqs[0] --in2 $R2_fqs[0] --out1 $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz --out2 $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.R1.fq.gz --unpaired2 $outdir/rawdata_qc/$sn/$sn.unpair.R2.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --correction --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report' > $outdir/rawdata_qc/$sn/$sn.fastp.log");
				$queue->addcommonds("$script/parse_fastp.pl --json $outdir/rawdata_qc/$sn/$sn.fastp.json --outdir $outdir/rawdata_qc/$sn --prefix $sn");
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
				
				foreach my $x(qw/rawdata cleandata/){
					foreach my $y(qw/R1 R2/){
						$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/$sn.$x.$y.fastqc");
					}
				}
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $R1_fqs[0] --outdir $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fastqc > $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $R2_fqs[0] --outdir $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fastqc > $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz --outdir $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fastqc > $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fq.gz --outdir $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fastqc > $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fastqc.log");
				
			}else{
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fq");
				foreach my $fq(@R1_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fq");
					}
				}
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fq");
				foreach my $fq(@R2_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fq");
					}
				}
				$queue->addcommonds("fastp --in1 $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fq --in2 $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fq --out1 $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz --out2 $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.R1.fq.gz --unpaired2 $outdir/rawdata_qc/$sn/$sn.unpair.R2.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --correction --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report' > $outdir/rawdata_qc/$sn/$sn.fastp.log");
				$queue->addcommonds("$script/parse_fastp.pl --json $outdir/rawdata_qc/$sn/$sn.fastp.json --outdir $outdir/rawdata_qc/$sn --prefix $sn");
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
				foreach my $x(qw/rawdata cleandata/){
					foreach my $y(qw/R1 R2/){
						$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/$sn.$x.$y.fastqc");
					}
				}
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fq --outdir $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fastqc &> $outdir/rawdata_qc/$sn/$sn.rawdata.R1.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fq --outdir $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fastqc &> $outdir/rawdata_qc/$sn/$sn.rawdata.R2.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fq.gz --outdir $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fastqc &> $outdir/rawdata_qc/$sn/$sn.cleandata.R1.fastqc.log");
				$queue->addcommonds("fastqc -t 10 --nogroup --extract $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fq.gz --outdir $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fastqc &> $outdir/rawdata_qc/$sn/$sn.cleandata.R2.fastqc.log");
			}
		}
		$queue->run();
		$results{rawdata_quality_control}++;
	}
	
	$all_dependencies{"rawdata_qc.merge"}++;
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_mode($job_mode);
	$queue->set_job_name("rawdata_qc.merge");
	$queue->set_job_disc("所有样本原始数据质量控制分析整合");
	my @p = ();
	foreach my $sn(sort keys %samples){
		push @p, "rawdata_qc.sample_$sn";
	}
	$queue->set_job_depend(join(",",@p));
	$queue->set_file_check("$outdir/rawdata_qc/cleandata.stat.xls,$outdir/rawdata_qc/rawdata.stat.xls,$outdir/rawdata_qc/rRNA_percentage.xls,$outdir/rawdata_qc/summary.filtering.xls,$outdir/rawdata_qc/summary.qc.xls");
	$queue->set_work_dir($workdir);

	my $cmd;
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/rawdata_qc/rawdata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.rawdata.stat.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/rawdata_qc/cleandata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.cleandata.stat.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/rawdata_qc/summary.filtering.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.filter_summary.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/rawdata_qc/summary.qc.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.qc_summary.xls";
	}
	$queue->addcommonds("$cmd");
	
	$queue->addcommonds("grep '^rRNA' $outdir/rawdata_qc/*/*.rfam_summary.xls | sed 's#$outdir/rawdata_qc/##g' | cut -f 2 -d '/' | sed 's/.rfam_summary.xls:rRNA//g' | cut -f 1,3 | sed '1 isample\trRNA_Rate(%)' > $outdir/rawdata_qc/rRNA_percentage.xls");
	
	#$queue->addcommonds("$script/qc_summary.pl --rawdata $outdir/rawdata_qc/rawdata.stat.xls -cleandata $outdir/rawdata_qc/cleandata.stat.xls -rRNA $outdir/rawdata_qc/rRNA_percentage.xls -o $outdir/rawdata_qc/QC_summary.xls");
	$queue->run();
	
	$results{rawdata_quality_control}++;
}


sub rawdata_qc{
	foreach my $sn(sort keys %samples){
		$all_dependencies{"rawdata_qc.sample_$sn"}++;
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("rawdata_qc.sample_$sn");
		$queue->set_job_disc("样本[$sn]原始数据质量控制分析");
		$queue->set_job_depend("");
		$queue->set_file_check("$outdir/rawdata_qc/$sn/$sn.fastp.html,$outdir/rawdata_qc/$sn/$sn.fastp.json");
		$queue->set_work_dir($workdir);
		$queue->addcommonds("module load fastp/0.21.0");
		$queue->addcommonds("module load seqtk/1.3");
		$queue->addcommonds("module load ucsc/1.0");
		$queue->addcommonds("module load blast/2.2.26");
		
		$queue->addcommonds("mkdir -p $outdir/rawdata_qc/$sn/");
		
		if($samples_type{$sn} eq "SE"){
			my @fqs = split /,/,$samples{$sn}{S};
			if(scalar(@fqs) == 1){
				$queue->addcommonds("fastp --in1 $fqs[0] --out1 $outdir/rawdata_qc/$sn/$sn.clean.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
			}else{
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.raw.fq");
				foreach my $fq(@fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.raw.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.raw.fq");
					}
				}
				$queue->addcommonds("fastp --in1 $outdir/rawdata_qc/$sn/$sn.raw.fq --out1 $outdir/rawdata_qc/$sn/$sn.clean.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
			}
		}else{
			my @R1_fqs = split /,/,$samples{$sn}{R1};
			my @R2_fqs = split /,/,$samples{$sn}{R2};
			
			if(scalar(@R1_fqs) != scalar(@R2_fqs)){
				die "[Error]: sample[$sn] has diff number of files in R1 and R2!\n";
			}elsif(scalar(@R1_fqs) == 1){
				$queue->addcommonds("fastp --in1 $R1_fqs[0] --in2 $R2_fqs[0] --out1 $outdir/rawdata_qc/$sn/$sn.clean.R1.fq.gz --out2 $outdir/rawdata_qc/$sn/$sn.clean.R2.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.R1.fq.gz --unpaired2 $outdir/rawdata_qc/$sn/$sn.unpair.R2.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --correction --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
			}else{
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.raw.R1.fq");
				foreach my $fq(@R1_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.raw.R1.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.raw.R1.fq");
					}
				}
				$queue->addcommonds("rm -f $outdir/rawdata_qc/$sn/$sn.raw.R2.fq");
				foreach my $fq(@R2_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/rawdata_qc/$sn/$sn.raw.R2.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/rawdata_qc/$sn/$sn.raw.R2.fq");
					}
				}
				$queue->addcommonds("fastp --in1 $outdir/rawdata_qc/$sn/$sn.raw.R1.fq --in2 $outdir/rawdata_qc/$sn/$sn.raw.R2.fq --out1 $outdir/rawdata_qc/$sn/$sn.clean.R1.fq.gz --out2 $outdir/rawdata_qc/$sn/$sn.clean.R2.fq.gz --unpaired1 $outdir/rawdata_qc/$sn/$sn.unpair.R1.fq.gz --unpaired2 $outdir/rawdata_qc/$sn/$sn.unpair.R2.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --correction --length_required 20 --thread 10 --json $outdir/rawdata_qc/$sn/$sn.fastp.json --html $outdir/rawdata_qc/$sn/$sn.fastp.html --report_title '$sn fastp report'");
				
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.clean.R1.fq.gz 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
				
				$queue->addcommonds("seqtk sample $outdir/rawdata_qc/$sn/$sn.clean.R1.fq 10000 > $outdir/rawdata_qc/$sn/$sn.sample.fq");
				$queue->addcommonds("fastqToFa $outdir/rawdata_qc/$sn/$sn.sample.fq $outdir/rawdata_qc/$sn/$sn.sample.fa");
				$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/rawdata_qc/$sn/$sn.sample.fa -o $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
				$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt > $outdir/rawdata_qc/$sn/$sn.rfam_blast.xls");
				$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/rawdata_qc/$sn/$sn.rfam_blast.txt -fasta $outdir/rawdata_qc/$sn/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/rawdata_qc/$sn/ --prefix $sn -maxEvalue 0.01");
			}
		}
		$queue->run();
		
		$results{rawdata_qc}++;
	}
	
	
	foreach my $sn(sort keys %samples){
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_name("样本[$sn]数据质控与统计分析");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("module load seqtk/1.2");
		$queue->addcommond("module load ucsc/1.0");
		$queue->addcommond("module load fastqc/0.11.4");
		$queue->addcommond("module load blast/2.2.26");
		
		$queue->addcommond("mkdir -p $outdir/fastq_qc/$sn");
		foreach my $x(qw/R1 R2/){
			my @fqs = split /,/,$samples{$sn}{$x};
			$queue->addcommond("rm -f $outdir/fastq_qc/$sn/$sn.$x.fq");
			foreach my $fq(@fqs){
				if($fq =~ /.gz$/){
					$queue->addcommond("zcat $fq >> $outdir/fastq_qc/$sn/$sn.$x.fq");
				}else{
					$queue->addcommond("cat  $fq >> $outdir/fastq_qc/$sn/$sn.$x.fq");
				}
			}
		}
		
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_qc/$sn/$sn.R1.fq\\t$outdir/fastq_qc/$sn/$sn.R2.fq\\n' > $outdir/fastq_qc/$sn/$sn.rawfq.list");
		$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_qc/$sn/$sn.rawfq.list > $outdir/fastq_qc/$sn/$sn.rawdata.stat");
		$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_qc/$sn/$sn.rawdata.stat > $outdir/fastq_qc/$sn/$sn.rawdata.stat.xls");
		$queue->addcommond("fastqc -t 10 --nogroup --extract $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq &> $outdir/fastq_qc/$sn/$sn.fastqc.log");
		
		$queue->addcommond("cutadapt -q 20 -e 0.1 -n 1 -m 20 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o $outdir/fastq_qc/$sn/$sn.trim.R1.fq -p $outdir/fastq_qc/$sn/$sn.trim.R2.fq $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq &> $outdir/fastq_qc/$sn/$sn.cutadapt.log");
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_qc/$sn/$sn.trim.R1.fq\\t$outdir/fastq_qc/$sn/$sn.trim.R2.fq\\n' > $outdir/fastq_qc/$sn/$sn.trimfq.list");
		$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_qc/$sn/$sn.trimfq.list > $outdir/fastq_qc/$sn/$sn.trimdata.stat");
		$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_qc/$sn/$sn.trimdata.stat > $outdir/fastq_qc/$sn/$sn.trimdata.stat.xls");
		$queue->addcommond("fastqc -t 10 --nogroup --extract $outdir/fastq_qc/$sn/$sn.trim.R1.fq $outdir/fastq_qc/$sn/$sn.trim.R2.fq &> $outdir/fastq_qc/$sn/$sn.trim.fastqc.log");
		
		my $rfam_seed = "/home/pub/database/blastdb/blast/Rfam.seed";
		$queue->addcommond("seqtk sample $outdir/fastq_qc/$sn/$sn.trim.R1.fq 10000 > $outdir/fastq_qc/$sn/$sn.sample.fq");
		$queue->addcommond("fastqToFa $outdir/fastq_qc/$sn/$sn.sample.fq $outdir/fastq_qc/$sn/$sn.sample.fa");
		$queue->addcommond("blastall -p blastn -d Rfam -i $outdir/fastq_qc/$sn/$sn.sample.fa -o $outdir/fastq_qc/$sn/$sn.sample_vs_Rfam.txt -v 5 -b 5 -F F -e 0.01 -a 10");
		$queue->addcommond("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/fastq_qc/$sn/$sn.sample_vs_Rfam.txt > $outdir/fastq_qc/$sn/$sn.sample_vs_Rfam.xls");
		$queue->addcommond("$script/rfam_blast_parse.pl -i $outdir/fastq_qc/$sn/$sn.sample_vs_Rfam.txt -fa $outdir/fastq_qc/$sn/$sn.sample.fa -s $rfam_seed -o $outdir/fastq_qc/$sn/$sn.sample.Rfam_table.xls -sumary $outdir/fastq_qc/$sn/$sn.sample.Rfam_summary.xls -maxEvalue 0.01");
		
		$queue->addcommond("gzip $outdir/fastq_qc/$sn/$sn.trim.R1.fq");
		$queue->addcommond("gzip $outdir/fastq_qc/$sn/$sn.trim.R2.fq");
		
		$queue->addcommond("# rm $outdir/fastq_qc/$sn/$sn.R1.fq $outdir/fastq_qc/$sn/$sn.R2.fq");
		
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("所有样本原始数据质量评估整合");
	$queue->set_work_dir($workdir);

	my $cmd;
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/fastq_qc/rawdata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_qc/$sn/$sn.rawdata.stat.xls";
	}
	$queue->addcommond("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/fastq_qc/trimdata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_qc/$sn/$sn.trimdata.stat.xls";
	}
	$queue->addcommond("$cmd");
	
	$queue->addcommond("grep '^rRNA' $outdir/fastq_qc/*/*.sample.Rfam_summary.xls | sed 's#$outdir/fastq_qc/##g' | cut -f 2 -d '/' | sed 's/.sample.Rfam_summary.xls:rRNA//g' | cut -f 1,3 | sed '1 isample\trRNA_percentage' > $outdir/fastq_qc/rRNA_percentage.xls");
	
	$queue->addcommond("# rm $outdir/fastq_qc/*.sample.fq $outdir/fastq_qc/*.rawdata.stat $outdir/fastq_qc/*.rawdata.stat.xls $outdir/fastq_qc/*.trimdata.stat.xls");
	
	$queue->run();
	$queue->wait();
}

sub run_rockhopper{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("rockhopper 索引构建");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/rockhopper/index"); 
	$queue->addcommond("module load bioawk/1.0");
	$queue->addcommond("module load samtools/1.6");
	$queue->addcommond("ln -sf $refer $outdir/rockhopper/index/genome.fa");
	$queue->addcommond("samtools faidx $outdir/rockhopper/index/genome.fa");
	$queue->addcommond("bioawk  -c fastx '{print \$name\"\\t\"length(\$seq)}' $outdir/rockhopper/index/genome.fa > $outdir/rockhopper/index/genome.length");
	$queue->addcommond("$script/gtf2ptt_rockhopper.pl -i $gtf  -l $outdir/rockhopper/index/genome.length -ref $outdir/rockhopper/index/genome.fa -o $outdir/rockhopper/index 2> $outdir/rockhopper/index/rock_index.log");
	$queue->run();
	$queue->wait();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("rockhopper 分析运行");
	$queue->set_work_dir($workdir);
	$queue->addcommond("genome_names=\$(dirname \$(ls $outdir/rockhopper/index/*/*.fna) | sed ':a;N;s/\\n/,/g;ta')");
	$queue->addcommond("mkdir -p $outdir/rockhopper/analysis");
	my $cmd = "java -Xmx24000m -cp /home/pub/software/Rockhopper/2.03/Rockhopper.jar Rockhopper -p 20 -z 0.2 -g \$genome_names -o $outdir/rockhopper/analysis";
	my @label;
	foreach my $sn(sort keys %samples){
		push(@label, $sn);
		$cmd .= " $outdir/fastq_qc/$sn/$sn.trim.R1.fq.gz"."%"."$outdir/fastq_qc/$sn/$sn.trim.R2.fq.gz";
	}
	$cmd .= " -L ".join(",", @label)." -v true -SAM";
	if($strand eq "yes"){
		$cmd .= " -s true -rf";
	}
	$queue->addcommond("$cmd");
	$queue->addcommond("# rm $outdir/rockhopper/analysis/*.sam");
	$queue->run();
	$queue->wait();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("rockhopper 结果整理");
	$queue->set_work_dir($workdir);
	$queue->addcommond("module load app/highcharts_tools");
	$queue->addcommond("module load bedtools/2.25.0");
	$queue->addcommond("module load app/universe_tools");
	$queue->addcommond("module load samtools/1.6");
	$queue->addcommond("mkdir -p $outdir/rockhopper/result");
	#$queue->addcommond("genomelist=\$(dirname \$(ls $outdir/rockhopper/index/*/*.fna) | sed ':a;N;s/\\n/,/g')");
	$queue->addcommond("genomelist=\$(dirname \$(cd $outdir/rockhopper/index/ && ls */*.fna) | sed ':a;N;s/\\n/,/g;ta')");
	$queue->addcommond("transcripts=\$(ls $outdir/rockhopper/analysis/*_transcripts.txt | sed ':a;N;s/\\n/,/g;ta')");
	$queue->addcommond("operons=\$(ls $outdir/rockhopper/analysis/*_operons.txt | sed ':a;N;s/\\n/,/g;ta')");
	$queue->addcommond("cd $outdir/rockhopper/result");
	$queue->addcommond("$script/rockhopper2bed.pl -t \$transcripts -seq \$genomelist -operon \$operons");
	$queue->addcommond("Rscript $script/operons_plot.R operon.xls");
	
	$queue->addcommond("awk '{if((\$3-\$2) >= 50 && (\$3-\$2) <= 500 ){print \$0}}' genome.predicted_RNA.bed > genome.sRNA.bed");
	$queue->addcommond("awk '{if((\$3-\$2) <  50 && (\$3-\$2) >  500 ){print \$0}}' genome.predicted_RNA.bed > genome.other_predicted_RNA.bed");
	$queue->addcommond("$script/bed2gtf_sRNA.pl -input $outdir/rockhopper/result/genome.predicted_RNA.bed -output $outdir/rockhopper/result/genome.predicted_RNA.gtf -str  predicted_RNA");
	$queue->addcommond("$script/bed2gtf_sRNA.pl -input $outdir/rockhopper/result/genome.sRNA.bed -output $outdir/rockhopper/result/genome.sRNA.gtf -str sRNA");
	$queue->addcommond("$script/bed2gtf_sRNA.pl -input $outdir/rockhopper/result/genome.other_predicted_RNA.bed -output $outdir/rockhopper/result/genome.other_predicted_RNA.gtf -str other_predicted_RNA");

	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed genome.feature.bed -s -name -fo genome.feature.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed genome.known_ncRNA.bed -s -name -fo genome.known_ncRNA.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed genome.mRNA.bed -s -name -fo genome.mRNA.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed genome.predicted_RNA.bed -s -name -fo genome.predicted_RNA.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed genome.sRNA.bed -s -name -fo genome.sRNA.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed UTR5.bed -s -name -fo UTR5.fa");
	$queue->addcommond("bedtools getfasta -fi $outdir/rockhopper/index/genome.fa -bed UTR3.bed -s -name -fo UTR3.fa");
	
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' genome.feature.fasta > genome.feature.fa");
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' genome.known_ncRNA.fasta > genome.known_ncRNA.fa");
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' genome.mRNA.fasta > genome.mRNA.fa");
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' genome.predicted_RNA.fasta > genome.predicted_RNA.fa");
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' UTR5.fasta > UTR5.fa");
	# $queue->addcommond("awk '{if(/^>/){if(match(\$0,\"^(>.+)::\",m)){print m[1]}else{print \$0}}else{print \$0}}' UTR3.fasta > UTR3.fa");
	
	$queue->addcommond("$script/FastaStats.pl genome.feature.fa &");
	$queue->addcommond("$script/FastaStats.pl genome.known_ncRNA.fa &");
	$queue->addcommond("$script/FastaStats.pl genome.mRNA.fa &");
	$queue->addcommond("$script/FastaStats.pl genome.sRNA.fa &");
	$queue->addcommond("$script/FastaStats.pl genome.predicted_RNA.fa &");
	$queue->addcommond("$script/FastaStats.pl UTR5.fa &");
	$queue->addcommond("$script/FastaStats.pl UTR3.fa &");
	$queue->addcommond("wait");
	$queue->addcommond("cd -");
	
	
	#$queue->addcommond("samtools view -bt rockhopper.genome.fa.fai -o $result/".$line[0].".bam $result/".$rock_sam.".sam && samtools sort $result/".$line[0].".bam $result/".$line[0].".bam.sorted && samtools index $result/".$line[0].".bam.sorted.bam");
	
	$queue->run();
	$queue->wait();
	
}

sub gene_annotation{
	$queue->set_job_cpu(40);
	$queue->set_job_mem(40);
	$queue->set_job_local(1);
	$queue->set_job_name("基因注释(mRNA)");
	$queue->set_work_dir($workdir);
	$queue->addcommond("module load app/annotation_tools");
	$queue->addcommond("module load app/annotation_meta");
	$queue->addcommond("module load EMBOSS/6.6.0");
	$queue->addcommond("module load diamond/v0.9.19.120");
	$queue->addcommond("module load blast+/2.3.0");
	$queue->addcommond("module load Prodigal/v2.6.3");
	$queue->addcommond("module load hmmer/3.1b2");
	
	$queue->addcommond("mkdir -p $outdir/gene_annotation");
	
	$queue->addcommond("mkdir -p $outdir/gene_annotation/basic");
	$queue->addcommond("cd $outdir/gene_annotation/basic");
	$queue->addcommond("annotation.pl -fasta $outdir/rockhopper/result/genome.mRNA.fa -chunksize 600 -nr /home/pub/database/blastdb/diamond/nr_bacteria -swissprot /home/pub/database/blastdb/diamond/swissprot -string /home/pub/database/blastdb/diamond/string -kegg /home/pub/database/blastdb/diamond/kegg_bacteria -job_prefix annot -max_cpus 200 -max_jobs 15");
	$queue->addcommond("cd -");
	
	$queue->addcommond("mkdir -p $outdir/gene_annotation/card");
	$queue->addcommond("cd $outdir/gene_annotation/card");
	$queue->addcommond("annot-card-meta-v1.0.pl -fasta $outdir/rockhopper/result/genome.mRNA.fa -ftype nuc -chunksize 600");
	$queue->addcommond("cd -");
	
	$queue->addcommond("mkdir -p $outdir/gene_annotation/cazy");
	$queue->addcommond("cd $outdir/gene_annotation/cazy");
	$queue->addcommond("transeq -sequence $outdir/rockhopper/result/genome.mRNA.fa -outseq protein.fa -frame 1 -table 11");
	$queue->addcommond("sed -i 's/_1\$//g' protein.fa");
	$queue->addcommond("annot-cazy-meta-v1.0.pl -fasta protein.fa -chunksize 600");
	$queue->addcommond("cd -");
	
	$queue->run();
}

sub sRNA_analysis{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_local($job_local);
	$queue->set_job_name("sRNA注释");
	$queue->set_work_dir($workdir);
	$queue->addcommond("module load app/annotation_tools");
	$queue->addcommond("module load blast+/2.3.0");
	$queue->addcommond("mkdir -p $outdir/sRNA_analysis/annotation");
	$queue->addcommond("cd $outdir/sRNA_analysis/annotation");
	$queue->addcommond("$script/sRNA_annot.sh $outdir/rockhopper/result/genome.sRNA.fa sRNA 0.001");
	$queue->addcommond("cd -");
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("sRNA靶基因预测");
	$queue->set_work_dir($workdir);
	$queue->addcommond("module load app/universe_tools");
	$queue->addcommond("module load app/venn_tools");
	$queue->addcommond("module load bedtools/2.25.0");
	#$queue->addcommond("module load RIsearch/2.1");
	#$queue->addcommond("module load RNAhybrid/2.1.2");
	#$queue->addcommond("module load ViennaRNA/2.4.12");
	$queue->addcommond("mkdir -p $outdir/sRNA_analysis/target");
	$queue->addcommond("cd $outdir/sRNA_analysis/target");
	$queue->addcommond("awk '{if(\$7-35<0){\$7=35};if(\$8-20<0){\$8=20}; if(\$6==\"+\"){print \$1\"\\t\"\$7-35\"\\t\"\$7+21\"\\t\"\$4\"\\t0\\t+\"}else{print \$1\"\\t\"\$8-20\"\\t\"\$8+36\"\\t\"\$4\"\\t0\\t-\"}}' $outdir/rockhopper/result/genome.mRNA.bed > genome.tar.bed");
	$queue->addcommond("fastalength $outdir/rockhopper/index/genome.fa | awk '{print \$2\"\\t0\\t\"\$1}' > genome.bed");
	$queue->addcommond("bedtools  intersect -a genome.tar.bed -b genome.bed  > genome.tar.range.bed");
	$queue->addcommond("bedtools  getfasta -name -fi $outdir/rockhopper/index/genome.fa -s -bed genome.tar.range.bed -fo genome.tar.fa");
	$queue->addcommond("$script/sRNA_target.pl -query $outdir/rockhopper/result/genome.sRNA.fa -target genome.tar.fa -thread 20 -method RIsearch,RNAhybrid -outdir ./");
	$queue->addcommond("awk '{if(NR==1){print \$0} else if(\$5 < -50 && \$6 < 0.01){print \$0}}' RNAhybrid.xls > RNAhybrid.clean.xls");
	$queue->addcommond("sed 1d RNAhybrid.clean.xls | awk '{print \$1\"\&\&\"\$2}' | sort -u > RNAhybrid.clean.list");
	$queue->addcommond("echo -e \"query_id\\tquery_start\\tquery_end\\ttarget_id\\ttarget_start\\ttarget_end\\tstrand\\tenergy\\tdetail\" > risearch.head");
	$queue->addcommond("awk '{if(\$8 < -50){print \$0}}' risearch.xls > risearch.clean.tmp");
	$queue->addcommond("cat risearch.head risearch.clean.tmp > RIsearch.clean.xls");
	$queue->addcommond("sed 1d RIsearch.clean.xls | awk '{print \$1\"\&\&\"\$4}' | sort -u > RIsearch.clean.list");
	$queue->addcommond("venn.pl -f RNAhybrid.clean.list,RIsearch.clean.list -l RNAhybrid,RIsearch -o Venn");
	$queue->addcommond("cd -");
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("sRNA二级结构预测");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/sRNA_analysis/structure");
	$queue->addcommond("cd $outdir/sRNA_analysis/structure");
	$queue->addcommond("$script/sRNA_structure.pl -query $outdir/rockhopper/result/genome.sRNA.fa -thread 20 -outdir ./");
	$queue->addcommond("cd -");
	$queue->run();
	
	$queue->wait();
}

sub refer_preparation{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("参考基因组格式化");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/refer");
	
	$queue->addcommond("module load rsem/1.3.1");
	$queue->addcommond("module load hisat/2.1.0");
	$queue->addcommond("module load star/2.5.4b");
	$queue->addcommond("module load rsem/1.3.1");
	$queue->addcommond("module load salmon/0.13.1");
	$queue->addcommond("module load kallisto/0.43.1");
	
	$queue->addcommond("ln -sf $refer $outdir/refer/genome.fa");
	$queue->addcommond("cat $gtf $outdir/rockhopper/result/genome.sRNA.gtf $outdir/rockhopper/result/genome.other_predicted_RNA.gtf | sort -k 1,1 -k 4,4n -k 5,5n -k 3,3 -u > $outdir/refer/gene.gtf");
	$queue->addcommond("ln -s $outdir/rockhopper/result/genome.feature.fa $outdir/refer/transcript.fa");
	
	$queue->addcommond("$script/gtf2bed.pl --input $outdir/refer/gene.gtf --output $outdir/refer/gene.bed");
	
	#$queue->addcommond("$script/gtf_g2t.pl -gtf $outdir/refer/gene.gtf  -out $outdir/refer/gene.g2t");
	#$queue->addcommond("$script/gtf_t2g2n.pl -gtf $outdir/refer/gene.gtf   -out $outdir/refer/gene.t2g2n");
	$queue->addcommond("$script/gtf_genename.pl -gtf $outdir/refer/gene.gtf -out $outdir/refer/gene_name.txt");
	$queue->addcommond("$script/gtf_genetype.pl -gtf $outdir/refer/gene.gtf -out $outdir/refer/gene_type.txt");
	#$queue->addcommond("cut -f 4 $outdir/rockhopper/result/genome.sRNA.bed | paste - - > $outdir/refer/gene_name.2");
	#$queue->addcommond("cut -f 4 $outdir/rockhopper/result/genome.sRNA.bed | awk '{print \$1\"\tsRNA\"}' > $outdir/refer/gene_type.2");
	#$queue->addcommond("cut -f 4 $outdir/rockhopper/result/genome.other_predicted_RNA.bed | paste - - > $outdir/refer/gene_name.3");
	#$queue->addcommond("cut -f 4 $outdir/rockhopper/result/genome.other_predicted_RNA.bed | awk '{print \$1\"\tother_predicted_RNA\"}' > $outdir/refer/gene_type.3");
	#$queue->addcommond("cat $outdir/refer/gene_name.1 $outdir/refer/gene_name.2 $outdir/refer/gene_name.3 > $outdir/refer/gene_name.txt");
	#$queue->addcommond("cat $outdir/refer/gene_type.1 $outdir/refer/gene_type.2 $outdir/refer/gene_name.3 > $outdir/refer/gene_type.txt");
	$queue->addcommond("cut -f 1 $outdir/refer/gene_name.txt | sed 1d | sort -u > $outdir/refer/gene_all.list");
	
	$queue->addcommond("mkdir -p $outdir/refer/genome_index_star");
	$queue->addcommond("cd $outdir/refer/");
	$queue->addcommond("STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $outdir/refer/genome_index_star --genomeFastaFiles $outdir/refer/genome.fa --sjdbGTFfile $outdir/refer/gene.gtf --genomeSAindexNbases 10 &> $outdir/refer/genome_index_star.log");
	$queue->addcommond("cd -");
	$queue->addcommond("hisat2-build $outdir/refer/genome.fa $outdir/refer/genome_index_hisat2");
	$queue->addcommond("bowtie2-build $outdir/refer/genome.fa $outdir/refer/genome_index_bowtie2");
	
	
	$queue->addcommond("rsem-prepare-reference --bowtie2 $outdir/refer/transcript.fa $outdir/refer/transcript_bowtie2_index");
	$queue->addcommond("salmon index -t $outdir/refer/transcript.fa -i $outdir/refer/transcript_salmon_index --type quasi -k 31");
	$queue->addcommond("kallisto index $outdir/refer/transcript.fa -i $outdir/refer/transcript_kallisto_index");
	
	#$queue->addcommond("rsem-prepare-reference --gtf $outdir/refer/gene.gtf --transcript-to-gene-map $outdir/refer/gene.g2t $outdir/refer/genome.fa $outdir/refer/genome_index_rsem");
	#$queue->addcommond("ln -s $outdir/refer/genome_index_rsem.transcripts.fa $outdir/refer/transcript.fa");
	# $queue->addcommond("rsem-prepare-reference --bowtie2 --gtf $outdir/refer/gene.gtf --transcript-to-gene-map $outdir/refer/gene.g2t $outdir/refer/genome.fa $outdir/refer/rsem_genome");
	# $queue->addcommond("bowtie2-build $- outdir/refer/rsem_genome.transcripts.fa $outdir/refer/bowtie2_transcript");
	
	$queue->run();
	$queue->wait();

}

sub alignment {
	foreach my $sample(sort keys %samples){
		
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_name("样本[$sample]比对基因组");
		$queue->set_work_dir($workdir);
 
		$queue->addcommond("module load hisat/2.1.0"); 
		$queue->addcommond("module load star/2.5.4b"); 
		$queue->addcommond("module load samtools/1.6"); 
		$queue->addcommond("mkdir -p $outdir/alignment");
		$queue->addcommond("cd $outdir/alignment");
		#$queue->addcommond("STAR --genomeDir $outdir/refer/genome_index_star --sjdbGTFfile $outdir/refer/gene.gtf --limitBAMsortRAM 40000000000 --runThreadN 20 --limitIObufferSize 500000000 --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --outFilterMultimapNmax 20 --outFilterMatchNminOverLread 0.66 --outFilterIntronMotifs None --outSJfilterReads All --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMstrandField intronMotif --outSAMattrRGline ID:$sample SM:$sample PL:illumina --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 0 --chimScoreDropMax 20 --chimScoreSeparation 10 --chimScoreJunctionNonGTAG -1 --quantMode TranscriptomeSAM --quantTranscriptomeBan IndelSoftclipSingleend --outReadsUnmapped Fastx --readFilesIn $outdir/fastq_qc/$sample.trim.R1.fq.gz $outdir/fastq_qc/$sample.trim.R2.fq.gz --readFilesCommand zcat --outFileNamePrefix  $sample.");
		#$queue->addcommond("mv $outdir/alignment/$sample.Aligned.sortedByCoord.out.bam $sample.sorted.bam");
		#$queue->addcommond("samtools index $sample.sorted.bam");
		
		if($strand eq "yes"){
			$queue->addcommond("hisat2 -q -x $outdir/refer/genome_index_hisat2 -1 $outdir/fastq_qc/$sample/$sample.trim.R1.fq.gz -2 $outdir/fastq_qc/$sample/$sample.trim.R2.fq.gz -p 10 -S $outdir/alignment/$sample.sam --dta --rna-strandness RF > $outdir/alignment/$sample.stdout 2> $outdir/alignment/$sample.stderr");
		}else{
			$queue->addcommond("hisat2 -q -x $outdir/refer/genome_index_hisat2 -1 $outdir/fastq_qc/$sample/$sample.trim.R1.fq.gz -2 $outdir/fastq_qc/$sample/$sample.trim.R2.fq.gz -p 10 -S $outdir/alignment/$sample.sam --dta > $outdir/alignment/$sample.stdout 2> $outdir/alignment/$sample.stderr");
		}
		$queue->addcommond("samtools view -S -b $outdir/alignment/$sample.sam | samtools sort - -\@ 10 -o $outdir/alignment/$sample.bam");
		$queue->addcommond("samtools index $outdir/alignment/$sample.bam");
		
		$queue->addcommond("mkdir -p $outdir/alignment/tmp");
		$queue->addcommond("java -jar -Xmx48g -Djava.io.tmpdir=$outdir/alignment/tmp -jar /home/pub/software/picard-tools/2.18.12/picard.jar CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT INPUT=$sample.bam OUTPUT=$sample.alignment_metrics.tsv REFERENCE_SEQUENCE=$outdir/refer/genome.fa");
		
		$queue->addcommond("java -Xmx8g -Djava.io.tmpdir=$outdir/alignment/tmp -jar /home/pub/software/picard-tools/2.18.12/picard.jar CollectInsertSizeMetrics I=$sample.bam O=$sample.insertSize.txt H=$sample.insertSize.pdf");

		#$queue->addcommond("java -Xmx8g -Djava.io.tmpdir=$outdir/alignment/tmp -jar /home/pub/software/picard-tools/2.18.12/picard.jar CollectRnaSeqMetrics VALIDATION_STRINGENCY=SILENT REF_FLAT=$outdir/refer/genome.fa RIBOSOMAL_INTERVALS=${ribsomal} STRAND_SPECIFICITY=${strand} INPUT=${prefix}_sorted.bam CHART_OUTPUT=${prefix}_3bias.pdf OUTPUT=${prefix}_RNA_metrics.tsv REFERENCE_SEQUENCE=$outdir/refer/genome.fa");
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(5);
	$queue->set_job_name("比对基因组结果统计分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("$script/alignment_hisat2.stat.pl $outdir/alignment/ $outdir/alignment/alignment.statistics.xls");
	
	$queue->run();
	$queue->wait();
}

sub explevel {
	foreach my $sample(sort keys %samples){
		$queue->set_job_cpu(8);
		$queue->set_job_mem(10);
		$queue->set_job_name("样本 $sample 表达量评估");
		$queue->addcommond("module load salmon/0.13.1");
		$queue->addcommond("module load rsem/1.3.1");
		$queue->addcommond("module load bowtie/2.3.3");
		$queue->addcommond("mkdir -p $outdir/explevel/salmon/");
		
		$queue->addcommond("salmon quant --index $outdir/refer/transcript_salmon_index --libType A --mates1 $outdir/fastq_qc/$sample/$sample.trim.R1.fq.gz --mates2 $outdir/fastq_qc/$sample/$sample.trim.R2.fq.gz -o $outdir/explevel/salmon/$sample");
		
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(5);
	$queue->set_job_name("生成所有样本的表达量表格");
	my $cmd = "$script/explevel_matrix_salmon.pl --g2n $outdir/refer/gene_name.txt --g2t $outdir/refer/gene_type.txt --prefix $outdir/explevel/matrix_gene";
	foreach my $sample(sort keys %samples){
		$cmd .= " $outdir/explevel/salmon/$sample/quant.sf";
	}
	$queue->addcommond("$cmd");
	$queue->addcommond("cut -f 1,2 $outdir/explevel/matrix_gene.count.xls > $outdir/explevel/gene.id2name.txt");
	$queue->run();
	$queue->wait();
}

sub sample_cluster{
	$queue->set_job_cpu(2);
	$queue->set_job_mem(5);
	$queue->set_job_name("基于表达的样本聚类分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/sample_cluster");
	$queue->addcommond("module load app/cluster_tools");
	$queue->addcommond("plot_pca_exp.pl -i $outdir/explevel/matrix_gene.tpm.txt -o $outdir/sample_cluster -pre Sample");
	foreach my $class(sort keys %meta){
		next if $class eq "Sample";
		$queue->addcommond("plot_pca_exp.pl -i $outdir/explevel/matrix_gene.tpm.txt -o $outdir/sample_cluster -m $file_meta -g $class -pre $class");
	}
	$queue->run();
	$queue->wait();
}

sub expdiff_analysis{
	foreach my $std(sort keys %group_info){
		my $class = [split /\./,$std]->[0];
		
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_name("差异表达基因鉴定与分析[$std]");
		$queue->set_work_dir($workdir);
		$queue->addcommond("module load app/expdiff_tools");
		$queue->addcommond("mkdir -p $outdir/expdiff_analysis");
		$queue->addcommond("expdiff_analysis.pl --count $outdir/explevel/matrix_gene.count.txt --norm $outdir/explevel/matrix_gene.tpm.txt --method edgeR --pval $de_pval --fdr $de_fdr --log2FC $de_log2fc --num_parallel 6 --group $outdir/group_info/$std.txt --outdir $outdir/expdiff_analysis/$std --prefix $class --id2name $outdir/explevel/gene.id2name.txt --go_asso $outdir/gene_annotation/basic/annotation/GO/transcript.GO.txt --kegg_asso $outdir/gene_annotation/basic/annotation/KEGG/transcript.KEGG.txt --species $species --class tpm"); #--reverse 
		$queue->run();
	}
	$queue->wait();
}

sub promoter{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->addcommond("module load app/sequence_tools");
	$queue->addcommond("module load bedtools/2.25.0");
	$queue->addcommond("module load bTSSfinder/1.2016");
	$queue->set_job_name("基因结构分析(promoter)");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/promoter");
	$queue->addcommond("cut -f 1-6 $outdir/rockhopper/result/genome.feature.bed | bedtools flank -i - -g $outdir/rockhopper/index/genome.length -l 1000 -r 0 -s > $outdir/promoter/gene.upstream.bed");
	$queue->addcommond("bedtools getfasta -fi $outdir/refer/genome.fa -bed $outdir/promoter/gene.upstream.bed -s -name -fo $outdir/promoter/gene.upstream.fasta");
	$queue->addcommond("reforamt_fasta_linewidth.pl $outdir/promoter/gene.upstream.fasta $outdir/promoter/gene.upstream.fa");
	$queue->addcommond("cd $outdir/promoter");
	$queue->addcommond("bTSSfinder -i gene.upstream.fa -o promoter -h 1");
	$queue->addcommond("cd -");
	$queue->run();
}

sub snp_indel{
	foreach my $sample(sort keys %samples){
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->addcommond("module load app/snp_tools");
		$queue->set_job_name("样本[$sample]的SNP/Indel分析");
		$queue->set_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/snp_indel");
		$queue->addcommond("rna_snp_GATK.pl --fq1 $outdir/fastq_qc/$sample/$sample.trim.R1.fq.gz --fq2 $outdir/fastq_qc/$sample/$sample.trim.R2.fq.gz --refer $outdir/refer/genome.fa --gtf $outdir/refer/gene.gtf --index $outdir/refer/genome_index_star --outdir $outdir/snp_indel --prefix $sample");
		$queue->run();
	}
	$queue->wait();
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->addcommond("module load app/snp_tools");
	$queue->set_job_name("所有样本的SNP/Indel统计分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("cd $outdir/snp_indel");
	$queue->addcommond("ls *.anno_snp.avinput | sed 's/.anno_snp.avinput//g' > sample.list");
	$queue->addcommond("snp_annovar_stat.pl --input sample.list --name $outdir/refer/gene_name.txt --outdir stat");
	$queue->run();
}

sub transcriptome_assesement{
	foreach my $sample(sort keys %samples){
		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_work_dir($workdir);
		$queue->set_job_name("转录组评估之 $sample 测序饱和度分析");
		$queue->addcommond("mkdir -p $outdir/transcriptome_assesement/$sample");
		$queue->addcommond("RPKM_saturation.py -i $outdir/alignment/$sample.bam -r $outdir/refer/gene.bed -o $outdir/transcriptome_assesement/$sample/$sample");
		$queue->addcommond("$script/saturation2plot.pl -in $outdir/transcriptome_assesement/$sample/$sample.eRPKM.xls -out $outdir/transcriptome_assesement/$sample/$sample");
		$queue->run();

		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_work_dir($workdir);
		$queue->set_job_name("转录组评估之 $sample 基因body分布评估");
		$queue->addcommond("mkdir -p $outdir/transcriptome_assesement/$sample");
		$queue->addcommond("head -n 1000 $outdir/refer/gene.bed > $outdir/refer/gene.tmp");
		$queue->addcommond("geneBody_coverage.py -i $outdir/alignment/$sample.bam  -r $outdir/refer/gene.tmp -o $outdir/transcriptome_assesement/$sample/$sample");
		$queue->run();

		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_work_dir($workdir);
		$queue->set_job_name("转录组评估之 $sample read重复度评估");
		$queue->addcommond("mkdir -p $outdir/transcriptome_assesement/$sample");
		$queue->addcommond("read_duplication.py -i $outdir/alignment/$sample.bam -o $outdir/transcriptome_assesement/$sample/$sample");
		$queue->run();
		
		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_work_dir($workdir);
		$queue->set_job_name("转录组评估之 $sample insert size评估");
		$queue->addcommond("inner_distance.py -i $outdir/alignment/$sample.bam -r $outdir/refer/gene.bed -o $outdir/transcriptome_assesement/$sample/$sample");
		$queue->run();
		
		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_work_dir($workdir);
		$queue->set_job_name("转录组评估之 $sample 链特异性建库评估");
		$queue->addcommond("infer_experiment.py -i $outdir/alignment/$sample.bam -r $outdir/refer/gene.bed > $outdir/transcriptome_assesement/$sample/$sample.infer_experiment.txt");
		$queue->run();
	}
	$queue->wait();
}

sub result_collection{
	$queue->set_job_cpu(2);
	$queue->set_job_mem(5);
	$queue->set_work_dir($workdir);
	$queue->set_job_name("结果整理");
	
	my ($num, $id, $dir);
	$num = 0;
	
	$queue->addcommond("mkdir -p $outdir/Result");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Rawdata_Stat";
	$queue->addcommond("mkdir -p $dir $dir/quality_control");
	$queue->addcommond("cp $outdir/fastq_qc/trimdata.stat.xls   $dir/");
	$queue->addcommond("cp $outdir/fastq_qc/rawdata.stat.xls    $dir/");
	$queue->addcommond("cp $outdir/fastq_qc/rRNA_percentage.xls $dir/");
	foreach my $sn ( sort keys %samples ) {
		foreach my $x( qw/R1 R2/ ){
			$queue->addcommond("cp -r $outdir/fastq_qc/$sn/$sn.$x\_fastqc $dir/quality_control/$sn.raw.$x");
			$queue->addcommond("cp -r $outdir/fastq_qc/$sn/$sn.trim.$x\_fastqc $dir/quality_control/$sn.trim.$x");
			$queue->addcommond("cp $outdir/fastq_qc/$sn/$sn.$x\_fastqc.html $dir/quality_control/$sn.raw.$x.html");
			$queue->addcommond("cp $outdir/fastq_qc/$sn/$sn.trim.$x\_fastqc.html $dir/quality_control/$sn.trim.$x.html");
		}
	}
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Alignment";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/alignment/alignment.statistics.xls $dir/");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Quality_Assessment";
	$queue->addcommond("mkdir -p $dir");
	foreach my $sn ( sort keys %samples ) {
		$queue->addcommond("mkdir -p $dir/$sn");
		$queue->addcommond("cp $outdir/transcriptome_assesement/$sn/$sn.DupRate_plot.pdf $dir/$sn");
		$queue->addcommond("cp $outdir/transcriptome_assesement/$sn/$sn.inner_distance_plot.pdf $dir/$sn");
		$queue->addcommond("cp $outdir/transcriptome_assesement/$sn/$sn.saturation.pdf $dir/$sn");
		$queue->addcommond("cp $outdir/transcriptome_assesement/$sn/$sn.geneBodyCoverage.curves.pdf $dir/$sn");
		$queue->addcommond("cp $outdir/transcriptome_assesement/$sn/$sn.infer_experiment.txt $dir/$sn");
	}
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.AnnotationStat";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/ORF/ $dir/ORF");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/NR/  $dir/NR");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/GO/  $dir/GO");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/KEGG/  $dir/KEGG");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/COG_KOG $dir/COG_KOG");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/annotstat.xls $dir/");
	$queue->addcommond("cp -r $outdir/gene_annotation/basic/annotation/transcript.annotation.table.xls $dir/");
	$queue->addcommond("mkdir -p $dir/CARD");
	$queue->addcommond("cp -r $outdir/gene_annotation/card/outdir/gene.card.annot.xls $outdir/gene_annotation/card/outdir/gene.card.blast.xls $dir/CARD");
	$queue->addcommond("mkdir -p $dir/CAZY");
	$queue->addcommond("cp -r $outdir/gene_annotation/cazy/outdir/gene.cazy.anno.xls $outdir/gene_annotation/cazy/outdir/gene.cazy.class.stat.xls $outdir/gene_annotation/cazy/outdir/gene.cazy.family.stat.xls $dir/CAZY");

	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.sRNA_Analysis";
	$queue->addcommond("mkdir -p $dir/01.sRNA_Prediction");
	$queue->addcommond("mkdir -p $dir/02.sRNA_Annotation");
	$queue->addcommond("mkdir -p $dir/03.sRNA_Secondstructure");
	$queue->addcommond("mkdir -p $dir/04.sRNA_Target");
	
	$queue->addcommond("cp $outdir/rockhopper/result/genome.predicted_RNA.bed $dir/01.sRNA_Prediction/genome.predicted_RNA.bed");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.predicted_RNA.fa  $dir/01.sRNA_Prediction/genome.predicted_RNA.fa ");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.predicted_RNA.fa.length.distribut.txt.pdf $dir/01.sRNA_Prediction/genome.predicted_RNA.length_distribut.pdf");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.predicted_RNA.gtf $dir/01.sRNA_Prediction/genome.predicted_RNA.gtf");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.sRNA.bed $dir/01.sRNA_Prediction/genome.sRNA.bed");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.sRNA.fa  $dir/01.sRNA_Prediction/genome.sRNA.fa ");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.sRNA.fa.length.distribut.txt.pdf $dir/01.sRNA_Prediction/genome.sRNA.length_distribut.pdf");
	$queue->addcommond("cp $outdir/rockhopper/result/genome.sRNA.gtf $dir/01.sRNA_Prediction/genome.sRNA.gtf");
	
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/All_rfam_stat.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/rfam.stat.pdf $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_annot_stat.pdf $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_annot_stat.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_rfam.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_SIPHI.annot.detail.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_SIPHI.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_sRNAMap.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_sRNATarBase.annot.detail.xls $dir/02.sRNA_Annotation");
	$queue->addcommond("cp $outdir/sRNA_analysis/annotation/sRNA_vs_sRNATarBase.xls $dir/02.sRNA_Annotation");
	
	$queue->addcommond("cp -r $outdir/sRNA_analysis/structure/RNAfold_pdf $dir/03.sRNA_Secondstructure");
	$queue->addcommond("cp -r $outdir/sRNA_analysis/structure/RNAfold.str.txt $dir/03.sRNA_Secondstructure");
	
	$queue->addcommond("cp -r $outdir/sRNA_analysis/target/RIsearch.clean.xls $dir/04.sRNA_Target/RIsearch.xls");
	$queue->addcommond("cp -r $outdir/sRNA_analysis/target/RNAhybrid.clean.xls $dir/04.sRNA_Target/RNAhybrid.xls");
	$queue->addcommond("cp -r $outdir/sRNA_analysis/target/Venn.pdf $dir/04.sRNA_Target");
	$queue->addcommond("cp -r $outdir/sRNA_analysis/target/Venn.xls $dir/04.sRNA_Target");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Expression_Analysis";
	$queue->addcommond("mkdir -p $dir/pca $dir/cor");
	$queue->addcommond("cp -r $outdir/explevel/matrix_gene.* $dir");
	$queue->addcommond("cp -r $outdir/sample_cluster/Sample.pca* $dir/pca");
	$queue->addcommond("cp -r $outdir/sample_cluster/Sample.cor* $dir/cor");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.DEG_Analysis";
	$queue->addcommond("mkdir -p $dir/01.DEG_Identification");
	$queue->addcommond("mkdir -p $dir/02.DEG_Visualisation");
	$queue->addcommond("mkdir -p $dir/03.DEG_Cluster");
	$queue->addcommond("mkdir -p $dir/04.DEG_Annotation/GO $dir/04.DEG_Annotation/KEGG");
	$queue->addcommond("mkdir -p $dir/05.DEG_Enrichment/GO $dir/05.DEG_Enrichment/KEGG");
	$queue->addcommond("mkdir -p $dir/06.DEG_IPATH");
	foreach my $std(sort keys %group_info){
		my $class = [split /\./,$std]->[0];
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.diffexpress_diff.xls $dir/01.DEG_Identification");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.diffexpress_total.xls $dir/01.DEG_Identification");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.volcano.pdf $dir/02.DEG_Visualisation");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.volcano.xls $dir/02.DEG_Visualisation");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.scatter.pdf $dir/02.DEG_Visualisation");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.scatter.xls $dir/02.DEG_Visualisation");
		$queue->addcommond("rm -f $dir/02.DEG_Visualisation/*.go_enrichment.scatter.pdf $dir/02.DEG_Visualisation/*.kegg_enrichment.scatter.pdf");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.heatmap.pdf $dir/03.DEG_Cluster");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.heatmap.xls $dir/03.DEG_Cluster");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.heatmap_subclusters.pdf $dir/03.DEG_Cluster");
		$queue->addcommond("cp -r $outdir/expdiff_analysis/$std/$class.*.heatmap_subclusters $dir/03.DEG_Cluster");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.up_down.level2-gobars.pdf $dir/04.DEG_Annotation/GO");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.up_down.level2-gobars.xls $dir/04.DEG_Annotation/GO");
		$queue->addcommond("cp -r $outdir/expdiff_analysis/$std/$class.*.kegg_annot $dir/04.DEG_Annotation/KEGG");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.go_enrichment.xls $dir/05.DEG_Enrichment/GO");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.go_enrichment.scatter.pdf $dir/05.DEG_Enrichment/GO");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.go_enrichment.bar.pdf $dir/05.DEG_Enrichment/GO");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.go_enrichment.relation.png $dir/05.DEG_Enrichment/GO");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.go_enrichment.relation.svg $dir/05.DEG_Enrichment/GO");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.kegg_enrichment.xls $dir/05.DEG_Enrichment/KEGG");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.kegg_enrichment.scatter.pdf $dir/05.DEG_Enrichment/KEGG");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.kegg_enrichment.bar.pdf $dir/05.DEG_Enrichment/KEGG");
		
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.antibiotic.pdf $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.metabolic.pdf $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.microbial.pdf $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.secondary.pdf $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.antibiotic.svg $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.metabolic.svg $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.microbial.svg $dir/06.DEG_IPATH");
		$queue->addcommond("cp $outdir/expdiff_analysis/$std/$class.*.ipath.secondary.svg $dir/06.DEG_IPATH");
	}
	
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Gene_Structure";
	$queue->addcommond("mkdir -p $dir/01.Operon");
	$queue->addcommond("mkdir -p $dir/02.TSS_and_TTS_prediction");
	$queue->addcommond("mkdir -p $dir/03.UTR");
	$queue->addcommond("cp $outdir/rockhopper/result/Fre_of_Operon_Length.pdf $dir/01.Operon");
	$queue->addcommond("cp $outdir/rockhopper/result/Fre_of_Operon_Size.pdf $dir/01.Operon");
	$queue->addcommond("cp $outdir/rockhopper/result/operon.xls $dir/01.Operon");
	
	$queue->addcommond("cp $outdir/rockhopper/result/TSS_and_TTS.xls $dir/02.TSS_and_TTS_prediction");
	
	$queue->addcommond("cp $outdir/rockhopper/result/TSS_and_TTS.xls $dir/03.UTR");
	$queue->addcommond("cp $outdir/rockhopper/result/UTR5.fa $dir/03.UTR");
	$queue->addcommond("cp $outdir/rockhopper/result/UTR3.fa $dir/03.UTR");
	$queue->addcommond("cp $outdir/rockhopper/result/UTR5.fa.length.distribut.txt.pdf $dir/03.UTR/UTR5.length_distribut.pdf");
	$queue->addcommond("cp $outdir/rockhopper/result/UTR3.fa.length.distribut.txt.pdf $dir/03.UTR/UTR3.length_distribut.pdf");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.SNP_Indel";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/snp_indel/stat/indel.* $dir/");
	$queue->addcommond("cp $outdir/snp_indel/stat/snp.*   $dir/");
	
	$queue->run();
	$queue->wait();
}

sub convert_num_to_id{
	my ($num) = @_;
	if($num < 10){
		return "0$num";
	}else{
		return "$num";
	}
}