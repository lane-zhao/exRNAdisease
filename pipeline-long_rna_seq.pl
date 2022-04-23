#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use List::Util qw/max min sum maxstr minstr shuffle/;
use BioQueue::Queue;

my $version  = "v2.0";
my $script   = "$Bin/script";
my $summary  = "$Bin/summary";
my $software = "$Bin/software";
my $database = "$Bin/database";
my $workdir  = getcwd;

my %species_info;
&load_species_info("$Bin/species.config.ini", \%species_info);

my ($help, $file_input, $outdir, $file_meta, $file_comp, $file_venn, $file_background, $strand, $report);
my ($species, $show, $refer_seq, $refer_gtf, $gene_go, $gene_kegg, $gene_annot, $tran_go, $tran_kegg, $tran_annot, $prot_seq, $prot_trans);
my ($index_hisat2, $index_bwa, $index_star, $index_salmon, $index_kallisto);
my ($method_exp, $de_method, $de_pval, $de_padj, $de_fc, $de_log2fc);
my ($class, $NR, $NT, $KEGG, $STRING, $SWISSPROT);
my ($kegg_species, $kegg_filter, $ppi_denovo, $ppi_taxonomy);
my ($run_trans, $run_index_hisat2, $run_index_star, $run_index_bwa, $run_index_kallisto, $run_index_salmon);
my ($run_novel, $run_snp, $run_as, $run_wgcna);
my ($annot_known, $no_ppi);
my (@exp_class, @exp_level);
my ($max_cpus, $max_mems, $max_jobs, $job_prefix, $job_mode, $queue_rerun, $queue_level, $partition);
my %bkg_info;
my (%meta_grp2sn,%meta_sn2grp);
my (%comp_grp2sn,%comp_sn2grp, %comp_param);
my (%sample_info);
my %all_dependencies;
my %results;

&GetOptions(
	"help!"           => \$help,
	"input:s"         => \$file_input,
	"meta:s"          => \$file_meta,
	"comp:s"          => \$file_comp,
	"venn:s"          => \$file_venn,
	"background:s"    => \$file_background,
	"outdir:s"        => \$outdir,
	"report!"         => \$report,
	
	"method_exp:s"    => \$method_exp,
	"de_method:s"     => \$de_method,
	"de_pval:f"       => \$de_pval,
	"de_padj:f"       => \$de_padj,
	"de_fc:f"         => \$de_fc,
	
	"species:s"       => \$species,
	"show!"           => \$show,
	
	"refer_seq:s"     => \$refer_seq,
	"refer_gtf:s"     => \$refer_gtf,
	"index_hisat2:s"  => \$index_hisat2,
	"index_star:s"    => \$index_star,
	"index_salmon:s"  => \$index_salmon,
	"index_kallisto:s"=> \$index_kallisto,
	
	"gene_go:s"       => \$gene_go,
	"gene_kegg:s"     => \$gene_kegg,
	"gene_annot:s"    => \$gene_annot,
	"tran_go:s"       => \$tran_go,
	"tran_kegg:s"     => \$tran_kegg,
	"tran_annot:s"    => \$tran_annot,
	"prot_seq:s"      => \$prot_seq,
	"prot_trans:s"    => \$prot_trans,
	
	"run_trans!"      => \$run_trans,
	"strand:s"        => \$strand,
	
	"run_novel!"      => \$run_novel,
	"run_snp!"        => \$run_snp,
	"run_as!"         => \$run_as,
	"run_wgcna!"      => \$run_wgcna,
	"no_ppi!"         => \$no_ppi,
	
	"class:s"         => \$class,
	"kegg_species:s"  => \$kegg_species,
	"kegg_filter!"    => \$kegg_filter,
	"ppi_denovo!"     => \$ppi_denovo,
	"ppi_taxonomy:i"  => \$ppi_taxonomy,
	
	"max_cpus:i"      => \$max_cpus,
	"max_mems:i"      => \$max_mems,
	"max_jobs:i"      => \$max_jobs,
	"partition:s"     => \$partition,
	"job_mode:s"      => \$job_mode,
	"job_prefix:s"    => \$job_prefix,
	"queue_level:s"   => \$queue_level,
	"rerun!"          => \$queue_rerun,
);

sub usage{
	my $program = basename($0);
	die "
Version  : $version
Usage    : $program [options]
Arguments:
  [Basic Options]:
    --help           <null>      print help information to screen
    --species        <string>    species name. [ optional ]
    --show           <null>      show the support species list
    
    --input          <string>    fastq file list. [ required ]
    --meta           <string>    meta information for grouping comparison. [ optional ]
    --comp           <string>    comparison configure file. [ optional ]
    --venn           <string>    venn configure file. [ optional ]
    --background     <string>    project background configure file. [ optional ]
    --outdir         <string>    output dir. default : \"out\".
    --strand         <string>    yes(strand-specific, RF, default) or no(non-strand-specific)
    --report         <null>      whether generate web report or not 
    
    --refer_seq      <string>    reference genome sequence in fasta fromat [ required ]
    --refer_gtf      <string>    reference gene annotation in gtf format [ required ]
    --index_hisat2   <string>    reference genome index in hisat2 format [ usefull if method_exp = A1, A2 OR A3 ]
    --index_star     <string>    reference genome index in star format [ usefull if method_exp = B1, B2 OR B3 ]
    --index_salmon   <string>    reference transcript index in salmon format [ usefull if method_exp = C1 ]
    --index_kallisto <string>    reference transcript index in salmon format [ usefull if method_exp = C2 ]
    
    --gene_annot     <string>    reference gene id to gene annotation [ optional ]
    --gene_go        <string>    reference gene id to go annotation file [ optional ]
    --gene_kegg      <string>    reference gene id to kegg annotation file [ optional ]
    --tran_annot     <string>    reference transcript id to transcript annotation [ optional ]
    --tran_go        <string>    reference transcript id to go annotation file [ optional ]
    --tran_kegg      <string>    reference transcript id to kegg annotation file [ optional ]
    --prot_seq       <string>    reference protein sequence file in fasta format [ optional ]
    --prot_trans     <string>    reference protein id to transcript id conversion file [ optional ]

  [Advanced Options]:
    --class          <string>    annotation class: animal,mammal,vertebrate,invertebrate,plant,fungi is available
    --method_exp     <string>    analysis method [alignment]->(assembly)->[explevel]
                                 default is A1
                                 A1 : [hisat2]->(stringtie)->[stringtie+ballgown]
                                 A2 : [hisat2]->(stringtie)->[featureCounts]
                                 A3 : [hisat2]->(stringtie)->[salmon]
                                 B1 : [star]->(stringtie)->[stringtie+ballgown]
                                 B2 : [star]->(stringtie)->[featureCounts]
                                 B2 : [star]->(stringtie)->[salmon]
                                 C1 : [salmon]
                                 C2 : [kallisto]
    --de_method     <string>    the default method of diff-exp analysis, default is edgeR
                                 available: edgeR | DESeq2 | DEGseq | limma | NOIseq | ttest | wilcox
                                 note: you should have biological replicates.
                                 a.edgeR only be used for comparison with/without biological replicates
                                 b.DESeq2 only be used for comparison with biological replicates
                                 with a fixed dispersion setting.
    --de_pval        <float>     the default threshold of de pval, default is 0
                                 if it set 0, this option will be missed
    --de_padj        <float>     the default threshold of de padj, default is 0.05
                                 if it set 0, this option will be missed
    --de_fc          <float>     the default threshold of fold change, default is 2
                                 if it set 0, this option will be missed

    --run_novel      <null>      identify novel transcripts or not
    --run_trans      <null>      run transcript analysis or not
    --run_as         <null>      run alternative splicing analysis or not
    --run_snp        <null>      run snp indel analysis or not
    --run_wgnca      <null>      run wgcna analysis or not
    

    --kegg_species   <string>    kegg species, default is ko
    --kegg_filter    <null>      kegg species are not hsa and mmu (filter HD & DD)
    --ppi_denovo     <null>      whether run ppi analysis in denovo mode or not.
    --ppi_taxonomy   <number>    taxonomy id of research species or similar species, 
                                 the id must exists in string db.

  [Pipeline Options]:
    --max_cpus       <number>    queue's max cpu number limitation, default is 500 [ optional ]
    --max_mems       <number>    queue's max mem size limitation[G], default is 1000 [optional]
    --max_jobs       <number>    queue's max job number limitation, default is 50  [ optional ]
    --partition      <string>    slurm or pbs partition, default is 'small'
    --queue_level    <string>    queue level : 'only_cmd', 'run_strict'. default is 'run_strict'
    --queue_rerun    <null>      rerun all pipeline
    --job_prefix     <string>    jjob name prefix, default is 'mrna' [optional]
    --job_mode       <string>    job run mode : 'local', 'slurm' and 'pbs'. default is pbs
    
###################
#PS1:<input> contain 3 colums, sample_id, read1.fq(gz), read2.fq(gz)
#    sample_1	/path/sample_1_R1.fq	/path/sample_1_R2.fq
#    sample_2	/path/sample_2_R1.fq	/path/sample_2_R2.fq
#    sample_3	/path/sample_3_R1.fq	/path/sample_3_R2.fq
#    ......
#    For different lines, samples name must be unique !!!
#
#PS2:[meta] contain at leat 2 colums, 1st is [Sample], 2nd, 3rd, 4th ... are [Group Standard]
#    #col_1      col_2        col_3
#    Sample      Group_STD_1  Group_STD_2
#    sample_1    g1           w
#    sample_2    g1           d
#    sample_3    g2           d
#    sample_4    g2           w
#    sample_5    g3           d
#    sample_6    g3           w
#    For different lines, samples name must be unique !!!
#
#PS3:[comp] 1st is Group Standard], 2nd is comparison, 3rd is expldiff method, 4th is expldiff pval, 5th is expldiff padj, 6th is expldiff fc
#    #STD         comparison          method     pval     padj     fc
#    Sample       sample_1,sample_3   edgeR      0        0.05    2
#    Group_STD_1  g1,g2               edgeR      0        0.05    2
#    Group_STD_1  g1,g3               deseq2     0        0.05    2
#    Group_STD_1  g2,g3               deseq2     0        0.01    1.5
#    Group_STD_2  w,d                 edgeR      0.05     0       1.5
#    For different lines, samples name must be unique !!!
#
#PS4:[venn] contain at leat 2 colums, 1st is [STD], 2nd is comparison, 3rd is expldiff method, 4th is expldiff pval, 5th is expldiff padj, 6th is expldiff fc
#    Sample:sample_1,sample_3    Group_STD_1:g1,g3    
#    Group_STD_1  g1,g3               edgeR      0        0.05    2
#    Group_STD_1  g1,g2,g3            deseq2     0        0.05    2
#    Group_STD_2  w,d                 edgeR      0.05     0       1.5
#    For different lines, samples name must be unique !!!
###################
\n";
}

&parsing_parameters();
my $queue = BioQueue::Queue->new(
	{
		'partition'   => $partition,
		'job_pref'    => $job_prefix,
		'job_mode'    => $job_mode,
		'max_cpus'    => $max_cpus,
		'max_jobs'    => $max_jobs,
		'max_mems'    => $max_mems,
		'queue_level' => $queue_level,
		'queue_rerun' => $queue_rerun,
	}
);

&parsing_sample_info;
&parsing_meta_and_comp_info;

&rawdata_quality_control();
&reference_preparation_known_1();
&genome_index();
&genome_alignment();
# &transcriptome_assesement();

# &transcript_assembly();
# &transcript_identification();
# &reference_preparation_total_1();
# &gene_annotation();
# &reference_preparation_known_2();
# &reference_preparation_total_2();

&explevel_analysis();
&expldiff_analysis();

# &wgcna_analysis();
# &snp_analysis();
# &splicing_analysis();

# &result_collection();
$queue->jointhreads();

# if($report){
	# &webpage_report();
# }

print " all jobs done!\n";

###########################################################################

sub parsing_parameters{
	&usage if ($help);
	
	unless(defined( $file_input )){
		print STDERR "Error:-input must be specified!\n";
		&usage;
	}
	
	if($species){
		if(exists $species_info{$species}){
			$refer_seq      = $species_info{$species}{refer_seq};
			$refer_gtf      = $species_info{$species}{refer_gtf};
			$index_hisat2   = $species_info{$species}{index_hisat2};
			$index_star     = $species_info{$species}{index_star};
			$index_salmon   = $species_info{$species}{index_salmon};
			$index_kallisto = $species_info{$species}{index_kallisto};
			
			$gene_annot     = $species_info{$species}{gene_annot};
			$gene_go        = $species_info{$species}{gene_go};
			$gene_kegg      = $species_info{$species}{gene_kegg};
			$tran_annot     = $species_info{$species}{tran_annot};
			$tran_go        = $species_info{$species}{tran_go};
			$tran_kegg      = $species_info{$species}{tran_kegg};
			
			$prot_seq       = $species_info{$species}{prot_seq};
			$prot_trans     = $species_info{$species}{prot_trans};
		}
	}

	if(!defined( $refer_seq )){
		print STDERR "Error: -refer_seq must be specified!\n";
		&usage;
	}elsif(!(-e $refer_seq)){
		print STDERR "Error: [$refer_seq] doses not exist!\n";
		&usage;
	}else{
		$refer_seq  = abs_path( $refer_seq  );
	}
	
	if(!defined( $refer_gtf )){
		print STDERR "Error: -refer_gtf must be specified!\n";
		&usage;
	}elsif(!(-e $refer_gtf)){
		print STDERR "Error: [$refer_gtf] doses not exist!\n";
		&usage;
	}else{
		$refer_gtf  = abs_path( $refer_gtf  );
	}
	
	# if(!defined( $gene_go )){
		# print STDERR "Warning: -gene_go is not specified!\n";
		# $annot_known = 1;
	# }elsif(!(-e $gene_go)){
		# print STDERR "Error: -gene_go [$gene_go] doses not exist!\n";
		# &usage;
	# }else{
		# $gene_go  = abs_path( $gene_go );
	# }
	
	# if(!defined( $gene_kegg )){
		# print STDERR "Warning: -gene_kegg is not specified!\n";
		# $annot_known = 1;
	# }elsif(!(-e $gene_kegg)){
		# print STDERR "Error: -gene_kegg [$gene_kegg] doses not exist!\n";
		# &usage;
	# }else{
		# $gene_kegg  = abs_path( $gene_kegg );
	# }
	
	# if(!defined( $gene_annot )){
		# print STDERR "Warning: -gene_annot is not specified!\n";
		# $annot_known = 1;
	# }elsif(!(-e $gene_annot)){
		# print STDERR "Error: -gene_annot [$gene_annot] doses not exist!\n";
		# &usage;
	# }else{
		# $gene_annot  = abs_path( $gene_annot );
	# }
	
	# if(!defined( $prot_seq )){
		# print STDERR "Warning: -prot_seq is not specified!\n";
		# $annot_known = 1;
	# }elsif(!(-e $prot_seq)){
		# print STDERR "Error: -prot_seq [$prot_seq] doses not exist!\n";
		# &usage;
	# }else{
		# $prot_seq  = abs_path( $prot_seq );
		# if(!defined( $prot_trans )){
			# print STDERR "Warning: -prot_trans is not specified!\n";
			# $annot_known = 1;
		# }elsif(!(-e $prot_trans)){
			# print STDERR "Error: -prot_trans [$prot_trans] doses not exist!\n";
			# &usage;
		# }else{
			# $prot_trans  = abs_path( $prot_trans );
		# }
	# }
	
	# if($run_trans){
		# if(!defined( $tran_go )){
			# print STDERR "Warning: -tran_go is not specified!\n";
			# $annot_known = 1;
		# }elsif(!(-e $tran_go)){
			# print STDERR "Error: -tran_go [$tran_go] doses not exist!\n";
			# &usage;
		# }else{
			# $tran_go  = abs_path( $tran_go );
		# }
		
		# if(!defined( $tran_kegg )){
			# print STDERR "Warning: -tran_kegg is not specified!\n";
			# $annot_known = 1;
		# }elsif(!(-e $tran_kegg)){
			# print STDERR "Error: -tran_kegg [$tran_kegg] doses not exist!\n";
		# &usage;
		# }else{
			# $tran_kegg  = abs_path( $tran_kegg );
		# }
		
		# if(!defined( $tran_annot )){
			# print STDERR "Warning: -tran_annot is not specified!\n";
			# $annot_known = 1;
		# }elsif(!(-e $tran_annot)){
			# print STDERR "Error: -tran_annot [$tran_annot] doses not exist!\n";
		# &usage;
		# }else{
			# $tran_annot  = abs_path( $tran_annot );
		# }
	# }
	
	# unless($ppi_taxonomy){
		# print STDERR "Error: -ppi_taxonomy must be specified!\n";
		# &usage;
	# }
	
	$queue_level ||= "run_strict";
	unless(grep /^$queue_level$/, qw/only_cmd run_strict/){
		print STDERR "Error: -queue_level must be one of (only_cmd run_strict)!\n";
		&usage;
	}
	
	$strand  ||= "no";
	if($strand ne "no" && $strand ne "yes"){
		print STDERR "Error: -strand must be yes or no!\n";
		&usage;
	}
	
	$method_exp ||= "A1";
	unless(grep /^$method_exp$/, qw/A1 A2 A3 B1 B2 B3 C1 C2/){
		print STDERR "Error: -method_exp must be one of (A1 A2 A3 B1 B2 B3 C1 C2)!\n";
		&usage;
	}
	
	if(grep(/^$method_exp$/, qw/A1 A2 A3/) && !defined( $index_hisat2 )){
		print STDERR "Warning: -index_hisat2 is not specified, the pipelien will construct it automatically!\n";
		$run_index_hisat2 = 1;
	}
	
	if(grep(/^$method_exp$/, qw/B1 B2 B3/) && !defined( $index_star )){
		print STDERR "Warning: -index_star is not specified, the pipelien will construct it automatically!\n";
		$run_index_star   = 1;
	}
	
	if(grep(/^$method_exp$/, qw/C1/) && !defined( $index_salmon )){
		print STDERR "Warning: -index_star is not specified, the pipelien will construct it automatically!\n";
		$run_index_salmon   = 1;
	}
	
	if(grep(/^$method_exp$/, qw/C2/) && !defined( $index_kallisto )){
		print STDERR "Warning: -index_kallisto is not specified, the pipelien will construct it automatically!\n";
		$run_index_kallisto   = 1;
	}
	
	$de_method ||= "edgeR";
	unless(grep /^$de_method$/, qw/edgeR DESeq2 DEGseq limma NOIseq ttest wilcox/){
		print STDERR "Error: -method_exp must be one of (edgeR DESeq2 DEGseq limma NOIseq ttest wilcox)!\n";
		&usage;
	}
	if($run_novel){
		#@exp_class = ("total","known","novel");
		@exp_class = ("total", "known");
	}else{
		@exp_class = ("known");
	}
	
	if($run_trans){
		@exp_level = ("gene", "transcript");
	}else{
		@exp_level = ("gene");
	}
	
	if($run_snp && !defined( $index_star ) ){
		$run_index_star = 1;
	}
	
	$de_fc     = 2    unless(defined $de_fc);
	$de_pval   = 0    unless(defined $de_pval);
	$de_padj   = 0.05 unless(defined $de_padj);
	if($de_fc > 0){
		$de_log2fc = log($de_fc)/log(2);
	}else{
		$de_log2fc = 0;
	}
	
	$NR        = "nr";
	$KEGG      = "kegg";
	$SWISSPROT = "swissprot";
	$STRING    = "string";

	if(defined($class)){
		my $flag=0;
		foreach my $x(qw/animal mammal vertebrate invertebrate plant fungi/){
			if($x eq $class){
				$flag++;
			}
		}
		unless($flag){
			die "class must be one of animal,mammal,vertebrate,invertebrate,plant,fungi\n";
		}else{
			$NR   = $NR."_".$class;
			$KEGG = $KEGG."_".$class;
		}
	}
	
	$kegg_species||= "ko";
	$partition   ||= "small";
	$job_prefix  ||= "mrna";
	$queue_rerun ||= 0;
	$job_mode    ||= "pbs";
	$max_cpus    ||= 500;
	$max_mems    ||= 1000;
	$max_jobs    ||= 50;
	$outdir      ||= "out";

	$outdir = abs_path( $outdir );
	system("mkdir -p $outdir");
	
	#&parsing_background();
}

sub load_species_info{
	my ($configfile, $hash) = @_;
	
	unless($configfile){
		print STDERR "species cfg must be specified!\n";
		&usage;
	}
	unless(-e $configfile){
		print STDERR "species cfg file $configfile does not exist!";
		&usage;
	}
	my $cfg = Config::IniFiles->new( -file => $configfile );
	my @sec = $cfg->Sections();
	foreach my $s (@sec) {
		my @pa = $cfg->Parameters($s);
		foreach my $p (@pa) {
			$$hash{$s}{$p} = $cfg->val( $s, $p );
		}
	}
}

sub parsing_sample_info{
	open  FIN,$file_input or die "can not open raw fq list file!";
	while(<FIN>){
		chomp;
		next if /^#/;
		my ( $sample, $fq1, $fq2 ) = split /\t/,$_;
		if( $sample !~ /^[a-zA-Z]/ ){
			warn "sample_name [$sample] must start with 大小写英文字母!";
		}
		if( $sample =~ /[^a-zA-Z0-9_\-\.]/ ){
			die "sample_name [$sample] must consisted by 大小写英文字母、数字、横杠、小数点以及下划线!";
		}
		if(exists $sample_info{$sample}){
			die "sample_name [$sample] is duplicate! please check input file[$file_input]!";
		}
		if($fq1 eq ""){
			die "Error: sample [$sample] must have fq data!\n";
		}
		if($fq2 ne ""){
			print STDERR "[INFO]: sample [$sample] is pair-end!\n";
			$sample_info{$sample}{R1}   = $fq1;
			$sample_info{$sample}{R2}   = $fq2;
			$sample_info{$sample}{type} = "PE";
		}else{
			print STDERR "[INFO]: sample [$sample] is single-end!\n";
			$sample_info{$sample}{S}    = $fq1;
			$sample_info{$sample}{type} = "SE";
		}
	}
	close FIN;
}

sub parsing_meta_and_comp_info{
	system("mkdir -p $outdir/group_info");
	my %www;
	my %all_std;
	if(defined $file_meta && defined $file_comp){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		my $head = <FIN>; chomp($head);
		unless($head =~ "^Sample\t"){
			die "Meta file format error[header], please check!";
		}
		my ($a,@grp_stds) = split /\t/,$head;
		foreach my $x(@grp_stds){
			if( $x !~ /^[a-zA-Z]/ ){
				warn "group standard [$x] must start with 大小写英文字母!";
			}
			if( $x =~ /[^a-zA-Z0-9_\-\.]/ ){
				die "group standard [$x] must consisted by 大小写英文字母、数字、横杠、小数点以及下划线!";
			}
		}
		while(<FIN>){
			chomp;
			my ($sn,@grps) = split /\t/,$_;
			foreach my $x(@grps){
				if( $x !~ /^[a-zA-Z]/ ){
					warn "group name [$x] must start with 大小写英文字母!";
				}
				if( $x =~ /[^a-zA-Z0-9_\-\.]/ ){
					warn "group name [$x] must consisted by 大小写英文字母、数字、横杠、小数点以及下划线!";
				}
			}
			unless(exists $sample_info{$sn}){
				die "sample_name[$sn] in meta file does not existe in fasta file!";
			}
			if(scalar @grps != scalar @grp_stds){
				die "meta file format error[content], please check!";
			}
			for(my $i = 0; $i < scalar(@grp_stds); $i++){
				if($grps[$i] ne "NA" && $grps[$i] ne ""){
					$meta_grp2sn{$grp_stds[$i]}{$grps[$i]}{$sn}++;
					$meta_sn2grp{$sn}{$grp_stds[$i]}=$grps[$i];
					$www{$grp_stds[$i]."__".$grps[$i]}{$sn}++;
				}
			}
			$meta_grp2sn{"Sample"}{$sn}{$sn}++;
			$meta_sn2grp{$sn}{"Sample"}= $sn;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			next if /^#/;
			my $name;
			my ($grp_std,$comps,$method,$cutoff_pval,$cutoff_padj,$cutoff_fc) = split /\t/,$_;
			
			unless(exists $meta_grp2sn{$grp_std}){
				die "group std name [$grp_std] does not exists in meta info file [$file_meta]!";
			}
			my @grps =split /,/,$comps;
			if( scalar(@grps) > 1 ){
				$name = $grp_std.".".join("_vs_",@grps);
			}else{
				die "the number of elements in comparison must larger than 1!";
			}
			
			if(exists $all_std{$name}){
				die "std[$name] was already exists!!!";
			}else{
				$all_std{$name}++;
			}

			foreach my $grp(@grps){
				die "group name [$grp] does not exists in group std [$grp_std] !!" unless exists $meta_grp2sn{$grp_std}{$grp};
				foreach my $sn(sort keys %{$meta_grp2sn{$grp_std}{$grp}}){
					$comp_grp2sn{$name}{$grp}{$sn}++;
					$comp_sn2grp{$sn}{$name} = $grp;
				}
			}
			
			open  FOUT,">$outdir/group_info/$name.cts";
			for(my $i=0; $i < (scalar(@grps)-1); $i++){
				for(my $j=$i+1; $j< scalar(@grps); $j++){
					print FOUT $grps[$i],"\t",$grps[$j],"\n";
				}
			}
			close FOUT;
			
			open  F1,">$outdir/group_info/$name.info";
			open  F2,">$outdir/group_info/$name.txt";
			print F1 "sample\tgroup\n";
			foreach my $grp(sort keys %{$comp_grp2sn{$name}}){
				foreach my $sn(sort keys %{$comp_grp2sn{$name}{$grp}}){
					print F1 $sn,"\t",$grp,"\n";
					print F2 $sn,"\t",$grp,"\n";
				}
			}
			close F1;
			close F2;
			
			my $flag_norep = 0;
			foreach my $grp(@grps){
				if( scalar keys %{$comp_grp2sn{$name}{$grp}} <= 1 ){
					$flag_norep++;
				}
			}
			
			unless($method){
				if($flag_norep){
					$method  = "edgeR";
				}else{
					$method  = "DESeq2";
				}
				$cutoff_pval = $de_pval;
				$cutoff_padj = $de_padj;
				$cutoff_fc   = $de_fc;
			}
			$comp_param{$name}{method} = $method;
			$comp_param{$name}{pval}   = $cutoff_pval;
			$comp_param{$name}{padj}   = $cutoff_padj;
			$comp_param{$name}{fc}     = $cutoff_fc;
		}
		close FIN;
		
	}elsif(defined $file_comp && !(defined $file_meta)){
		open FIN,"$file_input";
		while(<FIN>){
			chomp;
			my ($sn) = split /\t/,$_;
			$meta_grp2sn{"Sample"}{$sn}{$sn}++;
			$meta_sn2grp{$sn}{"Sample"}= $sn;
			#$www{"Sample"."__".$sn}{$sn}++;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			next if /^#/;
			my $name;
			my ($grp_std,$comps,$method,$cutoff_pval,$cutoff_padj,$cutoff_fc) = split /\t/,$_;
			
			unless(exists $meta_grp2sn{$grp_std}){
				die "group std name [$grp_std] does not exists in meta info file [$file_meta]!";
			}
			my @grps =split /,/,$comps;
			if( scalar(@grps) > 1 ){
				$name = $grp_std.".".join("_vs_",@grps);
			}else{
				die "the number of elements in comparison must larger than 1!";
			}
			
			if(exists $all_std{$name}){
				die "std[$name] was already exists!!!";
			}else{
				$all_std{$name}++;
			}

			foreach my $grp(@grps){
				die "group name [$grp] does not exists in group std [$grp_std] !!" unless exists $meta_grp2sn{$grp_std}{$grp};
				foreach my $sn(sort keys %{$meta_grp2sn{$grp_std}{$grp}}){
					$comp_grp2sn{$name}{$grp}{$sn}++;
					$comp_sn2grp{$sn}{$name} = $grp;
				}
			}
			
			open  FOUT,">$outdir/group_info/$name.cts";
			for(my $i=0; $i < (scalar(@grps)-1); $i++){
				for(my $j=$i+1; $j< scalar(@grps); $j++){
					print FOUT $grps[$i],"\t",$grps[$j],"\n";
				}
			}
			close FOUT;
			
			open  F1,">$outdir/group_info/$name.info";
			open  F2,">$outdir/group_info/$name.txt";
			print F1 "sample\tgroup\n";
			foreach my $grp(sort keys %{$comp_grp2sn{$name}}){
				foreach my $sn(sort keys %{$comp_grp2sn{$name}{$grp}}){
					print F1 $sn,"\t",$grp,"\n";
					print F2 $sn,"\t",$grp,"\n";
				}
			}
			close F1;
			close F2;
			
			my $flag_norep = 0;
			foreach my $grp(@grps){
				if( scalar keys %{$comp_grp2sn{$name}{$grp}} <= 1 ){
					$flag_norep++;
				}
			}
			
			unless($method){
				if($flag_norep){
					$method  = "edgeR";
				}else{
					$method  = "DESeq2";
				}
				$cutoff_pval = $de_pval;
				$cutoff_padj = $de_padj;
				$cutoff_fc   = $de_fc;
			}
			$comp_param{$name}{method} = $method;
			$comp_param{$name}{pval}   = $cutoff_pval;
			$comp_param{$name}{padj}   = $cutoff_padj;
			$comp_param{$name}{fc}     = $cutoff_fc;
		}
		close FIN;
		
	}elsif(defined $file_meta && !(defined $file_comp)){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		my $head = <FIN>; chomp($head);
		unless($head =~ "^Sample\t"){
			die "Meta file format error[header], please check!";
		}
		my ($a,@grp_stds) = split /\t/,$head;
		foreach my $x(@grp_stds){
			if( $x !~ /^[a-zA-Z]/ ){
				warn "group standard [$x] must start with 大小写英文字母!";
			}
			if( $x =~ /[^a-zA-Z0-9_\-\.]/ ){
				die "group standard [$x] must consisted by 大小写英文字母、数字、横杠、小数点以及下划线!";
			}
		}
		while(<FIN>){
			chomp;
			my ($sn,@grps) = split /\t/,$_;
			foreach my $x(@grps){
				if( $x !~ /^[a-zA-Z]/ ){
					warn "group name [$x] must start with 大小写英文字母!";
				}
				if( $x =~ /[^a-zA-Z0-9_\-\.]/ ){
					warn "group name [$x] must consisted by 大小写英文字母、数字、横杠、小数点以及下划线!";
				}
			}
			unless(exists $sample_info{$sn}){
				die "sample_name[$sn] in meta file does not existe in fasta file!";
			}
			if(scalar @grps != scalar @grp_stds){
				die "meta file format error[content], please check!";
			}
			for(my $i = 0; $i < scalar(@grp_stds); $i++){
				if($grps[$i] ne "NA" && $grps[$i] ne ""){
					$meta_grp2sn{$grp_stds[$i]}{$grps[$i]}{$sn}++;
					$meta_sn2grp{$sn}{$grp_stds[$i]}=$grps[$i];
					$www{$grp_stds[$i]."__".$grps[$i]}{$sn}++;
				}
			}
			$meta_grp2sn{"Sample"}{$sn}{$sn}++;
			$meta_sn2grp{$sn}{"Sample"}= $sn;
		}
		close FIN;
	}else{
		open FIN,"$file_input";
		while(<FIN>){
			chomp;
			my ($sn) = split /\t/,$_;
			$meta_grp2sn{"Sample"}{$sn}{$sn}++;
			$meta_sn2grp{$sn}{"Sample"}= $sn;
			#$www{"Sample"."__".$sn}{$sn}++;
		}
		close FIN;
	}
	
	if(scalar (keys %www) > 0){
		open FOUT,">$outdir/group_info/group_info.num";
		print FOUT "id\t".join("\t",sort keys %www)."\n";
		foreach my $sn(sort keys %meta_sn2grp){
			print FOUT $sn;
			foreach my $name(sort keys %www){
				if(exists $www{$name}{$sn}){
					print FOUT "\t",1;
				}else{
					print FOUT "\t",0;
				}
			}
			print FOUT "\n";
		}
		close FOUT;
	}
}

sub rawdata_quality_control{
	foreach my $sn(sort keys %sample_info){
		$all_dependencies{"rawdata_quality_control.sample.$sn"}++;
		$queue->set_job_cpu(20);
		$queue->set_job_mem(40);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("rawdata_quality_control.sample.$sn");
		$queue->set_job_disc("样本[$sn]原始数据质量控制分析");
		$queue->set_job_depend("");
		if($sample_info{$sn}{type} eq "SE"){
			$queue->set_file_check("$outdir/quality_control/clean/$sn.clean.fq.gz,$outdir/quality_control/fastp/$sn.fastp.html,$outdir/quality_control/stats/$sn.summary_qc.xls");
		}else{
			$queue->set_file_check("$outdir/quality_control/clean/$sn.clean.R1.fq.gz,$outdir/quality_control/clean/$sn.clean.R2.fq.gz,$outdir/quality_control/fastp/$sn.fastp.html,$outdir/quality_control/stats/$sn.summary_qc.xls");
		}
		
		$queue->set_work_dir($workdir);
		$queue->addcommonds("module load fastp/0.21.0");
		$queue->addcommonds("module load seqtk/1.3");
		$queue->addcommonds("module load ucsc/1.0");
		$queue->addcommonds("module load fastqc/0.11.8");
		$queue->addcommonds("module load blast/2.2.26");
		
		$queue->addcommonds("mkdir -p $outdir/quality_control/raw");
		$queue->addcommonds("mkdir -p $outdir/quality_control/clean");
		$queue->addcommonds("mkdir -p $outdir/quality_control/fastqc");
		$queue->addcommonds("mkdir -p $outdir/quality_control/fastp");
		$queue->addcommonds("mkdir -p $outdir/quality_control/rfam");
		$queue->addcommonds("mkdir -p $outdir/quality_control/stats");
		
		if($sample_info{$sn}{type} eq "SE"){
			my @fqs = split /,/,$sample_info{$sn}{S};
			
			if(scalar(@fqs) == 1){
				foreach my $fq(@fqs){
					if($fq =~ /.gz$/){
						$queue->addcommonds("cp $fq $outdir/quality_control/raw/$sn.raw.fq.gz");
					}else{
						$queue->addcommonds("gzip $fq -c > $outdir/quality_control/raw/$sn.raw.fq.gz");
					}
				}
			}else{
				$queue->addcommonds("rm -f $outdir/quality_control/raw/$sn.raw.fq");
				foreach my $fq(@fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/quality_control/raw/$sn.raw.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/quality_control/raw/$sn.raw.fq");
					}
				}
				$queue->addcommonds("gzip $outdir/quality_control/raw/$sn.raw.fq");
			}
			
			$queue->addcommonds("fastp --in1 $outdir/quality_control/raw/$sn.raw.fq.gz --out1 $outdir/quality_control/clean/$sn.clean.fq.gz --unpaired1 $outdir/quality_control/clean/$sn.unpair.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 20 --thread 16 --json $outdir/quality_control/fastp/$sn.fastp.json --html $outdir/quality_control/fastp/$sn.fastp.html --report_title '$sn fastp report' &> $outdir/quality_control/fastp/$sn.fastp.log");
			$queue->addcommonds("$script/parse_fastp.pl --json $outdir/quality_control/fastp/$sn.fastp.json --outdir $outdir/quality_control/stats --prefix $sn");
			
			$queue->addcommonds("seqtk sample $outdir/quality_control/clean/$sn.clean.fq.gz 10000 > $outdir/quality_control/rfam/$sn.sample.fq");
			$queue->addcommonds("fastqToFa $outdir/quality_control/rfam/$sn.sample.fq $outdir/quality_control/rfam/$sn.sample.fa");
			$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/quality_control/rfam/$sn.sample.fa -o $outdir/quality_control/rfam/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
			$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/quality_control/rfam/$sn.rfam_blast.txt > $outdir/quality_control/rfam/$sn.rfam_blast.xls");
			$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/quality_control/rfam/$sn.rfam_blast.txt -fasta $outdir/quality_control/rfam/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/quality_control/rfam/ --prefix $sn -maxEvalue 0.01");
			
			foreach my $x(qw/raw clean/){
				$queue->addcommonds("mkdir -p $outdir/quality_control/fastqc/$x/$sn.fastqc");
				$queue->addcommonds("fastqc -t 20 --nogroup --extract $outdir/quality_control/$x/$sn.$x.fq.gz --outdir $outdir/quality_control/fastqc/$x/$sn.fastqc   &> $outdir/quality_control/fastqc/$x/$sn.fastqc.log");
			}
		}else{
			my @R1_fqs = split /,/,$sample_info{$sn}{R1};
			my @R2_fqs = split /,/,$sample_info{$sn}{R2};
			
			if(scalar(@R1_fqs) != scalar(@R2_fqs)){
				die "[Error]: sample[$sn] has diff number of files in R1 and R2!\n";
			}elsif(scalar(@R1_fqs) == 1){
				foreach my $fq(@R1_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("cp $fq $outdir/quality_control/raw/$sn.raw.R1.fq.gz");
					}else{
						$queue->addcommonds("gzip $fq -c > $outdir/quality_control/raw/$sn.raw.R1.fq.gz");
					}
				}
				foreach my $fq(@R2_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("cp $fq $outdir/quality_control/raw/$sn.raw.R2.fq.gz");
					}else{
						$queue->addcommonds("gzip $fq -c > $outdir/quality_control/raw/$sn.raw.R2.fq.gz");
					}
				}
			}else{
				$queue->addcommonds("rm -f $outdir/quality_control/raw/$sn.raw.R1.fq");
				$queue->addcommonds("rm -f $outdir/quality_control/raw/$sn.raw.R2.fq");
				foreach my $fq(@R1_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/quality_control/raw/$sn.raw.R1.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/quality_control/raw/$sn.raw.R1.fq");
					}
				}
				foreach my $fq(@R2_fqs){
					$fq =~ s/ //g;
					if($fq =~ /.gz$/){
						$queue->addcommonds("zcat $fq >> $outdir/quality_control/raw/$sn.raw.R2.fq");
					}else{
						$queue->addcommonds("cat  $fq >> $outdir/quality_control/raw/$sn.raw.R2.fq");
					}
				}
				$queue->addcommonds("gzip $outdir/quality_control/raw/$sn.raw.R1.fq");
				$queue->addcommonds("gzip $outdir/quality_control/raw/$sn.raw.R2.fq");
			}
			$queue->addcommonds("fastp --in1 $outdir/quality_control/raw/$sn.raw.R1.fq.gz --in2 $outdir/quality_control/raw/$sn.raw.R2.fq.gz --out1 $outdir/quality_control/clean/$sn.clean.R1.fq.gz --out2 $outdir/quality_control/clean/$sn.clean.R2.fq.gz --unpaired1 $outdir/quality_control/clean/$sn.unpair.R1.fq.gz --unpaired2 $outdir/quality_control/clean/$sn.unpair.R2.fq.gz --compression 4 --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --correction --length_required 20 --thread 16 --json $outdir/quality_control/fastp/$sn.fastp.json --html $outdir/quality_control/fastp/$sn.fastp.html --report_title '$sn fastp report' &> $outdir/quality_control/fastp/$sn.fastp.log");
			$queue->addcommonds("$script/parse_fastp.pl --json $outdir/quality_control/fastp/$sn.fastp.json --outdir $outdir/quality_control/stats --prefix $sn");
			
			$queue->addcommonds("seqtk sample $outdir/quality_control/clean/$sn.clean.R1.fq.gz 10000 > $outdir/quality_control/rfam/$sn.sample.fq");
			$queue->addcommonds("fastqToFa $outdir/quality_control/rfam/$sn.sample.fq $outdir/quality_control/rfam/$sn.sample.fa");
			$queue->addcommonds("blastall -p blastn -d Rfam -i $outdir/quality_control/rfam/$sn.sample.fa -o $outdir/quality_control/rfam/$sn.rfam_blast.txt -v 5 -b 5 -F F -e 0.01 -a 10");
			$queue->addcommonds("$script/Blast2table -format 10 -expect 1E-2 -top $outdir/quality_control/rfam/$sn.rfam_blast.txt > $outdir/quality_control/rfam/$sn.rfam_blast.xls");
			$queue->addcommonds("$script/rfam_blast_parse.pl -blast $outdir/quality_control/rfam/$sn.rfam_blast.txt -fasta $outdir/quality_control/rfam/$sn.sample.fa -seed /mnt/nfs/users/pub/database/blastdb/blast/Rfam.seed -outdir $outdir/quality_control/rfam/ --prefix $sn -maxEvalue 0.01");
			
			foreach my $x(qw/raw clean/){
				foreach my $y(qw/R1 R2/){
					$queue->addcommonds("mkdir -p $outdir/quality_control/fastqc/$x/$sn.$y.fastqc");
					$queue->addcommonds("fastqc -t 20 --nogroup --extract $outdir/quality_control/$x/$sn.$x.$y.fq.gz --outdir $outdir/quality_control/fastqc/$x/$sn.$y.fastqc &> $outdir/quality_control/fastqc/$x/$sn.$y.fastqc.log");
				}
			}
		}
		$queue->run();
		$results{rawdata_quality_control}++;
	}
	
	$all_dependencies{"rawdata_quality_control.merge"}++;
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_mode($job_mode);
	$queue->set_job_name("rawdata_quality_control.merge");
	$queue->set_job_disc("所有样本原始数据质量控制分析整合");
	my @p = ();
	foreach my $sn(sort keys %sample_info){
		push @p, "rawdata_quality_control.sample.$sn";
	}
	$queue->set_job_depend(join(",",@p));
	$queue->set_file_check("$outdir/quality_control/stats/cleandata.stat.xls,$outdir/quality_control/stats/rawdata.stat.xls,$outdir/quality_control/stats/summary.filtering.xls,$outdir/quality_control/stats/summary.qc.xls");
	$queue->set_work_dir($workdir);

	my $cmd;
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/quality_control/stats/rawdata.stat.xls";
	foreach my $sn(sort keys %sample_info){
		$cmd .= " $outdir/quality_control/stats/$sn.rawdata.stat.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/quality_control/stats/cleandata.stat.xls";
	foreach my $sn(sort keys %sample_info){
		$cmd .= " $outdir/quality_control/stats/$sn.cleandata.stat.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/quality_control/stats/summary.filtering.xls";
	foreach my $sn(sort keys %sample_info){
		$cmd .= " $outdir/quality_control/stats/$sn.summary_filter.xls";
	}
	$queue->addcommonds("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/quality_control/stats/summary.qc.xls";
	foreach my $sn(sort keys %sample_info){
		$cmd .= " $outdir/quality_control/stats/$sn.summary_qc.xls";
	}
	$queue->addcommonds("$cmd");
	
	$queue->addcommonds("grep '^rRNA' $outdir/quality_control/rfam/*.rfam_summary.xls | sed 's#$outdir/quality_control/rfam/##g' | cut -f 2 -d '/' | sed 's/.rfam_summary.xls:rRNA//g' | sed 's/\\%//g'| cut -f 1,3 | sed '1 iSample_ID\trRNA_Rate(%)' > $outdir/quality_control/stats/rRNA_percentage.xls");
	
	$queue->addcommonds("module load app/table_tools");
	$queue->addcommonds("col_add.pl -t $outdir/quality_control/stats/summary.qc.xls -i $outdir/quality_control/stats/rRNA_percentage.xls -n 1 -headi T -headt T -fill T > $outdir/quality_control/stats/summary.quality_control.xls");
	
	$queue->addcommonds("source /mnt/nfs/users/pub/software/python/3.7.5/env/MultiQC-1.9/bin/activate");
	$queue->addcommonds("multiqc $outdir/quality_control/fastp --outdir $outdir/quality_control/fastp/multiqc --force --export");
	$queue->addcommonds("multiqc $outdir/quality_control/fastqc/raw --outdir $outdir/quality_control/fastqc/raw/multiqc --force --export");
	$queue->addcommonds("multiqc $outdir/quality_control/fastqc/clean --outdir $outdir/quality_control/fastqc/clean/multiqc --force --export");
	$queue->addcommonds("deactivate");
	
	$queue->run();
	$results{rawdata_quality_control}++;
}

sub reference_preparation_known_1 {
	$all_dependencies{"reference_preparation_known.part_1"}++;
	$queue->set_job_cpu(8);
	$queue->set_job_mem(10);
	$queue->set_job_mode($job_mode);
	$queue->set_job_name("reference_preparation_known.part_1");
	$queue->set_job_disc("基因注释文件整理[known].part_1");
	$queue->set_file_check("$outdir/refer/annot_known.refFlat,$outdir/refer/annot_known.gtf,$outdir/refer/annot_known.transcript.fa,$outdir/refer/annot_known.g2l");
	$queue->set_work_dir($workdir);
	
	$queue->addcommonds("module load ucsc/1.0");
	$queue->addcommonds("module load app/sequence_tools");
	$queue->addcommonds("module load app/gtf_tools");
	$queue->addcommonds("module load app/table_tools");
	$queue->addcommonds("module load bedtools/2.25.0");
	$queue->addcommonds("module load GTFtools/0.6.0");
	
	$queue->addcommonds("mkdir -p $outdir/refer");
	$queue->addcommonds("ln -sf $refer_seq $outdir/refer/genome.fa");
	$queue->addcommonds("$script/count.refer.pl $outdir/refer/genome.fa");
	$queue->addcommonds("ln -sf $refer_gtf $outdir/refer/annot_known.gtf");
	$queue->addcommonds("gtf_g2n.pl -gtf $outdir/refer/annot_known.gtf   -out $outdir/refer/annot_known.gene2name");
	$queue->addcommonds("gtf_t2g2n.pl -gtf $outdir/refer/annot_known.gtf -out $outdir/refer/annot_known.transcript2gene2name");
	$queue->addcommonds("cut -f 1,2 $outdir/refer/annot_known.transcript2gene2name > $outdir/refer/annot_known.transcript2gene");
	$queue->addcommonds("awk '{print \$2\"\\t\"\$1}' $outdir/refer/annot_known.transcript2gene2name > $outdir/refer/annot_known.gene2transcript");
	$queue->addcommonds("cut -f 1 $outdir/refer/annot_known.gene2name | sed 1d | sort -u  >  $outdir/refer/annot_known.gene.list");
	$queue->addcommonds("cut -f 1 $outdir/refer/annot_known.transcript2gene2name | sed 1d | sort -u  >  $outdir/refer/annot_known.transcript.list");
	$queue->addcommonds("gtfToGenePred $outdir/refer/annot_known.gtf $outdir/refer/annot_known.genePred");
	$queue->addcommonds("awk -v OFS='\\t' '{print \$1,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' $outdir/refer/annot_known.genePred >  $outdir/refer/annot_known.refFlat");
	
	$queue->addcommonds("gtf2bed.pl -input $outdir/refer/annot_known.gtf -out $outdir/refer/annot_known.exon.bed");
	$queue->addcommonds("shuf -n5000 $outdir/refer/annot_known.exon.bed | sortBed -i - > $outdir/refer/annot_known.exon_random.bed");
	$queue->addcommonds("bedtools getfasta -fi $outdir/refer/genome.fa -bed $outdir/refer/annot_known.exon.bed -fo $outdir/refer/annot_known.transcript.fa -name -split -s");
	
	$queue->addcommonds("gtftools-modified.py -l $outdir/refer/annot_known.g2l $outdir/refer/annot_known.gtf");
	$queue->addcommonds("cut -f 1,5 $outdir/refer/annot_known.g2l | sed 1d > $outdir/refer/annot_known.gene2length");
	$queue->addcommonds("fastalength $outdir/refer/annot_known.transcript.fa | awk '{print \$2\"\\t\"\$1}' > $outdir/refer/annot_known.transcript2length");
	
	$queue->run();
	$results{reference_preparation_known}++;
}

sub genome_index {
	if($run_index_hisat2){
		$all_dependencies{"index.hisat2"}++;
		$queue->set_job_cpu(20);
		$queue->set_job_mem(20);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("index.hisat2");
		$queue->set_job_disc("构建index[hisat2]");
		$queue->set_job_depend("reference_preparation_known.part_1");
		$queue->set_file_check("$outdir/refer/index_hisat2.1.ht2,$outdir/refer/index_hisat2.2.ht2,$outdir/refer/index_hisat2.6.ht2");
		$queue->set_work_dir($workdir);
		$queue->addcommonds("module load hisat/2.1.0");
		$queue->addcommonds("hisat2-build $outdir/refer/genome.fa $outdir/refer/index_hisat2 -p 20");
		$queue->run();
		$index_hisat2 = "$outdir/refer/index_hisat2";
	}
	if($run_index_star){
		$all_dependencies{"index.star"}++;
		$queue->set_job_cpu(20);
		$queue->set_job_mem(20);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("index.star");
		$queue->set_job_disc("构建index[star]");
		$queue->set_job_depend("reference_preparation_known.part_1");
		$queue->set_file_check("$outdir/refer/index_star/chrLength.txt");
		$queue->set_work_dir($workdir);
		$queue->addcommonds("module load star/2.7.2b");
		$queue->addcommonds("mkdir -p $outdir/refer/index_star");
		$queue->addcommonds("cd $outdir/refer/index_star");
		$queue->addcommonds("ln -s $outdir/refer/genome.fa genome.fa");
		$queue->addcommonds("ln -s $outdir/refer/annot_known.gtf gene.gtf");
		$queue->addcommonds("STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles genome.fa --sjdbGTFfile gene.gtf --sjdbOverhang 150");
		$queue->addcommonds("cd $workdir");
		$queue->run();
		$index_star = "$outdir/refer/index_star";
	}
}

sub genome_alignment {
	if(grep /^$method_exp$/, qw/A1 A2 A3/){
		foreach my $sn(sort keys %sample_info){
			$all_dependencies{"alignment.sample.$sn"}++;
			$queue->set_job_cpu(20);
			$queue->set_job_mem(40);
			$queue->set_job_mode($job_mode);
			$queue->set_job_name("alignment.sample.$sn");
			$queue->set_job_disc("样本[$sn]比对基因组(hisat2)");
			if($run_index_hisat2){
				$queue->set_job_depend("reference_preparation_known.part_1,index.hisat2,rawdata_quality_control.sample.$sn");
			}else{
				$queue->set_job_depend("reference_preparation_known.part_1,rawdata_quality_control.sample.$sn");
			}
			
			$queue->set_file_check("$outdir/alignment/$sn.bam,$outdir/alignment/$sn.bam.bai");
			$queue->set_work_dir($workdir);
			
			$queue->addcommonds("module load hisat/2.1.0");
			$queue->addcommonds("module load samtools/1.6");
			$queue->addcommonds("mkdir -p $outdir/alignment");
			if($strand eq "yes"){
				$queue->addcommonds("hisat2 -q -x $index_hisat2 -1 $outdir/quality_control/clean/$sn.clean.R1.fq.gz -2 $outdir/quality_control/clean/$sn.clean.R2.fq.gz -p 20 --novel-splicesite-outfile $outdir/alignment/$sn.splicesite.novel --novel-splicesite-infile $outdir/alignment/$sn.splicesite.novel -S $outdir/alignment/$sn.sam --dta --summary-file $outdir/alignment/$sn.hisat2_summary.txt --rna-strandness RF");
			}else{
				$queue->addcommonds("hisat2 -q -x $index_hisat2 -1 $outdir/quality_control/clean/$sn.clean.R1.fq.gz -2 $outdir/quality_control/clean/$sn.clean.R2.fq.gz -p 20 --novel-splicesite-outfile $outdir/alignment/$sn.splicesite.novel --novel-splicesite-infile $outdir/alignment/$sn.splicesite.novel -S $outdir/alignment/$sn.sam --dta --summary-file $outdir/alignment/$sn.hisat2_summary.txt");
			}
			$queue->addcommonds("samtools view -S -b $outdir/alignment/$sn.sam | samtools sort - -\@ 20 -o $outdir/alignment/$sn.bam");
			$queue->addcommonds("samtools index $outdir/alignment/$sn.bam -\@ 20");
			$queue->addcommonds("rm -rf $outdir/alignment/$sn.sam");
			
			$queue->addcommonds("source /mnt/nfs/users/pub/software/python/3.7.5/env/RSeQC-3.0.1/bin/activate");
			$queue->addcommonds("bam_stat.py -i $outdir/alignment/$sn.bam > $outdir/alignment/$sn.bam_stat.txt");
			#$queue->addcommonds("read_distribution.py -i $outdir/alignment/$sn.bam -r $outdir/refer/annot_known.exon.bed > $outdir/alignment/$sn.read_distribution.txt");
			$queue->addcommonds("deactivate");
			$queue->run();
			$results{genome_alignment}++;
		}
		close FIN;
		
		$all_dependencies{"alignment.merge"}++;
		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("alignment.merge");
		$queue->set_job_disc("比对基因组结果统计分析");
		my @p = ();
		foreach my $sn(sort keys %sample_info){
			push @p, "alignment.sample.$sn";
		}
		$queue->set_job_depend(join(",",@p));
		$queue->set_file_check("$outdir/alignment/alignment.statistics.xls");
		$queue->set_work_dir($workdir);
		#$queue->addcommonds("$script/hisat2.stat.pl $outdir/alignment/ $outdir/alignment/alignment.hisat2.xls");
		$queue->addcommonds("$script/bam_stat.summary.pl $outdir/alignment/ $outdir/alignment/alignment.statistics.xls");
		$queue->run();
		$results{genome_alignment}++;
	}elsif(grep /^$method_exp$/,  qw/B1 B2 B3/){
		foreach my $sn(sort keys %sample_info){
			$all_dependencies{"alignment.sample_$sn"}++;
			$queue->set_job_cpu(10);
			$queue->set_job_mem(20);
			$queue->set_job_mode($job_mode);
			$queue->set_job_name("alignment.sample.$sn");
			$queue->set_job_disc("样本[$sn]比对基因组(star2)");
			$queue->set_job_depend("rawdata_quality_control.sample.$sn");
			$queue->set_file_check("$outdir/alignment/$sn.bam,$outdir/alignment/$sn.bam.bai");
			$queue->set_work_dir($workdir);
			
			$queue->addcommonds("module load star/2.7.2b");
			$queue->addcommonds("module load samtools/1.6");
			$queue->addcommonds("mkdir -p $outdir/alignment/$sn");
			$queue->addcommonds("star --genomeDir $outdir/refer/index_star --sjdbGTFfile $outdir/refer/annot_known.gtf --limitBAMsortRAM 40000000000 --runThreadN 10 --limitIObufferSize 500000000 --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --outFilterMultimapNmax 20 --outFilterMatchNminOverLread 0.66 --outFilterIntronMotifs None --outSJfilterReads All --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMstrandField intronMotif --outSAMattrRGline ID:$sn SM:$sn PL:illumina --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 0 --chimScoreDropMax 20 --chimScoreSeparation 10 --chimScoreJunctionNonGTAG -1 --quantMode TranscriptomeSAM --quantTranscriptomeBan IndelSoftclipSingleend --outReadsUnmapped Fastx --readFilesIn $outdir/quality_control/$sn/$sn.cleandata.R1.fq.gz $outdir/quality_control/$sn/$sn.cleandata.R2.fq.gz --readFilesCommand zcat --outFileNamePrefix $outdir/alignment/$sn/$sn. > $outdir/alignment/$sn.stdout 2> $outdir/alignment/$sn.stderr");
			$queue->addcommonds("mv $outdir/alignment/$sn/$sn.Aligned.sortedByCoord.out.bam $outdir/alignment/$sn.bam");
			$queue->addcommonds("samtools index $outdir/alignment/$sn.bam -\@ 10");
			$queue->addcommonds("rm -rf $outdir/alignment/$sn.sam");
			
			$queue->addcommonds("source /mnt/nfs/users/pub/software/python/3.7.5/env/RSeQC-3.0.1/bin/activate");
			$queue->addcommonds("bam_stat.py -i $outdir/alignment/$sn.bam > $outdir/alignment/$sn.bam_stat.txt");
			#$queue->addcommonds("read_distribution.py -i $outdir/alignment/$sn.bam -r $outdir/refer/annot_known.exon.bed > $outdir/alignment/$sn.read_distribution.txt");
			$queue->addcommonds("deactivate");
			$queue->run();
		}
		close FIN;
		
		$all_dependencies{"alignment.merge"}++;
		$queue->set_job_cpu(2);
		$queue->set_job_mem(5);
		$queue->set_job_mode($job_mode);
		$queue->set_job_name("alignment.merge");
		$queue->set_job_disc("比对基因组结果统计分析");
		my @p = ();
		foreach my $sn(sort keys %sample_info){
			push @p, "alignment.sample_$sn";
		}
		$queue->set_job_depend(join(",",@p));
		$queue->set_file_check("$outdir/alignment/alignment.statistics.xls");
		$queue->set_work_dir($workdir);
		$queue->addcommonds("$script/bam_stat.summary.pl $outdir/alignment/ $outdir/alignment/alignment.statistics.xls");
		$queue->run();
	}
}

sub explevel_analysis {
	if($method_exp eq "A1"){
		foreach my $xxx(@exp_class){
			foreach my $sn(sort keys %sample_info){
				my $bam = "$outdir/alignment/$sn.bam";
				$all_dependencies{"explevel_$xxx.sample.$sn"}++;
				$queue->set_job_cpu(8);
				$queue->set_job_mem(10);
				$queue->set_job_mode($job_mode);
				$queue->set_job_name("explevel_$xxx.sample.$sn");
				$queue->set_job_disc("样本[$sn]表达量评估[$xxx]");
				#$queue->set_job_depend("reference_preparation_$xxx.part_2,alignment.sample.$sn");
				$queue->set_job_depend("alignment.sample.$sn");
				$queue->set_file_check("$outdir/explevel_analysis/$xxx/ballgown/$sn/$sn.gtf");
				$queue->set_work_dir($workdir);
				$queue->addcommonds("module load stringtie/1.3.3b");
				$queue->addcommonds("mkdir -p $outdir/explevel_analysis/$xxx/ballgown/$sn");
				$queue->addcommonds("stringtie -B -e -p 8 -G $outdir/refer/annot_$xxx.gtf -o $outdir/explevel_analysis/$xxx/ballgown/$sn/$sn.gtf $outdir/alignment/$sn.bam");
				$queue->run();
				$results{explevel_analysis}++;
			}
			
			$all_dependencies{"explevel_$xxx.matrix"}++;
			$queue->set_job_cpu(2);
			$queue->set_job_mem(5);
			$queue->set_job_mode($job_mode);
			$queue->set_job_name("explevel_$xxx.matrix");
			$queue->set_job_disc("生成表达量矩阵[$xxx]");
			my @p;
			foreach my $sn(sort keys %sample_info){
				push @p, "explevel_$xxx.sample.$sn";
			}
			$queue->set_job_depend(join(",",@p));
			$queue->set_file_check("$outdir/explevel_analysis/$xxx/gene/matrix_gene.count.xls,$outdir/explevel_analysis/$xxx/gene/matrix_gene.fpkm.xls");
			$queue->set_work_dir($workdir);
			$queue->addcommonds("module load app/explevel_tools");
			$queue->addcommonds("module load app/table_tools");
		
			$queue->addcommonds("mkdir -p $outdir/explevel_analysis/$xxx");
			$queue->addcommonds("$script/prepDE.py --input $outdir/explevel_analysis/$xxx/ballgown -g $outdir/explevel_analysis/$xxx/matrix_gene.count.raw.csv -t $outdir/explevel_analysis/$xxx/matrix_transcript.count.raw.csv --length 150");
			
			foreach my $yyy(@exp_level){
				$queue->addcommonds("mkdir -p $outdir/explevel_analysis/$xxx/$yyy");
				$queue->addcommonds("sed 's/,/\\t/g' $outdir/explevel_analysis/$xxx/matrix_$yyy.count.raw.csv > $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.raw.txt");
				$queue->addcommonds("convert_count_to_explevel.pl --count $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.raw.txt --length $outdir/refer/annot_$xxx.$yyy"."2length --run_TMM --prefix $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy --type $yyy");
				# $queue->addcommonds("col_merge.pl --x $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.xls --y $outdir/refer/annot_known.$yyy"."2annot.xls --by_x $yyy\_id --by_y $yyy\_id --head_x --head_y --output $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.annot.xls");
				# $queue->addcommonds("col_merge.pl --x $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls  --y $outdir/refer/annot_known.$yyy"."2annot.xls --by_x $yyy\_id --by_y $yyy\_id --head_x --head_y --output $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.annot.xls ");
				#$queue->addcommonds("col_merge.pl --x $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.tpm.xls   --y $outdir/refer/annot_known.$yyy"."2annot.xls --by_x $yyy\_id --by_y $yyy\_id --head_x --head_y --output $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.tpm.annot.xls  ");
			}
			$queue->run();
			$results{explevel_analysis}++;
			
			foreach my $yyy(@exp_level){
				$all_dependencies{"explevel_$xxx\_$yyy.distribution"}++;
				$queue->set_job_cpu(2);
				$queue->set_job_mem(5);
				$queue->set_job_mode($job_mode);
				$queue->set_job_name("explevel_$xxx\_$yyy.distribution");
				$queue->set_job_disc("基于[$yyy]表达量的分布分析[$xxx]");
				$queue->set_job_depend("explevel_$xxx.matrix");
				$queue->set_file_check("$outdir/explevel_analysis/$xxx/$yyy/explevel_distribution");
				$queue->set_work_dir($workdir);
				$queue->addcommonds("module load app/explevel_tools");
				$queue->addcommonds("mkdir -p $outdir/explevel_analysis/$xxx/$yyy/explevel_distribution");
				$queue->addcommonds("plot_exp_distribution.pl --input $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls --type fpkm --outdir $outdir/explevel_analysis/$xxx/$yyy/explevel_distribution --prefix Sample");
				$queue->run();
				$results{explevel_distribution}++;
			}
			
			foreach my $yyy(@exp_level){
				next if scalar keys %sample_info <= 2;
				$all_dependencies{"explevel_$xxx\_$yyy.cluster"}++;
				$queue->set_job_cpu(2);
				$queue->set_job_mem(5);
				$queue->set_job_mode($job_mode);
				$queue->set_job_name("explevel_$xxx\_$yyy.cluster");
				$queue->set_job_disc("基于[$yyy]表达量的聚类分析[$xxx]");
				$queue->set_job_depend("explevel_$xxx.matrix");
				$queue->set_file_check("$outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.xls,$outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls");
				$queue->set_work_dir($workdir);
				$queue->addcommonds("module load app/explevel_tools");
				$queue->addcommonds("mkdir -p $outdir/explevel_analysis/$xxx/$yyy/explevel_cluster");
				$queue->addcommonds("plot_exp_pca.pl -i $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls -o $outdir/explevel_analysis/$xxx/$yyy/explevel_cluster -pre Sample");
				foreach my $class(sort keys %meta_grp2sn){
					next if $class eq "Sample";
					$queue->addcommonds("plot_exp_pca.pl -i $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls -o $outdir/explevel_analysis/$xxx/$yyy/explevel_cluster -m $file_meta -g $class -pre $class");
				}
				$queue->run();
				$results{explevel_cluster}++;
			}
		}
		
	}elsif($method_exp eq "A2"){
		
	}elsif($method_exp eq "A3"){
		
	}
}

sub expldiff_analysis{
	return if scalar keys %sample_info <= 1;
	foreach my $std(sort keys %comp_grp2sn){
		next unless -e "$outdir/group_info/$std.cts";
		my $class = [split /\./,$std]->[0];
		if($method_exp eq "A1"){
			foreach my $xxx(@exp_class){
				foreach my $yyy(@exp_level){
					$all_dependencies{"expldiff_$xxx\_$yyy.$std"}++;
					$queue->set_job_cpu(15);
					$queue->set_job_mem(30);
					$queue->set_job_mode($job_mode);
					$queue->set_job_name("expldiff_$xxx\_$yyy.$std");
					$queue->set_job_disc("差异表达$yyy\[$xxx]鉴定与分析[$std]");
					$queue->set_job_depend("explevel_$xxx.matrix");
					$queue->set_file_check("$outdir/expldiff_analysis/$xxx/$yyy/$std/$std.de$yyy\_all.xls,$outdir/expldiff_analysis/$xxx/$yyy/$std/plot/$std.ScatterPlot.pdf,$outdir/expldiff_analysis/$xxx/$yyy/$std/heatmap/$std.HeatmapPlot.pdf");
					$queue->set_work_dir($workdir);
					
					my $kes = "";
					my $add_info = "";
					if($kegg_filter){
						$add_info .= " --remove_human";
					}
					if($ppi_denovo){
						$add_info .= " --ppi_denovo --ppi_prot $outdir/refer/annot_$xxx.protein.fa --ppi_prot2seq $outdir/refer/annot_$xxx.protein2$yyy";
					}
					if($xxx eq "known"){
						$kes = $kegg_species;
					}else{
						$kes = "ko";
					}
					$queue->addcommonds("module load app/expldiff_tools");
					$queue->addcommonds("mkdir -p $outdir/expldiff_analysis/$xxx/$yyy/$std");
					$queue->addcommonds("echo -e '$comp_param{$std}{method}\\t$comp_param{$std}{pval}\\t$comp_param{$std}{padj}\\t$comp_param{$std}{fc}' > $outdir/expldiff_analysis/$xxx/$yyy/$std/param.txt");
					
					$queue->addcommonds("run_deg_analysis.pl --count $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.xls --norm $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls --class $yyy --type fpkm --method $comp_param{$std}{method} --pval $comp_param{$std}{pval} --padj $comp_param{$std}{padj} --fc $comp_param{$std}{fc} --group $outdir/group_info/$std.txt --contrasts $outdir/group_info/$std.cts --outdir $outdir/expldiff_analysis/$xxx/$yyy/$std --prefix $class");
					
					#$queue->addcommonds("run_deg_analysis.pl --count $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.count.xls --norm $outdir/explevel_analysis/$xxx/$yyy/matrix_$yyy.fpkm.xls --class $yyy --type fpkm --method $comp_param{$std}{method} --pval $comp_param{$std}{pval} --padj $comp_param{$std}{padj} --fc $comp_param{$std}{fc} --group $outdir/group_info/$std.txt --contrasts $outdir/group_info/$std.cts --outdir $outdir/expldiff_analysis/$xxx/$yyy/$std --prefix $class --id2annot $outdir/refer/annot_$xxx.$yyy"."2annot.xls --id2go $outdir/refer/annot_$xxx.$yyy"."2go.txt --id2kegg $outdir/refer/annot_$xxx.$yyy"."2kegg.txt --rm_low_exp --run_volcano --run_scatter --run_maplot --run_heatmap --run_go --run_kegg --run_ipath --run_gsea --run_ppi --species $kes --ppi_taxonomy $ppi_taxonomy $add_info");
					$queue->run();
					$results{expldiff_analysis}++;
				}
			}
		}
	}
}




