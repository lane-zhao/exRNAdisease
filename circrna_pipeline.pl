#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use lib "$RealBin/PerlLib/";
use PBS::Queue;

my $BEGIN_TIME=time();
my $script   = "$Bin/script";
my $software = "$Bin/software";
my $database = "$Bin/database"; 

my $workdir = getcwd;

my %species_info;
&load_species_info("$Bin/species.config.ini", \%species_info);


my ($help, $sample_list, $outdir, $species, $show, $refer_seq, $refer_gtf, $index_bwa, $index_star, $mirna_seq, $gene2miRNA, $gene2go, $gene2kegg, $cpu, $step, $rerun, $group, $gene2name);
my $param_strand = "";
$step = "";

GetOptions(
	"help!"        => \$help,
	"input:s"      => \$sample_list,
	"outdir:s"     => \$outdir,
	"species:s"    => \$species,
	"group:s"      => \$group,
	"list!"        => \$show,
	
	"refer_seq:s"  => \$refer_seq,
	"refer_gtf:s"  => \$refer_gtf,
	"index_bwa:s"  => \$index_bwa,
	"index_star:s" => \$index_star,
	
	"gene2name:s"  => \$gene2name,
	"gene2go:s"    => \$gene2go,
	"gene2kegg:s"  => \$gene2kegg,
	
	"mirna_seq:s"  => \$mirna_seq,
	"gene2miRNA:s" => \$gene2miRNA,
	
	"cpu:i"    => \$cpu,
	
	"resun"    => \$rerun,
);
	
&usage if ($help);

if($species){
	if(exists $species_info{$species}){
		$refer_seq  = $species_info{$species}{refer_seq};
		$refer_gtf  = $species_info{$species}{refer_gtf};
		$index_bwa  = $species_info{$species}{index_bwa};
		$index_star = $species_info{$species}{index_star};
		
		$gene2name  = $species_info{$species}{gene2name};
		$gene2go    = $species_info{$species}{gene2go};
		$gene2kegg  = $species_info{$species}{gene2kegg};
		
		$gene2miRNA = $species_info{$species}{gene2miRNA};
		$mirna_seq  = $species_info{$species}{mirna_seq};
	}
}

## initialization ##
unless(defined( $sample_list )){
	print STDERR "Error:-input must be specified!\n";
	&usage;
}
unless(defined( $refer_seq )){
	print STDERR "Error:-refer_seq must be specified!\n";
	&usage;
}
unless(defined( $refer_gtf )){
	print STDERR "Error:-refer_gtf must be specified!\n";
	&usage;
}

$cpu ||= 4;
$outdir ||= "result";
$outdir = abs_path( $outdir );
system("mkdir -p $outdir");

unless(defined $group){
	open FOUT,">$outdir/group.info";
	print FOUT "id\tgroup\n";
	open FIN,"$sample_list";
	while(<FIN>){
		chomp;
		my ($a) = split /\t/,$_;
		print FOUT $a,"\t",$a,"\n";
	}
	close FOUT;
	close FIN;
	$group = "$outdir/group.info";
}

$group      = abs_path( $group      ) if ( defined( $group      ) );

#$refer_seq  = abs_path( $refer_seq  ) if ( defined( $refer_seq  ) );
#$refer_gtf  = abs_path( $refer_gtf  ) if ( defined( $refer_gtf  ) );
#$index_bwa  = abs_path( $index_bwa  ) if ( defined( $index_bwa  ) );
#$index_star = abs_path( $index_star ) if ( defined( $index_star ) );

$sample_list = abs_path( $sample_list) if ( defined( $sample_list) );

$gene2name  = abs_path( $gene2name  ) if ( defined( $gene2name  ) );
$gene2go    = abs_path( $gene2go    ) if ( defined( $gene2go    ) );
$gene2kegg  = abs_path( $gene2kegg  ) if ( defined( $gene2kegg  ) );

$gene2miRNA = abs_path( $gene2miRNA ) if ( defined( $gene2miRNA ) );
$mirna_seq  = abs_path( $mirna_seq  ) if ( defined( $mirna_seq  ) );


$rerun ||= 0;

my $queue = PBS::Queue->new(
	{
		'cluster_queue'  => "zh",
		'pbs_queue_name' => "circ",
		'max_threads'    => 200,
		'max_jobs'       => 20,
		'run_local'      => 1,
		'runmod'         => $rerun,
	}
);

unless(defined($index_bwa)){
	&refer_index_bwa();
	$index_bwa = "$outdir/refer/bwa/genome.fa";
}
unless(defined($index_star)){
	&refer_index_star();
	$index_star = "$outdir/refer/star";
}

my %samples;
&fastq_qc();
&identification();
&explevel();
&diffexp();
&diffcluster();
&diffhost();
$queue->jointhreads();

print " all jobs done!\n";

###########################################################################

sub refer_index_bwa {
	$queue->set_ppn(2);
	$queue->set_memory('5G');
	$queue->set_pbs_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/refer/bwa");
	$queue->addcommond("cd $outdir/refer/bwa");
	$queue->addcommond("ln -s $refer_seq genome.fa");
	$queue->addcommond("bwa index genome.fa");
	$queue->addcommond("cd $workdir");
	$queue->qsub();
}

sub refer_index_star {
	$queue->set_ppn(2);
	$queue->set_memory('5G');
	$queue->set_pbs_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/refer/star");
	$queue->addcommond("cd $outdir/refer/stat");
	$queue->addcommond("ln -s $refer_seq genome.fa");
	$queue->addcommond("ln -s $refer_gtf genome.gtf");
	$queue->addcommond("STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles genome.fa --sjdbGTFfile genome.gtf --sjdbOverhang 150");
	$queue->addcommond("cd $workdir");
	$queue->qsub();
}

sub fastq_qc{
	$queue->set_ppn(1);
	$queue->set_memory('5G');
	$queue->set_pbs_work_dir($workdir);
	
	$queue->addcommond("fastq_qc.pl -i $sample_list -o $outdir/cleandata -n 6 -z");
	$queue->get_log()->write("原始数据质控分析开始");
	$queue->qsub();
	$queue->wait();
}

sub identification {
	open FIN,"$outdir/cleandata/cleandata.fq.list";
	while(<FIN>){
		chomp;
		my ($sample,$type,$fq1,$fq2) = split /\t/,$_;
		$samples{$sample}++;
		my @p = split /\t/,$_;
		if(scalar @p == 4){
			my ($sample,$type,$fq1,$fq2) = @p;
			$queue->set_ppn(15);
			$queue->set_memory('20G');
			$queue->set_pbs_work_dir($workdir);
			$queue->addcommond("mkdir -p $outdir/identification/CIRCexplorer/$sample");
			$queue->addcommond("STAR --runThreadN 20 --genomeDir $index_star --readFilesIn $fq1 $fq2 --outFileNamePrefix $outdir/identification/CIRCexplorer/$sample/ --chimSegmentMin 10");
			$queue->addcommond("star_parse.py $outdir/identification/CIRCexplorer/$sample/Chimeric.out.junction $outdir/identification/CIRCexplorer/$sample/fusion_junction.txt");
			$queue->addcommond("sh $script/gtf2genePred.sh $refer_gtf $outdir/identification/CIRCexplorer");
			$queue->addcommond("CIRCexplorer.py -j $outdir/identification/CIRCexplorer/$sample/fusion_junction.txt -g $refer_seq -r $outdir/identification/CIRCexplorer/genome.genePred -o $outdir/identification/CIRCexplorer/$sample");
			$queue->get_log()->write("\[$sample\] CIRCexplorer分析开始");
			$queue->qsub();
			
			$queue->set_ppn(15);
			$queue->set_memory('20G');
			$queue->set_pbs_work_dir($workdir);
			$queue->addcommond("mkdir -p $outdir/identification/CIRI");
			$queue->addcommond("bwa mem -M -t 20 $index_bwa $fq1 $fq2 > $outdir/identification/CIRI/$sample.sam 2> $outdir/identification/CIRI/$sample.bwa.log");
			$queue->addcommond("CIRI2.pl -I $outdir/identification/CIRI/$sample.sam -O $outdir/identification/CIRI/$sample.ciri -F $refer_seq -T 20 -A $refer_gtf -G $outdir/identification/CIRI/$sample.ciri.log");
			$queue->get_log()->write("\[$sample\] CIRI分析开始");
			$queue->qsub();
		}elsif(scalar @p == 3){
			my ($sample,$type,$fq1) = @p;
			die "not support SE now!!";
		}
		
	}
	close FIN;
	$queue->wait();
	
	foreach my $sample(sort keys %samples){
		$queue->set_ppn(2);
		$queue->set_memory('5G');
		$queue->set_pbs_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/identification/Intersect/");
		$queue->addcommond("ln -s $outdir/identification/CIRCexplorer/$sample\_circ.txt $outdir/identification/Intersect/$sample.circ.xls");
		$queue->addcommond("ln -s $outdir/identification/CIRI/$sample.ciri $outdir/identification/Intersect/$sample.ciri.xls");
		$queue->addcommond("cut -f 1,2,3,4 $outdir/identification/Intersect/$sample.circ.xls > $outdir/identification/Intersect/$sample.circ");
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.ciri.xls | awk -F\"\\t\" 'OFS=\"\\t\"{print \$2,\$3,\$4,\$1}' > $outdir/identification/Intersect/$sample.ciri");
		$queue->addcommond("bedtools intersect -a $outdir/identification/Intersect/$sample.ciri -b $outdir/identification/Intersect/$sample.circ -wo | awk 'OFS=\"\\t\"{print \$4,\$3-\$2-\$9}'| awk '{if(\$2==0){print \$1}}' > $outdir/identification/Intersect/$sample.circRNA.list");
		$queue->addcommond("row_select.pl -input $outdir/identification/Intersect/$sample.ciri.xls -list $outdir/identification/Intersect/$sample.circRNA.list -col 1 -skip 1 -output $outdir/identification/Intersect/$sample.circRNA.txt");
		$queue->addcommond("$script/CIRI_format.pl $outdir/identification/Intersect/$sample.circRNA.txt $outdir/identification/Intersect/$sample.circRNA.xls $gene2name");
		
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.circRNA.xls | cut -f 9 | sort | uniq -c | sed 's/^[ ][ ]*//g' | awk '{print \$2\"\\t\"\$1}' > $outdir/identification/Intersect/$sample.type.tmp_1");
		$queue->addcommond("echo -e \"type\\tnumber\" > $outdir/identification/Intersect/$sample.type.tmp_2");
		$queue->addcommond("cat $outdir/identification/Intersect/$sample.type.tmp_2 $outdir/identification/Intersect/$sample.type.tmp_1 > $outdir/identification/Intersect/$sample.type.xls");
		$queue->addcommond("rm $outdir/identification/Intersect/$sample.type.tmp_2 $outdir/identification/Intersect/$sample.type.tmp_1");
		$queue->addcommond("Highchart.pl -type pie -t $outdir/identification/Intersect/$sample.type.xls -width 900 -piemark \"u:Intergenic transcript\" -title \"circRNA types distribution\" -pdfcompress no -color_the google");
		$queue->addcommond("mv $outdir/identification/Intersect/$sample.type.xls.pdf $outdir/identification/Intersect/$sample.type.pdf");
		$queue->addcommond("rm $outdir/identification/Intersect/$sample.type.xls.html");
		
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.circRNA.xls | cut -f 2 | sort | uniq -c | sed 's/^[ ][ ]*//g' | awk '{print \$2\"\\t\"\$1}' > $outdir/identification/Intersect/$sample.chr.tmp_1");
		$queue->addcommond("echo -e \"chromosome\\tcircRNA_number\" > $outdir/identification/Intersect/$sample.chr.tmp_2");
		$queue->addcommond("cat $outdir/identification/Intersect/$sample.chr.tmp_2 $outdir/identification/Intersect/$sample.chr.tmp_1 > $outdir/identification/Intersect/$sample.chr.xls");
		$queue->addcommond("rm $outdir/identification/Intersect/$sample.chr.tmp_1 $outdir/identification/Intersect/$sample.chr.tmp_2");
		
		$queue->addcommond("Rscript $script/chr_distribution.R $outdir/identification/Intersect/$sample.chr.xls $outdir/identification/Intersect/$sample.chr.pdf");
		
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.circRNA.xls | awk '{print \$2\"\\t\"\$3-1\"\\t\"\$4\"\\t\"\$1\"\\t\"0\"\\t\"\$11}' > $outdir/identification/Intersect/$sample.circRNA.bed");
		$queue->addcommond("bedtools getfasta -fi $refer_seq -bed $outdir/identification/Intersect/$sample.circRNA.bed -name -s -fo $outdir/identification/Intersect/$sample.circRNA.fa");
		
		$queue->get_log()->write("\[$sample\] circRNA鉴定分析开始");
		$queue->qsub();
	}
	$queue->wait();
	
	$queue->set_ppn(2);
	$queue->set_memory('10G');
	$queue->set_pbs_work_dir($workdir);
	$queue->addcommond("cd $outdir/identification/Intersect/");
	foreach my $sample(sort keys %samples){
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.circRNA.xls | awk '{print \$2\"\\t\"\$3-1\"\\t\"\$4\"\\t\"\$1\"\\t\"0\"\\t\"\$11}' > $sample.bed");
		$queue->addcommond("head -n 1 $outdir/identification/Intersect/$sample.circRNA.xls > head");
		$queue->addcommond("sed 1d $outdir/identification/Intersect/$sample.circRNA.xls > $sample.temp");
		$queue->addcommond("");
	}
	$queue->addcommond("cat *.bed | sort -u | bedtools sort -i - > circRNA.region");
	$queue->addcommond("cat *.temp | sort -u | cat head - > circRNA.xls");
	$queue->addcommond("bedtools getfasta -fi $refer_seq -bed circRNA.region -name -s -fo circRNA.fa");
	$queue->addcommond("rm *.bed *.temp circRNA.region head");
	
	$queue->get_log()->write("cirRNA序列提取分析开始");
	$queue->qsub();
	
	
	$queue->set_ppn(2);
	$queue->set_memory('5G');
	$queue->set_pbs_work_dir($workdir);
	$queue->addcommond("$script/star.stat.pl $outdir/identification/CIRCexplorer/ $outdir/identification/alignment.statistics.xls");
	$queue->get_log()->write("序列比对结果统计开始");
	$queue->qsub();
	$queue->wait();
}


sub explevel {
	foreach my $sample(sort keys %samples){
		$queue->set_ppn(2);
		$queue->set_memory('10G');
		$queue->set_pbs_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/explevel/");
		$queue->addcommond("cd $outdir/explevel/");
		$queue->addcommond("$script/circRNA_exp_sample.pl $outdir/cleandata/cleandata.fq.list $outdir/identification/alignment.statistics.xls $sample $outdir/identification/Intersect/$sample.circRNA.txt");
		$queue->get_log()->write("[$sample] 表达量评估开始");
		$queue->qsub();
	}
	$queue->wait();
	
	
	open F1,">$outdir/explevel/circRNA.count.list";
	open F2,">$outdir/explevel/circRNA.srpbm.list";
	foreach my $sample(sort keys %samples){
		print F1 "$sample\t$outdir/explevel/$sample.count.xls\n";
		print F2 "$sample\t$outdir/explevel/$sample.srpbm.xls\n";
	}
	close F1;
	close F2;
	
	$queue->set_ppn(1);
	$queue->set_memory('2G');
	$queue->set_pbs_work_dir($workdir);
	$queue->addcommond("$script/circRNA_exp_merge.pl $outdir/explevel/circRNA.count.list $outdir/explevel/circRNA.count.xls");
	$queue->addcommond("$script/circRNA_exp_merge.pl $outdir/explevel/circRNA.srpbm.list $outdir/explevel/circRNA.srpbm.xls");
	$queue->get_log()->write("表达矩阵生成开始");
	$queue->qsub();
	$queue->wait();
	
	if( scalar(keys %samples) > 1){
		$queue->set_ppn(1);
		$queue->set_memory('2G');
		$queue->set_pbs_work_dir($workdir);
		$queue->addcommond("$script/plot_sample2sample.pl -i $outdir/explevel/circRNA.count.xls -c $outdir/explevel/circRNA.srpbm.xls -o $outdir/explevel/ -m $group -g group -pre sample");
		$queue->get_log()->write("表达聚类分析开始");
		$queue->qsub();
		$queue->wait();
	}
}


sub diffexp {
	if( scalar(keys %samples) > 1){
		$queue->set_ppn(2);
		$queue->set_memory('10G');
		$queue->set_pbs_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/diffexp/");
		$queue->addcommond("cd $outdir/diffexp/");
		$queue->addcommond("sed 1d $group | awk '{print \$2\"\\t\"\$1}' > $outdir/diffexp/group.txt");
		if($group && -e $group){
			$queue->addcommond("run_DE_analysis.pl --count $outdir/explevel/circRNA.count.xls --fpkm $outdir/explevel/circRNA.srpbm.xls --method edgeR --group $outdir/diffexp/group.txt --num_parallel 2 --pval 0.05 --log2FC 1 --output $outdir/diffexp/");
		}else{
			$queue->addcommond("run_DE_analysis.pl --count $outdir/explevel/circRNA.count.xls --fpkm $outdir/explevel/circRNA.srpbm.xls --method edgeR --num_parallel 2 --pval 0.05 --log2FC 1 --output $outdir/diffexp/");
		}
		
		$queue->get_log()->write("差异表达分析开始");
		$queue->qsub();
	}
}

sub diffcluster {
	if( scalar(keys %samples) > 1){
		$queue->set_ppn(2);
		$queue->set_memory('10G');
		$queue->set_pbs_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/diffcluster/");
		$queue->addcommond("cd $outdir/diffcluster/");
		$queue->addcommond("cat $outdir/diffexp/*.sde_all.list | sort -u > $outdir/diffcluster/total.de.list");
		$queue->addcommond("row_select.pl -input $outdir/explevel/circRNA.srpbm.xls -list $outdir/diffcluster/total.de.list -col 1 -skip 1 -output $outdir/diffcluster/total.de.matrix");
		$queue->addcommond("Rscript $script/heatmap.r $outdir/diffcluster/total.de.matrix $outdir/diffcluster/");
		$queue->get_log()->write("差异cirRNA聚类分析开始");
		$queue->qsub();
	}
}

sub diffhost {
	if( scalar(keys %samples) > 1){
		foreach my $list(glob("$outdir/diffexp/*gene.sde_all.list")){
			(my $name = basename($list)) =~ s/.gene.sde_all.list//;
			$queue->set_ppn(2);
			$queue->set_memory('10G');
			$queue->set_pbs_work_dir($workdir);
			$queue->addcommond("mkdir -p $outdir/diffhost/");
			$queue->addcommond("cd $outdir/diffhost/");
			$queue->addcommond("row_select.pl -i $outdir/sequence/circRNA.xls -l $list -skip 1 | cut -f 11 | sed 's/;/\\n/g' | sort -u >  $outdir/diffhost/$name.host_gene.list");
			$queue->addcommond("kegg_enrich_human.pl --de $outdir/diffhost/$name.host_gene.list");
			$queue->addcommond("go_enrich_human.pl   --de $outdir/diffhost/$name.host_gene.list");
			
		}
		$queue->get_log()->write("cirRNA序列提取分析开始");
		$queue->qsub();
	}
}



sub load_species_info{
	my ($configfile,$hash) = @_;
	
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

sub usage{
	die(qq/
Usage: $0 [options] <arguments>
    <clean_data_fq.lst> is a basic input file: absolute path containing the clean data with fastq type(eg: trimPairFq.list)
    <out_dir> all result files will be placed in this direction.
Arguments:
  Input Options:
    -input          <string>    clean data fq file list. [ required ]
    -group          <string>    group info [ optional ]
    -species        <string>    species name, support list could be check by list arguments. [ optional ]
    -list           <NA>        show the support species list	
    -outdir         <string>    output dir, default is "result".
    -refer_seq      <string>    reference sequence in fast fromat [ required ]
    -refer_gtf      <string>    reference annotation in gtf format [ required ]
    -index_bwa      <string>    bwa index of reference [ optional ]
    -index_star     <string>    star index of reference [ optional ]
  
    -gene2name      <string>    gene id to gene name conversion [ optional ]
    -gene2go        <string>    go annotation file [ optional ]
    -gene2kegg      <string>    kegg annotation file [ optional ]

    -mirna_seq      <string>    miRNA sequence file in fasta format [ optional ]
    -gene2miRNA     <string>    gene and miRNA interaction information [ optional ]
	
    -cpu            <string>    thread number, default is 4
    -rerun          <NA>        rerun whole pipeline
    -only_cmd       <NA>        only print cmd
    -help           <NA>        output help information to screen

###################
#PS1:<clean_data_fq.lst> contain 3 or 4 colums, sample_name SE(or PE) read1.fq read2.fq
#    Sample_1	PE	Sample_1_trim1.fq	Sample_1_trim2.fq	
#    Sample_2	PE	Sample_2_trim1.fq	Sample_2_trim2.fq
#    ......
#    For different lines, samples name must be unique !!!
#PS2:[group] contain 2 colums, sample_name replicated_samples_name
#    id	group
#    Sample_1	Group_1
#    Sample_2	Group_1
#    Sample_3	Group_1
#    Sample_4	Group_2
#    Sample_5	Group_2
#    Sample_6	Group_2
#    ......
#    For each group, samples are biological replicates !!!
###################		
\n/);
}	

