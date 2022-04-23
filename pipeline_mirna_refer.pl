#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Spec::Functions qw(rel2abs);
use PBS::Queue;

my ($input, $outdir, $refer, $refer_index, $library, $group_file, $rfam_seed, $miRbase, $hairpin, $species, $list_species, $plant, $plant_db, $animal, $animal_db, $target, $anno, $class, $go,  $kegg, $gene, $help);
my ($max_cpus, $max_jobs, $queue_name, $job_prefix, $job_local, $queue_rerun);
my (%samples, %groups);
my $sample_num;
my $group_num;
my $script   = "$Bin/script";
my $database = "$Bin/database"; 
my $software = "$Bin/software"; 
my $workdir  = getcwd;

GetOptions(
	"help!"         => \$help,
	"input=s"       => \$input,
	"outdir=s"      => \$outdir,
	"refer:s"       => \$refer,
	"refer_index:s" => \$refer_index,
	"group:s"       => \$group_file,
	"library:s"     => \$library,
	
	"rfam_seed:s"   => \$rfam_seed,
	"miRbase:s"     => \$miRbase,
	"hairpin:s"     => \$hairpin,
	"species:s"     => \$species,
	"list!"         => \$list_species,
	"class:s"       => \$class,

	"animal!"       => \$animal,
	"animal_db:s"   => \$animal_db,
	"plant!"        => \$plant,
	"plant_db:s"    => \$plant_db,
	
	"target:s"      => \$target,
	"anno:s"        => \$anno,
	"go:s"          => \$go,
	"kegg:s"        => \$kegg,
	"gene:s"        => \$gene,
	
	"max_cpus:i"    => \$max_cpus,
	"max_jobs:i"    => \$max_jobs,
	"queue_name:s"  => \$queue_name,
	"rerun!"        => \$queue_rerun,
	"job_prefix:s"  => \$job_prefix,
	"local!"        => \$job_local,
);


&Parsing_Parameters();

my $queue = PBS::Queue->new(
	{
		'queue_name'  => $queue_name,
		'job_prefix'  => $job_prefix,
		'max_cpus'    => $max_cpus,
		'max_jobs'    => $max_jobs,
		'job_local'   => $job_local,
		'queue_rerun' => $queue_rerun,
	}
);

&Load_Sample_Info;
&Load_Group_Info;
&Fastq_QC();
&Uniq();
&Align_to_rfam();
&Bais();
&Mapping2genome();
&Expression();
&Family();
&Venn();
&DiffExpression();
&DiffCluster();
&DiffTarget();
&DiffTargetEnrichment();
&result_collection();


################################################################################

sub Parsing_Parameters{
	&usage if ($help);
	&list_species if($list_species);
	
	unless(defined( $input )){
		print STDERR "Error: -input must be specified!\n";
		&usage;
	}
	unless(defined( $refer )){
	print STDERR "Error: -refer must be specified!\n";
		&usage;
	}

	$class     ||= "animal";
	
	$library   ||= "true-seq";
	
	$rfam_seed ||= "$Bin/database/Rfam/Rfam.seed";
	$miRbase   ||= "$Bin/database/miRBase/mature.fa";
	$hairpin   ||= "$Bin/database/miRBase/hairpin.fa";
	$plant_db  ||= "$Bin/database/miRBase/Plant";
	$animal_db ||= "$Bin/database/miRBase/Animal";
	$outdir    ||= "out";
	
	$queue_name  ||= "yx";
	$job_prefix  ||= "miRNA";
	$queue_rerun ||= 0;
	$job_local   ||= 1;
	$max_cpus    ||= 200;
	$max_jobs    ||= 30;
	
	$outdir = abs_path( $outdir );
	system("mkdir -p $outdir");
}

sub Load_Sample_Info{
	open  FIN,$input or die "can not open raw fq list file!";
	while (<FIN>){
		chomp;
		my @s = split /\s+/,$_;
		next if @s != 3;
		if(exists $samples{$s[0]}){
			die "sample_name $s[0] is duplicate! please check input file[$input]!";
		}
		$samples{$s[0]}{data} = $s[1];
		$samples{$s[0]}{abbr} = $s[2];
	}
}

sub Load_Group_Info{
	unless(defined $group_file){
		open FOUT,">$outdir/group.tmp";
		print FOUT "id\tgroup\n";
		foreach my $sn(sort keys %samples){
			print $sn,"\t",$sn,"\n";
		}
		close FOUT;
	}else{
		open  FIN,"<$group_file" or die "can not open group file [$group_file]!";
		my $head = <FIN>; chomp($head);
		unless($head eq "id\tgroup"){
			die "group file format error!";
		}
		while(<FIN>){
			chomp;
			my ($sample,$group) = split /\t/,$_;
			unless(exists $samples{$sample}){
				die "sample_name[$sample] in group_file does not existe in input!";
			}
			$groups{$group}{$sample}++;
		}
		close FIN;
		
		my %hash;
		open FOUT,">$outdir/group.info";
		open F2,">$outdir/group.tmp";
		print F2 "id\tgroup\n";
		foreach my $group(sort keys %groups){
			foreach my $sample(sort keys %{$groups{$group}}){
				print FOUT $group,"\t",$sample,"\n";
				print F2 $sample,"\t",$group,"\n";
				$hash{$sample}++;
				if($hash{$sample} > 1){
					die "In this version, one sample only belongs to one group!!";
				}
			}
		}
		close FOUT;
		close F2;
		$group_file = "$outdir/group.info";
	}
	
}

sub Fastq_QC{
	foreach my $sn(sort keys %samples){
		my $sn_dir = "$outdir/fastq_QC/$sn";
		my $sn_fq  = "$sn_dir/$sn.fq";
		$queue->set_job_cpu(2);
		$queue->set_job_mem(2);
		$queue->set_job_name("样本 $sn 数据质量评估与剪切");
		$queue->set_work_dir($workdir);
		$queue->addcommond("mkdir -p $sn_dir");
		my @fqs = split /,/,$samples{$sn}{data};
		$queue->addcommond("rm -f $sn_dir/$sn.fq");
		foreach my $fq(@fqs){
			if($fq =~ /.gz$/){
				$queue->addcommond("zcat $fq >> $sn_dir/$sn.fq");
			}else{
				$queue->addcommond("cat  $fq >> $sn_dir/$sn.fq");
			}
		}
		$queue->addcommond("$script/qc/trim_fastq.pl -input $sn_dir/$sn.fq -out $sn_dir/$sn.cut.fq -type cut -start 1 -len 51");
		$queue->addcommond("$software/fastx_toolkit/0.0.14/bin/fastx_quality_stats -i $sn_dir/$sn.cut.fq -o $sn_dir/$sn.stat -Q 33");
		$queue->addcommond("$software/fastx_toolkit/0.0.14/bin/fastq_quality_boxplot_graph.sh -i $sn_dir/$sn.stat -o $sn_dir/$sn.stat.qual.ps -p");
		$queue->addcommond("ps2pdf $sn_dir/$sn.stat.qual.ps $sn_dir/$sn.stat.qual.pdf ");
		$queue->addcommond("echo -e '$sn\\t$sn_dir/$sn.fq\\n' > $sn_dir/$sn.rawfq.list");
		$queue->addcommond("java -jar $script/qa/FastqStat.jar -i $sn_dir/$sn.rawfq.list > $sn_dir/$sn.rawdata.stat");
		$queue->addcommond("cut -f 1-3,11-14 $sn_dir/$sn.rawdata.stat > $sn_dir/$sn.rawdatashow.xls");
		
		if($library eq "exosome"){
			#$queue->addcommond("cutadapt -n 1 -q 20 -e 0.1 -m 18 -a AAAAAAAAAA -o $sn_dir/$sn.trim.fq $sn_dir/$sn.cut.fq &> $sn_dir/$sn.cutadapt.log");
			$queue->addcommond("$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $sn_dir/$sn.cut.fq -o $sn_dir/$sn.filtered.fq -a AAAAAAAAAA -Q 33 -l 18 -v > $sn_dir/$sn.qc.txt");
			
		}elsif($library eq "true-seq"){
			#$queue->addcommond("cutadapt -n 1 -q 20 -e 0.1 -m 18 -a TGGAATTCTCGGGTGCCAAGG -o $sn_dir/$sn.trim.fq $sn_dir/$sn.cut.fq &> $sn_dir/$sn.cutadapt.log");
			$queue->addcommond("$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $sn_dir/$sn.cut.fq -o $sn_dir/$sn.filtered.fq -a TGGAATTCTCGGGTGCCAAGG -Q 33 -l 18 -v > $sn_dir/$sn.qc.txt");
		}elsif($library eq "VAHTS"){
			#$queue->addcommond("cutadapt -n 1 -q 20 -e 0.1 -m 18 -a AGATCGGAAGAGCAC -g TACAGTCCGACGATC -o $sn_dir/$sn.trim.fq $sn_dir/$sn.cut.fq &> $sn_dir/$sn.cutadapt.log");
			$queue->addcommond("$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $sn_dir/$sn.cut.fq -o $sn_dir/$sn.filtered.fq -a AGATCGGAAGAGCAC -Q 33 -l 18 -v > $sn_dir/$sn.qc.txt");
		}
		$queue->addcommond("cd $sn_dir");
		$queue->addcommond("$script/qc/DynamicTrim.pl $sn.filtered.fq -h 20 -bwa");
		$queue->addcommond("$script/qc/fastq_to_fasta.pl -i $sn.filtered.fq.trimmed -o $sn.filtered.trimmed.fa");
		$queue->addcommond("$script/qc/trim_seq.pl -f $sn.filtered.trimmed.fa -o $sn.filtered.trimmed.len.fa -min 18 -max 32 >> $sn.qc.txt");
		$queue->addcommond("$script/qc/fasta_clipping_histogram_unique.pl $sn.filtered.trimmed.len.fa $sn.length_distribution.png > $sn.length_distribution.xls");
		# $queue->addcommond("$script/qc/fastq_to_fasta.pl -i $sn.trim.fq -o $sn.trim.fa");
		# $queue->addcommond("$script/qc/trim_seq.pl -f $sn.trim.fa -o $sn.trim.len.fa -min 18 -max 32 >> $sn.qc.txt");
		# $queue->addcommond("$script/qc/fasta_clipping_histogram_unique.pl $sn.trim.len.fa $sn.length_distribution.png > $sn.length_distribution.xls");
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("所有数据质量评估整合");
	$queue->set_work_dir($workdir);
	
	$queue->addcommond("cd $outdir/fastq_QC");
	$queue->addcommond("rm -rf qc.lst");

	foreach my $sn(sort keys %samples){
		$queue->addcommond("echo $outdir/fastq_QC/$sn/$sn.qc.txt >> qc.lst");
	}
	$queue->addcommond("$script/qc/get_miRNA_fq_data_table.pl qc.lst > QC_stat.xls");
	
	my $cmd = "$script/public/cat_files.pl -head 1 -outfile $outdir/fastq_QC/rawdatashow.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn/$sn.rawdatashow.xls";
	}
	$queue->addcommond("$cmd");
	$queue->run();
	$queue->wait();
}

sub Uniq{
	$queue->set_job_cpu(2);
	$queue->set_job_mem(10);
	$queue->set_job_name("选取uniq读段");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Uniq");
	my $cmd = "$script/uniq/get_uniqueFastas_cfg.pl $outdir/Uniq";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn/$sn.filtered.trimmed.len.fa $sn $samples{$sn}{abbr}";
	}
	$cmd .= " > $outdir/Uniq/Uniq.cfg.ini";
	$queue->addcommond("$cmd");
	$queue->addcommond("cd $outdir/Uniq");
	$queue->addcommond("$script/uniq/uniqueFastas.pl -i Uniq.cfg.ini -u unique.fasta -table table.xls -uniq_all -a");
	$queue->run();
	$queue->wait();
}

sub Align_to_rfam{
	$queue->set_job_cpu(40);
	$queue->set_job_mem(40);
	$queue->set_job_name("比对Rfam数据库");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Rfam");
	my $cmd = "cat";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn/$sn.filtered.trimmed.len.uniq.fa";
	}
	$cmd .= " > $outdir/Rfam/merge.unique.fasta";
	$queue->addcommond("$cmd");
	$queue->addcommond("cd $outdir/Rfam");
	#$queue->addcommond("split -l 200000 -a 3 merge.unique.fasta X_");
	$queue->addcommond("blastall -d Rfam -i merge.unique.fasta -o merge.unique.blast.out -p blastn -v 5 -b 5 -F F -e 0.01 -a 40");
	$queue->addcommond("cp $rfam_seed ./");
	my $rfam_seed_basename = basename $rfam_seed;
	$queue->addcommond("$script/rfam/rfam_blast_parse_merge.pl -i merge.unique.blast.out -s $rfam_seed_basename -o blast_table.xls -sumary -merg $outdir/Rfam/merge.unique.fasta -config $outdir/Uniq/Uniq.cfg.ini -unmatch rfam_umatch.list -mirna rfam_miRNA.list -pie pie_out");
	$queue->addcommond("$script/rfam/chooseseq.pl $outdir/Rfam/merge.unique.fasta rfam_miRNA.list > rfam_miRNA.fa");
	$queue->addcommond("$script/rfam/chooseseq.pl $outdir/Rfam/merge.unique.fasta rfam_umatch.list > rfam_unmatched.fa");
	$queue->addcommond("cat rfam_miRNA.fa rfam_unmatched.fa > rfam_trimed.fa");
	$queue->run();
	$queue->wait();
}

sub Bais{
	$queue->set_job_cpu(5);
	$queue->set_job_mem(10);
	$queue->set_job_name("Bais分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Bais");
	$queue->addcommond("cd $outdir/Bais");
	$queue->addcommond("ln -s $outdir/Rfam/merge.unique.fasta ./");
	$queue->addcommond("$script/bias/bias.pl -cfg $outdir/Uniq/Uniq.cfg.ini -fa merge.unique.fasta -min 18 -max 32");
	$queue->run();
	$queue->wait();
}

sub Mapping2genome{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("Mapping2genome");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Mapping2genome");
	$queue->addcommond("cd $outdir/Mapping2genome");
	$queue->addcommond("$script/public/fasta_clean_x01_and_angle_bracket.pl $refer $outdir/Mapping2genome/ref.fa");
	unless(defined($refer_index)){
		$refer_index = "$outdir/Mapping2genome/ref_index";
		$queue->addcommond("bowtie-build $outdir/Mapping2genome/ref.fa ref_index");
	}
	$queue->addcommond("cp $outdir/Rfam/rfam_trimed.fa ./ ");
	$queue->addcommond("rm -rf reads_vs_genome.all.arf");
	$queue->addcommond("$software/mirdeep/2.0.0.7/mapper.pl rfam_trimed.fa -p $refer_index -t reads_vs_genome.all.arf -c -j -q -v -o 30");
	$queue->addcommond("$script/genome/sRNA_genome_distribution_merge.pl -fa rfam_trimed.fa -config $outdir/Uniq/Uniq.cfg.ini -arf reads_vs_genome.all.arf -bar");
	$queue->addcommond("awk '\$0!~/x[1,2,3,4,5]\\t/' reads_vs_genome.all.arf > reads_vs_genome.arf");
	$queue->run();
	$queue->wait();
	
	$refer = "$outdir/Mapping2genome/ref.fa";
}

sub Expression{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("miRNA鉴定与表达量评估");
	$queue->set_work_dir($workdir);
	if ( $class eq 'animal' ) {
		if ( defined($species) ) {
			$queue->addcommond("mkdir -p $outdir/Expression");
			$queue->addcommond("cd $outdir/Expression");
			$queue->addcommond("grep \'>$species\' $miRbase | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $miRbase mature_$species.list > mature_$species.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.fa > mature_$species.dna.fa");
			$queue->addcommond("grep \">\" $miRbase | grep -v \'>$species\' | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.other.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $miRbase mature_$species.other.list > mature_$species.other.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.other.fa > mature_$species.other.dna.fa");
			$queue->addcommond("grep \'>$species\' $hairpin | awk \'{ print \$1 }\' | sed \'s/>//\' > hairpin_$species.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $hairpin hairpin_$species.list > hairpin_$species.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' hairpin_$species.fa > hairpin_$species.dna.fa");
			#$queue->addcommond("$software/mirdeep/2.0.0.7/miRDeep2.pl $outdir/Rfam/rfam_trimed.fa $refer $outdir/Mapping2genome/reads_vs_genome.arf mature_$species.dna.fa mature_$species.other.dna.fa hairpin_$species.dna.fa -c -d 2> report.log");
			$queue->addcommond("$software/mirdeep2/0.1.0/bin/miRDeep2.pl $outdir/Rfam/rfam_trimed.fa $refer $outdir/Mapping2genome/reads_vs_genome.arf mature_$species.dna.fa mature_$species.other.dna.fa hairpin_$species.dna.fa 2> report.log");
			$queue->addcommond("id=\$\($script/public/get_miRDeep_id.pl ./\)");
			$queue->addcommond("$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/Uniq/Uniq.cfg.ini");
			$queue->addcommond("$script/exp_dir/get_mirdeep_count.pl mirdeep/miRNAs_expressed_all_samples.xls > known_miR_count.xls");
			$queue->addcommond("sed 1d mirdeep/predicted_total_sumary_table.xls |sed '\$d' |sed '\$d' |sed 's/#provisional_id/miRNA_ID/' |cut -f 1,3- > novel_miR_count.xls");
			$queue->addcommond("cp mirdeep/miRNAs_expressed_all_samples_normalization.xls known_miR_norm.xls");
			$queue->addcommond("sed 's/#provisional_id/miRNA_ID/' mirdeep/novo_mature_normalization.xls > novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args known_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i known_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre known_miRNA");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i novel_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre novel_miRNA");
		}elsif($animal) {
			$queue->addcommond("mkdir -p $outdir/Expression");
			$queue->addcommond("cd $outdir/Expression");
			$queue->addcommond("cp $database/miRBase/Animal/* ./");
			#$queue->addcommond("$software/mirdeep/2.0.0.7/miRDeep2.pl outdir/Rfam/rfam_trimed.fa $refer $outdir/Mapping2genome/reads_vs_genome.arf animal.mature.dna.fa animal.mature.other.dna.fa animal.hairpin.dna.fa 2> report.log");
			$queue->addcommond("$software/mirdeep2/0.1.0/bin/miRDeep2.pl outdir/Rfam/rfam_trimed.fa $refer $outdir/Mapping2genome/reads_vs_genome.arf animal.mature.dna.fa animal.mature.other.dna.fa animal.hairpin.dna.fa 2> report.log");
			$queue->addcommond("id=\$\($script/public/get_miRDeep_id.pl ./\)");
			$queue->addcommond("$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/Uniq/Uniq.cfg.ini");
			$queue->addcommond("$script/exp_dir/get_mirdeep_count.pl mirdeep/miRNAs_expressed_all_samples.xls > known_miR_count.xls");
			$queue->addcommond("sed 1d mirdeep/predicted_total_sumary_table.xls |sed '\$d' |sed '\$d' |sed 's/#provisional_id/miRNA_ID/' |cut -f 1,3- > novel_miR_count.xls");
			$queue->addcommond("cp mirdeep/miRNAs_expressed_all_samples_normalization.xls known_miR_norm.xls");
			$queue->addcommond("sed 's/#provisional_id/miRNA_ID/' mirdeep/novo_mature_normalization.xls > novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args known_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i known_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre known_miRNA");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i novel_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre novel_miRNA");
		}
	}elsif($class eq 'plant'){
		if ( defined($species) ) {
			$queue->addcommond("mkdir -p $outdir/Expression");
			$queue->addcommond("cd $outdir/Expression");
			$queue->addcommond("grep \'>$species\' $miRbase | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $miRbase mature_$species.list > mature_$species.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.fa > mature_$species.dna.fa");
			$queue->addcommond("grep \">\" $miRbase | grep -v \'>$species\' | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.other.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $miRbase mature_$species.other.list > mature_$species.other.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.other.fa > mature_$species.other.dna.fa");
			$queue->addcommond("grep \'>$species\' $hairpin | awk \'{ print \$1 }\' | sed \'s/>//\' > hairpin_$species.list");
			$queue->addcommond("$script/rfam/chooseseq.pl $hairpin hairpin_$species.list > hairpin_$species.fa");
			$queue->addcommond("sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' hairpin_$species.fa > hairpin_$species.dna.fa");
			$queue->addcommond("ln -sf $refer ./");
			$queue->addcommond("ln -sf $outdir/Rfam/rfam_trimed.fa ./");
			$queue->addcommond("$software/mirdeep2/0.1.0/bin/quantifier.pl -p hairpin_$species.dna.fa -m mature_$species.dna.fa -r rfam_trimed.fa -y dir -k");
			$queue->addcommond("$software/mirdeep2/0.1.0/bin/parse_mappings.pl $outdir/Mapping2genome/reads_vs_genome.arf -a 0 -b 18 -c 25 -i 5 > reads_vs_genome.arf");
			# $queue->addcommond("$software/mirdeep/2.0.0.7/quantifier.pl -p hairpin_$species.dna.fa -m mature_$species.dna.fa -r rfam_trimed.fa -y dir -k");
			# $queue->addcommond("$software/mirdeep/2.0.0.7/parse_mappings.pl $outdir/Mapping2genome/reads_vs_genome.arf -a 0 -b 18 -c 25 -i 5 > reads_vs_genome.arf");
			$queue->addcommond("$script/exp_dir/parse_arf.pl -rfam rfam_trimed.fa -mrd $outdir/Expression/expression_analyses/expression_analyses_dir/miRBase.mrd -arf reads_vs_genome.arf");
			$queue->addcommond("$software/mireap/0.2/bin/mireap.pl -i filtered.fa -m map.txt -r $refer -t Nov");
			$queue->addcommond("$script/exp_dir/parse_mireap_result.pl -exp miRNAs_expressed_all_samples_dir.csv -aln mireap-Nov.aln -gff mireap-Nov.gff -config $outdir/Uniq/Uniq.cfg.ini -o mireap");
			$queue->addcommond("mkdir -p mireap/pdfs");
			$queue->addcommond("cd mireap/pdfs");
			$queue->addcommond("$software/ViennaRNA/2.1.6h/RNAfold < $outdir/Expression/mireap/novel_miR_pre.fa > NovmiR_pre.str");
			$queue->addcommond("for i in *ps; do ps2pdf \$i \$i.pdf; done");
			$queue->addcommond("cd -");
			$queue->addcommond("cp mireap/novel_miR_count.xls ./");
			$queue->addcommond("cp mireap/novel_miR_norm.xls ./");
			$queue->addcommond("$script/exp_dir/get_mireap_exp.pl -exp $outdir/Expression/miRNAs_expressed_all_samples_dir.csv -mrd $outdir/Expression/expression_analyses/expression_analyses_dir/miRBase.mrd -config $outdir/Uniq/Uniq.cfg.ini");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args known_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i known_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre known_miRNA");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i novel_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre novel_miRNA");
		}elsif($plant) {
			$queue->addcommond("mkdir -p $outdir/Expression");
			$queue->addcommond("cd $outdir/Expression");
			$queue->addcommond("ln -sf $refer ./");
			$queue->addcommond("ln -sf $outdir/Rfam/rfam_trimed.fa ./");
			$queue->addcommond("cp $plant_db/* ./");
			$queue->addcommond("$software/mireap/0.2/bin/quantifier.pl -p plant.hairpin.dna.fa -m plant.mature.dna.fa -r rfam_trimed.fa -y dir -k");
			$queue->addcommond("$software/mireap/0.2/bin/parse_mappings.pl $outdir/Mapping2genome/reads_vs_genome.arf -a 0 -b 18 -c 25 -i 5 > reads_vs_genome.arf");
			$queue->addcommond("$script/exp_dir/parse_arf.pl -rfam rfam_trimed.fa -mrd $outdir/Expression/expression_analyses/expression_analyses_dir/miRBase.mrd -arf reads_vs_genome.arf");
			$queue->addcommond("$software/mireap/0.2/bin/mireap.pl -i filtered.fa -m map.txt -r $refer -t Nov");
			$queue->addcommond("$script/exp_dir/parse_mireap_result.pl -exp miRNAs_expressed_all_samples_dir.csv -aln mireap-Nov.aln -gff mireap-Nov.gff -config $outdir/Uniq/Uniq.cfg.ini -o mireap");
			$queue->addcommond("mkdir -p mireap/pdfs");
			$queue->addcommond("cd mireap/pdfs");
			$queue->addcommond("$software/ViennaRNA/2.1.6h/RNAfold < $outdir/Expression/mireap/novel_miR_pre.fa > NovmiR_pre.str");
			$queue->addcommond("for i in *ps; do ps2pdf \$i \$i.pdf; done");
			$queue->addcommond("cd -");
			$queue->addcommond("cp mireap/novel_miR_count.xls ./");
			$queue->addcommond("cp mireap/novel_miR_norm.xls ./");
			$queue->addcommond("$script/exp_dir/get_mireap_exp.pl -exp $outdir/Expression/miRNAs_expressed_all_samples_dir.csv -mrd $outdir/Expression/expression_analyses/expression_analyses_dir/miRBase.mrd -config $outdir/Uniq/Uniq.cfg.ini");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args novel_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/top_ten_miRNA.r --args known_miR_norm.xls");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i known_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre known_miRNA");
			$queue->addcommond("Rscript $script/exp_dir/plot_pca_exp.pl -i novel_miR_norm.xls -o ./ -m $outdir/group.tmp -g group -pre novel_miRNA");
		}
    }
	
	$queue->run();
	$queue->wait();

}

sub Venn{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("miRNA venn 分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Venn");
	
	$queue->addcommond("cd $outdir/Venn");
	$queue->addcommond("ln -sf $outdir/Expression/known_miR_count.xls");
	$queue->addcommond("ln -sf $outdir/Expression/novel_miR_count.xls");
	$queue->addcommond("$script/common/common.pl -count known_miR_count.xls -o known -min 5");
	$queue->addcommond("$script/common/common.pl -count novel_miR_count.xls -o novel -min 5");
	$queue->run();
	$queue->wait();
}

sub Family{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("miRNA family分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Family");
	
	if ( $plant || $animal ) {
		$queue->addcommond("echo 'Family analysis is jumped (under development)!'");
	}else{
		if($class eq 'animal'){
			$queue->addcommond("cd $outdir/Family");
			$queue->addcommond("cp $outdir/Expression/miRNAs_expressed_all_samples*.csv ./miRNAs_expressed.xls");
			$queue->addcommond("$script/family/known_miR_family.pl miRNAs_expressed.xls");
			$queue->addcommond("$script/family/novel_miR_family.pl $outdir/Expression/mature_$species.dna.fa known_miR_family.xls $outdir/Expression/mirdeep/novo_mature_seq.fa > novel_miR_family.xls");
		}elsif($class eq "plant"){
			$queue->addcommond("cd $outdir/Family");
			$queue->addcommond("cp $outdir/Expression/miRNAs_expressed_all_samples*.csv ./miRNAs_expressed.xls");
			$queue->addcommond("$script/family/known_miR_family.pl miRNAs_expressed.xls");
			$queue->addcommond("$script/family/novel_miR_family.pl $outdir/Expression/mature_$species.dna.fa known_miR_family.xls $outdir/Expression/mireap/novel_miR_mature.fa > novel_miR_family.xls");
		}
	}
	$queue->run();
	$queue->wait();
}

sub DiffExpression{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("差异表达量计算");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/Diffexpression");
	$queue->addcommond("mkdir -p $outdir/Diffexpression/known");
	$queue->addcommond("mkdir -p $outdir/Diffexpression/novel");
	$queue->addcommond("cd $outdir/Diffexpression");
	$queue->addcommond("cp $outdir/Expression/known_miR_*.xls known/");
	$queue->addcommond("cp $outdir/Expression/novel_miR_*.xls novel/");
	
	if ( defined($group_file) ) {
		$queue->addcommond("cd $outdir/Diffexpression/known");
		$queue->addcommond("$script/diff/run_DEGseq_for_rep_DE.pl -group $group_file -tpm known_miR_norm.xls -count known_miR_count.xls -qvalue 0.05 -bypvalue");
		$queue->addcommond("cd $outdir/Diffexpression/novel");
		$queue->addcommond("$script/diff/run_DEGseq_for_rep_DE.pl -group $group_file -tpm novel_miR_norm.xls -count novel_miR_count.xls -qvalue 0.05 -bypvalue");
	}else{
		$queue->addcommond("cd $outdir/Diffexpression/known");
		$queue->addcommond("$script/diff/run_DEGseq.pl -tpm known_miR_norm.xls -count known_miR_count.xls -qvalue 0.05 -bypvalue");
		$queue->addcommond("cd $outdir/Diffexpression/novel");
		$queue->addcommond("$script/diff/run_DEGseq.pl -tpm novel_miR_norm.xls -count novel_miR_count.xls -qvalue 0.05 -bypvalue");
	}
	
	
	$queue->run();
	$queue->wait();
}

sub DiffCluster{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("差异表达miRNA的聚类分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/DiffCluster");
	$queue->addcommond("mkdir -p $outdir/DiffCluster/known");
	$queue->addcommond("mkdir -p $outdir/DiffCluster/novel");
	$queue->addcommond("cd $outdir/DiffCluster");
	$queue->addcommond("cp $outdir/DiffCluster/known_miR_*.xls known/");
	$queue->addcommond("cp $outdir/DiffCluster/novel_miR_*.xls novel/");
	
	my $known_DEG_dir;
	my $novel_DEG_dir;

	if ( defined($group_file) ) {
		$known_DEG_dir = "$outdir/Diffexpression/known/replicates_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/replicates_exp";
	}else {
		$known_DEG_dir = "$outdir/Diffexpression/known/samples_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/samples_exp";
	}
	
	$queue->addcommond("cd $outdir/DiffCluster/known");
	$queue->addcommond("cat $known_DEG_dir/*DE.list | sort -u > all_known_DE.list");
	$queue->addcommond("$script/public/tabletools_select.pl -i all_known_DE.list -t $known_DEG_dir/../known_miR_norm.xls -n 1 -head T > all_known_DE.matrix");
	$queue->addcommond("Rscript $script/cluster/heatmap.r all_known_DE.matrix known_cluster");
	
	$queue->addcommond("cd $outdir/DiffCluster/novel");
	$queue->addcommond("cat $novel_DEG_dir/*DE.list | sort -u > all_novel_DE.list");
	$queue->addcommond("$script/public/tabletools_select.pl -i all_novel_DE.list -t $novel_DEG_dir/../novel_miR_norm.xls -n 1 -head T > all_novel_DE.matrix");
	$queue->addcommond("Rscript $script/cluster/heatmap.r all_novel_DE.matrix novel_cluster");
	
	$queue->run();
	$queue->wait();

}

sub DiffTarget{
	return unless ($target);
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("差异表达miRNA的靶基因鉴定分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/DiffTarget");

	$queue->addcommond("$script/public/fasta_trimdup.pl $target > $outdir/DiffTarget/target.fa");

	if ( $class eq 'animal' ) {
		my $known_miRNA_fa;
		if ( defined($species) ) {
			$known_miRNA_fa = "$outdir/Expression/mature_$species.dna.fa";
		}elsif ($animal) {
			$known_miRNA_fa = "$outdir/Expression/animal.mature.dna.fa";
		}
		my $novel_miRNA_fa = "$outdir/Expression/mirdeep/novo_mature_seq.fa";
		$queue->addcommond("cd $outdir/DiffTarget");
		$queue->addcommond("$script/public/chooseseq.pl $known_miRNA_fa $outdir/DiffCluster/known/all_known_DE.list > all_known_DEM.fa");
		$queue->addcommond("$script/public/chooseseq.pl $novel_miRNA_fa $outdir/DiffCluster/novel/all_novel_DE.list > all_novel_DEM.fa");

		$queue->addcommond("$software/miRanda/3.3a/bin/miranda all_known_DEM.fa $outdir/DiffTarget/target.fa -out known_targets_out -sc 150 -en -10 -quiet");
		$queue->addcommond("$software/miRanda/3.3a/bin/miranda all_novel_DEM.fa $outdir/DiffTarget/target.fa -out novel_targets_out -sc 150 -en -10 -quiet");
		$queue->addcommond("$script/target/mi_result.pl -i known_targets_out -o known_targets.xls");
		$queue->addcommond("$script/target/mi_result.pl -i novel_targets_out -o novel_targets.xls");
		if ($anno) {
			$queue->addcommond("$script/public/tabletools_add.pl -i $anno -t known_targets.xls -n 2 -headi T -headt T > known_targets_anno.xls");
			$queue->addcommond("$script/public/tabletools_add.pl -i $anno -t novel_targets.xls -n 2 -headi T -headt T > novel_targets_anno.xls");
		}else{
			$queue->addcommond("cp known_targets.xls known_targets_anno.xls");
			$queue->addcommond("cp novel_targets.xls novel_targets_anno.xls");
		}
	}elsif( $class eq 'plant' ){
		my $known_miRNA_fa;
		if ( defined($species) ) {
			$known_miRNA_fa = "$outdir/Expression/mature_$species.dna.fa";
		}
		elsif ($plant) {
			$known_miRNA_fa = "$outdir/Expression/plant.mature.dna.fa";
		}
		my $novel_miRNA_fa = "$outdir/Expression/mireap/novel_miR_mature.fa";
		$queue->addcommond("cd $outdir/Target");
		$queue->addcommond("$script/public/chooseseq.pl $known_miRNA_fa $outdir/DiffCluster/known/all_known_DE.list > all_known_DEM.fa");
		$queue->addcommond("$script/public/chooseseq.pl $novel_miRNA_fa $outdir/DiffCluster/novel/all_novel_DE.list > all_novel_DEM.fa");

		$queue->addcommond("$software/psRobot/1.2/psRobot_tar -s all_known_DEM.fa -t $outdir/DiffTarget/target.fa -o known_targets_out -ts 5 -p 15 -gn 1");
		$queue->addcommond("$software/psRobot/1.2/psRobot_tar -s all_novel_DEM.fa -t $outdir/DiffTarget/target.fa -o novel_targets_out -ts 5 -p 15 -gn 1");
		$queue->addcommond("$script/target/tr_psRobot_result2table.pl known_targets_out known_targets.xls");
		$queue->addcommond("$script/target/tr_psRobot_result2table.pl novel_targets_out novel_targets.xls");
		if ($anno) {
			$queue->addcommond("$script/public/tabletools_add.pl -i $anno -t known_targets.xls -n 2 -headi T -headt T > known_targets_anno.xls");
			$queue->addcommond("$script/public/tabletools_add.pl -i $anno -t novel_targets.xls -n 2 -headi T -headt T > novel_targets_anno.xls");
		}else{
			$queue->addcommond("cp known_targets.xls known_targets_anno.xls");
			$queue->addcommond("cp novel_targets.xls novel_targets_anno.xls");
		}
	}
	$queue->run();
	$queue->wait();
}

sub DiffTargetEnrichment{
	return unless ($target);
	return unless ($anno);
	$queue->set_job_cpu(10);
	$queue->set_job_mem(20);
	$queue->set_job_name("差异表达miRNA的靶基因功能富集分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/DiffTargetEnrichment");
	
	my $known_DEG_dir;
	my $novel_DEG_dir;
	
	if ( defined($group_file) ) {
		$known_DEG_dir = "$outdir/Diffexpression/known/replicates_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/replicates_exp";
	}else {
		$known_DEG_dir = "$outdir/Diffexpression/known/samples_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/samples_exp";
	}
	
	$queue->addcommond("cd $outdir/DiffTargetEnrichment");
	$queue->addcommond("cut -f 1,5 $outdir/DiffTarget/known_targets_anno.xls | sed 1d | sort -u > all_known_degtar.xls");
	$queue->addcommond("cut -f 1,5 $outdir/DiffTarget/novel_targets_anno.xls | sed 1d | sort -u > all_novel_degtar.xls");
	$queue->addcommond("for i in \$(find $known_DEG_dir -name \"*DE.list\"); do $script/public/tabletools_select.pl -i \$i -t all_known_degtar.xls -n 1 -head F | cut -f 2 | sort -u > \$i.known.tmp; done");
	$queue->addcommond("for i in \$(find $novel_DEG_dir -name \"*DE.list\"); do $script/public/tabletools_select.pl -i \$i -t all_novel_degtar.xls -n 1 -head F | cut -f 2 | sort -u > \$i.novel.tmp; done");
	$queue->addcommond("mv $known_DEG_dir/*known.tmp ./");
	$queue->addcommond("mv $novel_DEG_dir/*novel.tmp ./");
	$queue->addcommond("rename 'DE.list.known.tmp' 'KnownTar_DE.list' *");
	$queue->addcommond("rename 'DE.list.novel.tmp' 'NovelTar_DE.list' *");
	
	if($kegg){
		$queue->addcommond("for i in *_DE.list; do kegg_enrichment.pl -de \$i -all $gene -asso $kegg --pre \$(basename \$i) -pval 0.05 -org $species; done");
	}
	
	if($go){
		$queue->addcommond("for i in *_DE.list; do go_enrichment.pl -de \$i -all $gene -asso $go --pre \$(basename \$i) -pval 0.05; done");
	}

	$queue->run();
	$queue->wait();
	
}

sub result_collection{
	$queue->set_job_cpu(5);
	$queue->set_job_mem(5);
	$queue->set_job_name("分析结果整理");
	$queue->set_work_dir($workdir);
	
	my ($num, $id, $dir);
	$num = 0;
	
	$queue->addcommond("mkdir -p $outdir/Result");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.RawdataQC";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/fastq_QC/rawdatashow.xls $dir/rawdata_stat.xls");
	foreach my $sn ( sort keys %samples ) {
		$queue->addcommond("cp $outdir/fastq_QC/$sn/$sn.stat.qual.pdf $dir/$sn.qual.pdf");
    }
	$queue->addcommond("cp $outdir/fastq_QC/QC_stat.xls $dir/qc_stat.xls");
	foreach my $sn ( sort keys %samples ) {
		$queue->addcommond("cp $outdir/fastq_QC/$sn/$sn.length_distribution.png $dir/$sn.length_distribution.png");
		$queue->addcommond("cp $outdir/fastq_QC/$sn/$sn.length_distribution.xls $dir/$sn.length_distribution.xls");
    }
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Unique";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/Uniq/*.pdf $dir");
	$queue->addcommond("cp $outdir/Uniq/sumary.xls $dir");

	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Bias";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/Bais/*per.xls $dir");
	$queue->addcommond("cp $outdir/Bais/*.pdf $dir");

	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Rfam";
	$queue->addcommond("mkdir -p $dir");
	foreach my $sn ( sort keys %samples ) {
		$queue->addcommond("cp $outdir/Rfam/$sn\_rfam_stat.xls $dir/$sn.rfam_stat.xls");
		$queue->addcommond("cp $outdir/Rfam/pie_out/$sn\_rfam.pdf $dir/$sn.rfam_stat.pdf");
    }
	$queue->addcommond("cp $outdir/Rfam/All_rfam_stat.xls $dir/all.rfam_stat.xls");
	$queue->addcommond("cp $outdir/Rfam/pie_out/all_rfam.pdf $dir/all.rfam_stat.pdf");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Mapping2genome";
	$queue->addcommond("mkdir -p $dir");
	foreach my $sn ( sort keys %samples ) {
		$queue->addcommond("cp $outdir/Mapping2genome/$sn\_map.pdf $dir/$sn.map_stat.pdf");
		$queue->addcommond("cp $outdir/Mapping2genome/$sn\_map_stat.xls $dir/$sn.map_stat.xls");
    }
	$queue->addcommond("cp $outdir/Mapping2genome/All_map_stat.xls $dir/all.map_stat.xls");
	$queue->addcommond("cp $outdir/Mapping2genome/All_map.pdf $dir/all.map_stat.pdf");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Expression";
	$queue->addcommond("mkdir -p $dir/known $dir/novel");
	$queue->addcommond("mkdir -p $dir/known/2d_structure");
	$queue->addcommond("mkdir -p $dir/novel/2d_structure");
	
	$queue->addcommond("cp $outdir/Expression/known_miR_norm.xls  $dir/known/known_miRNA.norm.xls");
	$queue->addcommond("cp $outdir/Expression/known_miR_count.xls $dir/known/known_miRNA.count.xls");
	$queue->addcommond("cp $outdir/Expression/known_miR_norm.xls.top10.pdf $dir/known/known_miRNA.top10.pdf");
	$queue->addcommond("cp -r $outdir/Expression/pdf*/*.pdf $dir/known/2d_structure");
	if ( $class eq 'animal' ) {
		$queue->addcommond("cp $outdir/Expression/mirdeep/miRBase.mrd $dir/known/known_miRNA.mrd");
		if ( defined($species) ) {
			$queue->addcommond("cp $outdir/Expression/mature_$species.dna.fa  $dir/known/known_mature.fa");
			$queue->addcommond("cp $outdir/Expression/hairpin_$species.dna.fa $dir/known/known_hairpin.fa");
		}elsif($animal) {
			$queue->addcommond("cp $outdir/Expression/animal.mature.dna.fa  $dir/known/known_mature.fa");
			$queue->addcommond("cp $outdir/Expression/animal.hairpin.dna.fa $dir/known/known_hairpin.fa");
		}
	}elsif($class eq 'plant'){
		$queue->addcommond("cp $outdir/Expression/expression_analyses/expression_analyses_dir/miRBase.mrd $dir/known/known_miRNA.mrd");
		if ( defined($species) ) {
			$queue->addcommond("cp $outdir/Expression/mature_$species.dna.fa $dir/known/known_mature.fa");
			$queue->addcommond("cp $outdir/Expression/hairpin_$species.dna.fa $dir/known/known_hairpin.fa");
		}elsif($plant) {
			$queue->addcommond("cp $outdir/Expression/plant.mature.dna.fa $dir/known/known_mature.fa");
			$queue->addcommond("cp $outdir/Expression/plant.hairpin.dna.fa $dir/known/known_hairpin.fa");
		}
	}
	$queue->addcommond("cp $outdir/Expression/novel_miR_norm.xls  $dir/novel/novel_miRNA.norm.xls");
	$queue->addcommond("cp $outdir/Expression/novel_miR_count.xls $dir/novel/novel_miRNA.count.xls");
	$queue->addcommond("cp $outdir/Expression/novel_miR_norm.xls.top10.pdf $dir/novel/novel_miRNA.top10.pdf");
	if ( $class eq 'animal' ) {
		$queue->addcommond("mv $dir/known/2d_structure/[^$species]* $dir/novel/2d_structure");
		$queue->addcommond("cp $outdir/Expression/mirdeep/novo_mature_seq.fa $dir/novel/novel_mature.fa");
		$queue->addcommond("cp $outdir/Expression/mirdeep/novo_precursor_seq.fa $dir/novel/novel_hairpin.fa");
		$queue->addcommond("cp $outdir/Expression/mirdeep/predict.mrd $dir/novel/novel_miRNA.mrd");
	}elsif($class eq 'plant'){
		$queue->addcommond("cp $outdir/Expression/mireap/pdfs/*.pdf $dir/novel_miRNA_structure");
		$queue->addcommond("cp $outdir/Expression/mireap/novel_miR_pre.fa $dir/novel_hairpin.fa");
		$queue->addcommond("cp $outdir/Expression/mireap/novel_miR_mature.fa $dir/novel_mature.fa");
		$queue->addcommond("cp $outdir/Expression/mireap-Nov.aln $dir/novel_miR.aln");
	}
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Family";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/Family/family.species.xls $dir");
	$queue->addcommond("cp $outdir/Family/known_miR_family.xls $dir/known_miRNA_family.xls");
	$queue->addcommond("cp $outdir/Family/novel_miR_family.xls $dir/novel_miRNA_family.xls");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.Venn";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp $outdir/Venn/*_vs_* $dir");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.DiffExpresson";
	$queue->addcommond("mkdir -p $dir/known $dir/novel");
	my $known_DEG_dir;
	my $novel_DEG_dir;
	if ( defined($group_file) ) {
		$known_DEG_dir = "$outdir/Diffexpression/known/replicates_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/replicates_exp";
	}else {
		$known_DEG_dir = "$outdir/Diffexpression/known/samples_exp";
		$novel_DEG_dir = "$outdir/Diffexpression/novel/samples_exp";
	}
	$queue->addcommond("cp $known_DEG_dir/*pdf $known_DEG_dir/*xls $known_DEG_dir/*list $dir/known");
	$queue->addcommond("cp $novel_DEG_dir/*pdf $novel_DEG_dir/*xls $novel_DEG_dir/*list $dir/novel");
	
	$num = $num + 1;
	$id  = &convert_num_to_id($num);
	$dir = "$outdir/Result/$id.DiffCluster";
	$queue->addcommond("mkdir -p $dir");
	$queue->addcommond("cp -r $outdir/DiffCluster/known/known_cluster $dir/known");
	$queue->addcommond("cp -r $outdir/DiffCluster/novel/novel_cluster $dir/novel");
	
	if($target){
		$num = $num + 1;
		$id  = &convert_num_to_id($num);
		$dir = "$outdir/Result/$id.DiffTarget";
		$queue->addcommond("mkdir -p $dir");
		$queue->addcommond("cp $outdir/DiffTarget/*targets_anno.xls $dir");
	}
	
	if($target && $anno){
		$num = $num + 1;
		$id  = &convert_num_to_id($num);
		$dir = "$outdir/Result/$id.DiffTargetEnrichment";
		$queue->addcommond("mkdir -p $dir");
		if($go){
			$queue->addcommond("cp $outdir/DiffTargetEnrichment/*go_enrichment.xls $dir");
			$queue->addcommond("cp $outdir/DiffTargetEnrichment/*go_enrichment*.pdf $dir");
			$queue->addcommond("cp $outdir/DiffTargetEnrichment/*go_enrichment*.png $dir");
		}
		if($kegg){
			$queue->addcommond("cp $outdir/DiffTargetEnrichment/*kegg_enrichment.xls $dir");
			$queue->addcommond("cp $outdir/DiffTargetEnrichment/*kegg_enrichment*.pdf $dir");
		}
	}
	
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



sub list_species {
    open FIN,"$Bin/database/miRBase/organisms.txt";
	<FIN>;
	while(<FIN>){
		my ($a,$b,$c) = split /\t/,$_;
		print "$a\t$c\n";
	}
	exit;
}

sub usage{
	my $program = basename($0);
    die "
Usage   :
	$program [options] <input.fq.lst> <out_dir>
Author  :
	guantao.zheng
Contact :
    guantao.zheng\@origin-gene.com
Update  : 2018.10.22 1st version complete.
Version : 
    1.0
Discription:
	This program is used for miRNA analysis.
	<input.fq.lst> are basic inputs, in default, results will be placed in <out_dir>.
Options: 
  Input Options:
    -input        FILE    input file, contains three cols: sample_name, fq_file, sample_abbr. required
    -refer        FILE    Reference fasta file. required
    -group        FILE    group file, contains two cols: sample_name, group_name.
    -library      FILE    library type: true-seq, VAHTS or exosome is OK. default is true-seq

    -target       FILE    Target fasta file for target prediction, 3UTR, EST, cDNA.
    -anno         FILE    Annotation file from Ensemble
                          1st col should be same with target ID !!!
                          2nd is the gene name which are sampel with go and kegg file !!!
	
    -go           FILE    gene go association file, 1st is gene name, 2nd is go ids, separate by ';'
    -kegg         FILE	  gene kegg association file, 1st is gene name, 2nd is kegg gene ids.
    -gene         FILE	  all genes list of the species

  Output Options:
    -outdir       STRING  output dir, default is out

  Analysis Options:
    -class        STRING  object class: 'animal' or 'plant', default is 'animal'
	
    -species      STRING  Specie name, use miRbase data of this species. Don't set this if species is unknown
    -list         STRING  Use \"$program -list\" to see wheather the supported species in the default database.
    -animal       NULL    Use miRbase data of all animal. Set this if species is unknown, and object is animal.
    -plant        NULL    Use miRbase data of all plant. Set this if species is unknown, and object is plant.
	
    -animal_db    DIR     The directory contains \"animal.mature.dna.fa animal.mature.other.dna.fa animal.hairpin.dna.fa\",
                          default \"$Bin/database/miRBase/Animal\".
    -plant_db     DIR     The directory contains \"plant.mature.dna.fa plant.hairpin.dna.fa\",
                          default \"$Bin/database/miRBase/Plant\"
						  
    -rfam_seed    FILE    Rfam Seed File, default \"$Bin/database/Rfam/Rfam.seed\".
    -miRbase      FILE    miRbase file, default \"$Bin/database/miRBase/mature.fa\".
    -hairpin      FILE    hairpin file, default \"$Bin/database/miRBase/hairpin.fa\".

  Task Options:
    -max_cpus     INI    max cpu number limitation, default is 200 [ optional ]
    -max_jobs     INI    max job number limitation, default is 20  [ optional ]
    -queue_name   STRING queue name, default is yxsw
    -local        NULL   local model
    -rerun        NULL   rerun model

    -help         NULL    Help, display help information to screen.

#################################################################################################
input.fq.lst file format(TAB separate):
sample_1  /path/Sample_HZS_1.fq  ZS1
sample_2  /path/Sample_HZS_2.fq  ZS2
sample_3  /path/Sample_HZS_3.fq  ZS3
#################################################################################################
group file format(TAB separate): header is needed !!
id	group
sample_1	cond_A
sample_2	cond_A
sample_3	cond_A
sample_4	cond_B
sample_5	cond_B
sample_6	cond_B
#################################################################################################
anno file format(TAB separate)
Transcript	Genes	Names	Description
ENST00000456328	ENSG00000223972	DDX11L1	DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1
ENST00000450305	ENSG00000223972	DDX11L1	DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1
ENST00000488147	ENSG00000227232	WASH7P	WAS protein family homolog 7 pseudogene
#################################################################################################
\n"
}

