#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use Cwd qw/getcwd abs_path/;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script $RealBin);

my $script   = "$Bin/script";
my $version  = "v1.0";
my $database = "$Bin/database"; 
my $software = "$Bin/software"; 

my ($help, $file_input, $outdir, $file_meta,$refer, $refer_index, $library, $rfam_seed, $miRbase, $hairpin, $species,$list_species,$plant, $animal, $class,$novel );
my %samples;

GetOptions(
	"help!"         => \$help,
	"input=s"       => \$file_input,
	"meta:s"        => \$file_meta,
	"library:s"     => \$library,
	"outdir:s"      => \$outdir,

	"refer:s"       => \$refer,
	"refer_index:s" => \$refer_index,
	"rfam_seed:s"   => \$rfam_seed,
	"miRbase:s"     => \$miRbase,
	"hairpin:s"     => \$hairpin,

	"novel!"        => \$novel,
	"class:s"       => \$class,
	"species:s"     => \$species,
	"list!"         => \$list_species,
	"animal!"       => \$animal,
	"plant!"        => \$plant,
);

sub usage{
Arguments:
my $program = basename($0);
    die "
  [Input Options]:
    --help         <null>      print help information to screen
    --input        <string>    input file, contains three cols: sample_name, fq_file, sample_abbr. required
    --meta         <string>    meta information for grouping comparison. [ optional ]
    --library      <string>    library type: true-seq, VAHTS or exosome is OK. default is VAHTS
    --refer        <string>    reference genome fasta file, [ required ].
    --refer_index  <string>    reference genome bowtie index, [ optional ].
    --outdir       <string>    output dir, default is out
    --novel        <null>      identify novel miRNA or not
    --class        <string>    object class: animal, mammal, plant, vertebrate, invertebrate or virus. default is 'animal'
    --species      <string>    Specie name, use miRbase data of this species. Do not set this if species is unknown
    --list         <string>    Use \"$program -list\" to see wheather the supported species in the default database.
    --animal       <null>      Use miRbase data of all animal. Set this if species is unknown, and object is animal.
    --plant        <null>      Use miRbase data of all plant. Set this if species is unknown, and object is plant.
    --rfam_seed    <string>    Rfam Seed File, default \"$Bin/database/Rfam/Rfam.seed\".
    --miRbase      <string>    miRbase file, default \"$Bin/database/miRBase/mature.fa\".
    --hairpin      <string>    hairpin file, default \"$Bin/database/miRBase/hairpin.fa\".
\n"    
}

$outdir      ||= "out";
$outdir = abs_path( $outdir );
system("mkdir -p $outdir");
my $workdir  = getcwd;
system("mkdir -p $workdir/cmd");
my (@exp_class);

sub parsing_parameters{
	        &usage if ($help);
			$class     ||= "animal";
			$library   ||= "VAHTS";

			@exp_class = qw/known/;
			if($novel){
				@exp_class = qw/known novel/;
			}

			$rfam_seed ||= "$Bin/database/Rfam/Rfam.seed";
			$miRbase   ||= "$Bin/database/miRBase/mature.fa";
			$hairpin   ||= "$Bin/database/miRBase/hairpin.fa";
}

&parsing_parameters();

open FIN, $file_input or die "can not open raw fq list file!";
while (<FIN>){
    chomp ;
    my @s = split /\s+/,$_;
    next if @s != 3;
    if(exists $samples{$s[0]}){
        die "sample_name $s[0] is duplicate! please check input file[$file_input]!";
    }
    $samples{$s[0]}{data} = $s[1];
    $samples{$s[0]}{abbr} = $s[2];
}
close FIN;

&rawdata_quality_control();
&rawdata_qc_merge();
&obtain_uniq_seq();
&align_to_rfam();
&mapping_to_genome();
&structure_bais();
&identification();
&explevel_analysis();

sub rawdata_quality_control {
	my $pm = Parallel::ForkManager->new(10);
	foreach my $sn(sort keys %samples){
	$pm->start and next;
		open FOUT,">$workdir/cmd/run.quality_control_${sn}.sh";
		print FOUT "mkdir -p $outdir/rawdata_qc/$sn\n";
		my @fqs = split /,/,$samples{$sn}{data};
		print FOUT "rm -f $outdir/rawdata_qc/$sn/$sn.fq\n";

		foreach my $fq(@fqs){
			if($fq =~ /.gz$/){
				print FOUT "zcat $fq >> $outdir/rawdata_qc/$sn/$sn.fq\n";
			}else{
				print FOUT "cat  $fq >> $outdir/rawdata_qc/$sn/$sn.fq\n";
			}
		}

		print FOUT "$script/qc/trim_fastq.pl -input $outdir/rawdata_qc/$sn/$sn.fq -out $outdir/rawdata_qc/$sn/$sn.cut.fq -type cut -start 1 -len 51\n";
		print FOUT "$software/fastx_toolkit/0.0.14/bin/fastx_quality_stats -i $outdir/rawdata_qc/$sn/$sn.cut.fq -o $outdir/rawdata_qc/$sn/$sn.stat -Q 33\n";
		print FOUT "$software/fastx_toolkit/0.0.14/bin/fastq_quality_boxplot_graph.sh -i $outdir/rawdata_qc/$sn/$sn.stat -o $outdir/rawdata_qc/$sn/$sn.stat.qual.ps -p\n";
		print FOUT "ps2pdf $outdir/rawdata_qc/$sn/$sn.stat.qual.ps $outdir/rawdata_qc/$sn/$sn.stat.qual.pdf\n";
		print FOUT "echo -e '$sn\\t$outdir/rawdata_qc/$sn/$sn.fq\\n' > $outdir/rawdata_qc/$sn/$sn.rawfq.list\n";
		print FOUT "java -jar $script/qa/FastqStat.jar -i $outdir/rawdata_qc/$sn/$sn.rawfq.list > $outdir/rawdata_qc/$sn/$sn.rawdata.stat\n";
		print FOUT "cut -f 1-3,11-14 $outdir/rawdata_qc/$sn/$sn.rawdata.stat > $outdir/rawdata_qc/$sn/$sn.rawdatashow.xls\n";
	
		if($library eq "exosome"){
			print FOUT "$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $outdir/rawdata_qc/$sn/$sn.cut.fq -o $outdir/rawdata_qc/$sn/$sn.filtered.fq -a AAAAAAAAAA -Q 33 -l 18 -v > $outdir/rawdata_qc/$sn/$sn.qc.txt\n";
		}elsif($library eq "VAHTS"){
			print FOUT "$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $outdir/rawdata_qc/$sn/$sn.cut.fq -o $outdir/rawdata_qc/$sn/$sn.filtered.fq -a AGATCGGAAGAGCAC -Q 33 -l 18 -v > $outdir/rawdata_qc/$sn/$sn.qc.txt\n";
		}elsif($library eq "true-seq"){
			print FOUT "$software/fastx_toolkit/0.0.14/bin/fastx_clipper -i $outdir/rawdata_qc/$sn/$sn.cut.fq -o $outdir/rawdata_qc/$sn/$sn.filtered.fq -a TGGAATTCTCGGGTGCCAAGG -Q 33 -l 18 -v > $outdir/rawdata_qc/$sn/$sn.qc.txt\n";
		}

		print FOUT "cd $outdir/rawdata_qc/$sn\n";
		print FOUT "$script/qc/DynamicTrim.pl $sn.filtered.fq -h 20 -bwa\n";
		print FOUT "$script/qc/fastq_to_fasta.pl -i $sn.filtered.fq.trimmed -o $sn.filtered.trimmed.fa\n";
		print FOUT "$script/qc/trim_seq.pl -f $sn.filtered.trimmed.fa -o $sn.filtered.trimmed.len.fa -min 18 -max 32 >> $sn.qc.txt\n";
		print FOUT "$script/qc/fasta_clipping_histogram_unique.pl $sn.filtered.trimmed.len.fa $sn.length_distribution.png > $sn.length_distribution.xls\n";
		close FOUT;
		system("sh $workdir/cmd/run.quality_control_${sn}.sh > $workdir/cmd/run.quality_control_${sn}.o 2> $workdir/cmd/run.quality_control_${sn}.e ");
	$pm->finish;
	}
	$pm->wait_all_children;
}

sub rawdata_qc_merge {
	open FOUT,">$workdir/cmd/run.rawdata_qc.merge.sh";
	print FOUT "cd $outdir/rawdata_qc\n";
	print FOUT "rm -rf qc.lst\n";
	
	foreach my $sn(sort keys %samples){
		print FOUT "echo $outdir/rawdata_qc/$sn/$sn.qc.txt >> $outdir/rawdata_qc/qc.lst\n";
	}
	print FOUT "$script/qc/get_miRNA_fq_data_table.pl $outdir/rawdata_qc/qc.lst > $outdir/rawdata_qc/QC_stat.xls\n";

	my $cmd = "$script/public/cat_files.pl -head 1 -outfile $outdir/rawdata_qc/rawdatashow.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.rawdatashow.xls";
	}
	print FOUT $cmd;
	close FOUT;
	system("sh $workdir/cmd/run.rawdata_qc.merge.sh > $workdir/cmd/run.rawdata_qc.merge.o 2> $workdir/cmd/run.rawdata_qc.merge.e ");
}


sub obtain_uniq_seq {
	open FOUT,">$workdir/cmd/run.obtain_uniq_seq.sh";
	print FOUT "mkdir -p $outdir/uniq\n";

	my $cmd = "$script/uniq/get_uniqueFastas_cfg.pl $outdir/uniq";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.filtered.trimmed.len.fa $sn $samples{$sn}{abbr}";
	}
	$cmd .= " > $outdir/uniq/uniq.cfg.ini";
	print FOUT $cmd."\n";

	print FOUT "cd $outdir/uniq\n";
	print FOUT "$script/uniq/uniqueFastas.pl -i uniq.cfg.ini -u unique.fasta -table table.xls -uniq_all -a\n";
	close FOUT;
	system("sh $workdir/cmd/run.obtain_uniq_seq.sh > $workdir/cmd/run.obtain_uniq_seq.o 2> $workdir/cmd/run.obtain_uniq_seq.e ");
}


sub align_to_rfam {
	open FOUT,">$workdir/cmd/run.align_to_rfam.sh";
	print FOUT "module load blast/2.2.26\n";
	print FOUT "mkdir -p $outdir/rfam\n";

	my $cmd = "cat";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/rawdata_qc/$sn/$sn.filtered.trimmed.len.uniq.fa";
	}
	$cmd .= " > $outdir/rfam/merge.unique.fasta";
	print FOUT $cmd."\n";
	
	print FOUT "cd $outdir/rfam\n";
	print FOUT "blastall -d Rfam -i merge.unique.fasta -o merge.unique.blast.out -p blastn -v 5 -b 5 -F F -e 0.01 -a 30\n";

	print FOUT "cp $rfam_seed ./\n";
	my $rfam_seed_basename = basename $rfam_seed;
	print FOUT "$script/rfam/rfam_blast_parse_merge.pl -i merge.unique.blast.out -s $rfam_seed_basename -o blast_table.xls -sumary -merg $outdir/rfam/merge.unique.fasta -config $outdir/uniq/uniq.cfg.ini -unmatch rfam_umatch.list -mirna rfam_miRNA.list -pie pie_out\n";
	print FOUT "$script/rfam/chooseseq.pl $outdir/rfam/merge.unique.fasta rfam_miRNA.list > rfam_miRNA.fa\n";
	print FOUT "$script/rfam/chooseseq.pl $outdir/rfam/merge.unique.fasta rfam_umatch.list > rfam_unmatched.fa\n";
	print FOUT "cat rfam_miRNA.fa rfam_unmatched.fa > rfam_trimed.fa\n";
	close FOUT;
	system("sh $workdir/cmd/run.align_to_rfam.sh > $workdir/cmd/run.align_to_rfam.o 2> $workdir/cmd/run.align_to_rfam.e ");
}

sub mapping_to_genome {
	open FOUT,">$workdir/cmd/run.mapping_to_genome.sh";
	print FOUT "module load bowtie/0.12.9\n";
	print FOUT "module load mirdeep/2.0.0.7\n";	
	print FOUT "mkdir -p $outdir/mapping\n";
	print FOUT "cd $outdir/mapping\n";
	print FOUT "$script/public/fasta_clean_x01_and_angle_bracket.pl $refer $outdir/mapping/ref.fa\n";

	unless(defined($refer_index)){
		$refer_index = "$outdir/mapping/ref_index";
		print FOUT "bowtie-build $outdir/mapping/ref.fa ref_index\n";
	}

	print FOUT "cp $outdir/rfam/rfam_trimed.fa ./ \n";
	print FOUT "rm -rf reads_vs_genome.all.arf\n";
	print FOUT "mapper.pl rfam_trimed.fa -p $refer_index -t reads_vs_genome.all.arf -c -j -q -v -o 30\n";
	print FOUT "$script/genome/sRNA_genome_distribution_merge.pl -fa rfam_trimed.fa -config $outdir/uniq/uniq.cfg.ini -arf reads_vs_genome.all.arf -bar\n";
	print FOUT "awk '\$0!~/x[1,2,3,4,5]\\t/' reads_vs_genome.all.arf > reads_vs_genome.arf\n";
	$refer = "$outdir/mapping/ref.fa";
	close FOUT;  
	system("sh $workdir/cmd/run.mapping_to_genome.sh > $workdir/cmd/run.mapping_to_genome.o 2> $workdir/cmd/run.mapping_to_genome.e ");
}


sub structure_bais {
	open FOUT,">$workdir/cmd/run.structure_bais.sh";
	print FOUT "mkdir -p $outdir/bais\n";
	print FOUT "cd $outdir/bais\n";
	print FOUT "ln -s $outdir/rfam/merge.unique.fasta ./\n";
	print FOUT "$script/bias/bias.pl -cfg $outdir/uniq/uniq.cfg.ini -fa merge.unique.fasta -min 18 -max 32\n";
	close FOUT;
	system("sh $workdir/cmd/run.structure_bais.sh > $workdir/cmd/run.structure_bais.o 2> $workdir/cmd/run.structure_bais.e ");
}


sub identification {
	open FOUT,">$workdir/cmd/run.identification.sh";
	print FOUT "module load bowtie/0.12.9\n";
	print FOUT "module load mirdeep/2.0.0.7\n";
	print FOUT "module load ViennaRNA/2.1.8\n";
	print FOUT "module load randfold/2.0.1\n";
	print FOUT "mkdir -p $outdir/identification\n";
	print FOUT "cd $outdir/identification\n";
	if ( $class eq 'animal' ) {
 		if ( defined($species) ) {
			print FOUT "grep \'>$species\' $miRbase | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.list\n";
			print FOUT "$script/rfam/chooseseq.pl $miRbase mature_$species.list > mature_$species.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.fa > mature_$species.dna.fa\n";
			print FOUT "grep \">\" $miRbase | grep -v \'>$species\' | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.other.list\n";
			print FOUT "$script/rfam/chooseseq.pl $miRbase mature_$species.other.list > mature_$species.other.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.other.fa > mature_$species.other.dna.fa\n";
			print FOUT "grep \'>$species\' $hairpin | awk \'{ print \$1 }\' | sed \'s/>//\' > hairpin_$species.list\n";
			print FOUT "$script/rfam/chooseseq.pl $hairpin hairpin_$species.list > hairpin_$species.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' hairpin_$species.fa > hairpin_$species.dna.fa\n";
			print FOUT "miRDeep2.pl $outdir/rfam/rfam_trimed.fa $refer $outdir/mapping/reads_vs_genome.arf mature_$species.dna.fa mature_$species.other.dna.fa hairpin_$species.dna.fa -c -d 2> report.log\n";
			print FOUT "id=\$\($script/public/get_miRDeep_id.pl ./\)\n";
			print FOUT "$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/uniq/uniq.cfg.ini\n";
		}elsif($animal){
			print FOUT "cp $database/miRBase/Animal/* ./\n";
			print FOUT "miRDeep2.pl $outdir/rfam/rfam_trimed.fa $refer $outdir/mapping/reads_vs_genome.arf animal.mature.dna.fa animal.mature.other.dna.fa animal.hairpin.dna.fa 2> report.log\n";
			print FOUT "id=\$\($script/public/get_miRDeep_id.pl ./\)\n";
			print FOUT "$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/uniq/uniq.cfg.ini\n";
		}

	}elsif($class eq 'plant'){
		if ( defined($species) ) {
			print FOUT "grep \'>$species\' $miRbase | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.list\n";
			print FOUT "$script/rfam/chooseseq.pl $miRbase mature_$species.list > mature_$species.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.fa > mature_$species.dna.fa\n";
			print FOUT "grep \">\" $miRbase | grep -v \'>$species\' | awk \'{ print \$1 }\' | sed \'s/>//\' > mature_$species.other.list\n";
			print FOUT "$script/rfam/chooseseq.pl $miRbase mature_$species.other.list > mature_$species.other.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' mature_$species.other.fa > mature_$species.other.dna.fa\n";
			print FOUT "grep \'>$species\' $hairpin | awk \'{ print \$1 }\' | sed \'s/>//\' > hairpin_$species.list\n";
			print FOUT "$script/rfam/chooseseq.pl $hairpin hairpin_$species.list > hairpin_$species.fa\n";
			print FOUT "sed \'s/^>\\(\\S*\\)\\s*.*\$/>\\1/g;/^[^>]/ s/U/T/g\' hairpin_$species.fa > hairpin_$species.dna.fa\n";
			print FOUT "miRDeep2.pl $outdir/rfam/rfam_trimed.fa $refer $outdir/mapping/reads_vs_genome.arf mature_$species.dna.fa mature_$species.other.dna.fa hairpin_$species.dna.fa -c -d 2> report.log\n";
			print FOUT "id=\$\($script/public/get_miRDeep_id.pl ./\)\n";
			print FOUT "$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/uniq/uniq.cfg.ini\n";
		}elsif($plant) {
			print FOUT "cp $database/miRBase/Plant/* ./\n";
			print FOUT "miRDeep2.pl $outdir/rfam/rfam_trimed.fa $refer $outdir/mapping/reads_vs_genome.arf plant.mature.dna.fa plant.mature.other.dna.fa plant.hairpin.dna.fa 2> report.log\n";
			print FOUT "id=\$\($script/public/get_miRDeep_id.pl ./\)\n";
			print FOUT "$script/exp_dir/convent_mirdeep_result.pl -id \$id -o mirdeep -config $outdir/uniq/uniq.cfg.ini\n";
		}
	}

	close FOUT;
	system("sh $workdir/cmd/run.identification.sh > $workdir/cmd/run.identification.o 2> $workdir/cmd/run.identification.e ");
}

sub explevel_analysis {
	foreach my $xxx(@exp_class){
		open FOUT,">$workdir/cmd/run.explevel_${xxx}.matrix.sh";
		print FOUT "mkdir -p $outdir/explevel_analysis/$xxx\n";
		if($xxx eq "known"){
			print FOUT "$script/exp_dir/get_mirdeep_count.pl $outdir/identification/mirdeep/miRNAs_expressed_all_samples.xls > $outdir/explevel_analysis/known/matrix_miRNA.count.xls\n";
			print FOUT "cp $outdir/identification/mirdeep/miRNAs_expressed_all_samples_normalization.xls $outdir/explevel_analysis/known/matrix_miRNA.tpm.xls\n";
			print FOUT "cd $outdir/explevel_analysis/known/\n";
			print FOUT "Rscript $script/exp_dir/top_ten_miRNA.r --args matrix_miRNA.tpm.xls\n";
		}else{
			print FOUT "sed 1d $outdir/identification/mirdeep/predicted_total_sumary_table.xls |sed '\$d' |sed '\$d' |sed 's/#provisional_id/miRNA_ID/' | cut -f 1,3- > $outdir/explevel_analysis/novel/matrix_miRNA.count.xls\n";
			print FOUT "sed 's/#provisional_id/miRNA_ID/' $outdir/identification/mirdeep/novo_mature_normalization.xls > $outdir/explevel_analysis/novel/matrix_miRNA.tpm.xls\n";
			print FOUT "cd $outdir/explevel_analysis/novel/\n";
			print FOUT "Rscript $script/exp_dir/top_ten_miRNA.r --args matrix_miRNA.tpm.xls\n";
		}
		close FOUT;
		system("sh $workdir/cmd/run.explevel_${xxx}.matrix.sh > $workdir/cmd/run.explevel_${xxx}.matrix.o 2> $workdir/cmd/run.explevel_${xxx}.matrix.e ");
	}
}
