$stamp = `date|awk '{print \$3_\$2_\$6}'`;
$job   = $$;
chomp ($stamp,$job);

$re_cut     = "HindIII.fa";
$oligoMatch = "/mnt/BioApps/UCSC/oligoMatch";
$bedtools   = "/mnt/BioApps/bedtools/bin/bedtools";
@chr_list   = `cat ../scripts/chr.list`;
$chrom_num  = 0; map{chomp $_; $chrom_num++; $chrom_id{$_} = $chrom_num;}@chr_list;
@param_file = `cat amp.simulation.config.txt`;
$i = 0;
while ($i <= $#param_file){
	chomp $param_file[$i];
	$param_file[$i] =~ s/\s+//g;
	@paramData = split(/=/,$param_file[$i]);
	$param{@paramData[0]} = @paramData[1];
	undef @paramData;
	$i++;
}
$amplification_reads = int($param{'reads'}*$param{'hetero'});
$normal_reads   = int($param{'reads'}*(1-$param{'hetero'}));
$amp_put_end = $param{'amp_put'}+1;

open (out,">Amplification_Simulation_$param{'amp_chrom'}\_$param{'amp_start'}\_$param{'amp_end'}\_$stamp.$job.sh");
print out "#HiCnv script to call CNVs from Hi-C data\n";
print out "#Job performed on $stamp; Job id: $job\n\n";
print out "Rscript ../scripts/duplication.r $param{'amp_chrom'} $param{'amp_chrom'} $param{'amp_start'} $param{'amp_end'}  $param{'amp_put'} $amp_put_end $param{'cnv_copy'} $param{'sim_rest'}\n\n";
if ($param{'sim_rest'} eq "yes"){
	map {
		chomp $_;
		print out "$oligoMatch ../scripts/$re_cut Normal.$_.genome.fasta Normal.$_.genome.$re_cut\_cut.bed\n";
		print out "Rscript ../scripts/calculate.distance.r Normal.$_.genome.$re_cut\_cut.bed Normal $chrom_id{$_}\n";
		if ($_ ne $param{'amp_chrom'}){
			print out "shuf -n $param{'reads'} IMR90.$_.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$_.Normal.re.shuffled.interaction.txt\n";	
		}
		else {
			print out "shuf -n $normal_reads IMR90.$_.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$_.Normal.re.shuffled.interaction.txt\n";
		}
		print out "perl ../scripts/generatePairReads.pl IMR90.$_.Normal.re.shuffled.interaction.txt IMR90.$_.Normal.re.shuffled.bed Normal.$_.genome.fasta Normal $_ $bedtools\n";
	}@chr_list;
	print out "$oligoMatch ../scripts/$re_cut Duplicated.$_.genome.fasta Duplicated.$_.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Duplicated.$_.genome.$re_cut\_cut.bed Duplicated $chrom_id{$_}\n";
	print out "shuf -n $amplification_reads IMR90.$_.Duplicated.re.interaction.txt|sort -g -k 2 > IMR90.$_.Duplicated.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$_.Duplicated.re.shuffled.interaction.txt IMR90.$_.Duplicated.re.shuffled.bed Duplicated.$_.genome.fasta Duplicated $_ $bedtools\n\n";
}
else {
	print out "$oligoMatch ../scripts/$re_cut Normal.$param{'amp_chrom'}.genome.fasta Normal.$param{'amp_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Normal.$param{'amp_chrom'}.genome.$re_cut\_cut.bed Normal $chrom_id{$param{'amp_chrom'}}\n";
	print out "shuf -n $normal_reads IMR90.$param{'amp_chrom'}.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$param{'amp_chrom'}.Normal.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'amp_chrom'}.Normal.re.shuffled.interaction.txt IMR90.$param{'amp_chrom'}.Normal.re.shuffled.bed Normal.$param{'amp_chrom'}.genome.fasta Normal $param{'amp_chrom'} $bedtools\n";
	print out "$oligoMatch ../scripts/$re_cut Duplicated.$param{'amp_chrom'}.genome.fasta Duplicated.$param{'amp_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Duplicated.$param{'amp_chrom'}.genome.$re_cut\_cut.bed Duplicated $chrom_id{$param{'amp_chrom'}}\n";
	print out "shuf -n $amplification_reads IMR90.$param{'amp_chrom'}.Duplicated.re.interaction.txt|sort -g -k 2 > IMR90.$param{'amp_chrom'}.Duplicated.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'amp_chrom'}.Duplicated.re.shuffled.interaction.txt IMR90.$param{'amp_chrom'}.Duplicated.re.shuffled.bed Duplicated.$param{'amp_chrom'}.genome.fasta Duplicated $param{'amp_chrom'} $bedtools\n\n";
}
print out "cat *.re.shuffled.bed.re_R1.fastq > Amplification_Simulation_$param{'amp_chrom'}\_$param{'amp_start'}\_$param{'amp_end'}\_$stamp.$job.re_R1.fastq\n";
print out "cat *.re.shuffled.bed.re_R2.fastq > Amplification_Simulation_$param{'amp_chrom'}\_$param{'amp_start'}\_$param{'amp_end'}\_$stamp.$job.re_R2.fastq\n";
print out "mkdir FinalFastqs\n";
print out "mkdir IntermediateFiles\n";
print out "mv Amplification_Simulation_$param{'amp_chrom'}\_$param{'amp_start'}\_$param{'amp_end'}\_$stamp.$job.re_*.fastq FinalFastqs/\n";
print out "mv *.fastq IntermediateFiles/\n";
print out "mv *.bed IntermediateFiles/\n";
print out "mv *.interaction.txt IntermediateFiles/\n";
print out "mv *.fasta IntermediateFiles/\n";
print out "mv *.fai IntermediateFiles/\n";
close out;
`chmod 755 Amplification_Simulation_$param{'amp_chrom'}\_$param{'amp_start'}\_$param{'amp_end'}\_$stamp.$job.sh`;
undef %chrom_id;
undef %param;
undef @param_file;
