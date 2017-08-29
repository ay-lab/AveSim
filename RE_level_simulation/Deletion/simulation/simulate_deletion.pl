$stamp = `date|awk '{print \$3_\$2_\$6}'`;
$job   = $$;
chomp ($stamp,$job);

$re_cut     = "HindIII.fa";
$oligoMatch = "/mnt/BioApps/UCSC/oligoMatch";
$bedtools   = "/mnt/BioApps/bedtools/bin/bedtools";
@chr_list   = `cat ../scripts/chr.list`;
$chrom_num  = 0; map{chomp $_; $chrom_num++; $chrom_id{$_} = $chrom_num;}@chr_list;
@param_file = `cat del.simulation.config.txt`;
$i = 0;
while ($i <= $#param_file){
	chomp $param_file[$i];
	$param_file[$i] =~ s/\s+//g;
	@paramData = split(/=/,$param_file[$i]);
	$param{@paramData[0]} = @paramData[1];
	undef @paramData;
	$i++;
}
$deletion_reads = int($param{'reads'}*$param{'hetero'});
$normal_reads   = int($param{'reads'}*(1-$param{'hetero'}));

open (out,">Deletion_Simulation_$param{'del_chrom'}\_$param{'del_start'}\_$param{'del_end'}\_$stamp.$job.sh");
print out "#HiCnv script to call CNVs from Hi-C data\n";
print out "#Job performed on $stamp; Job id: $job\n\n";
print out "Rscript ../scripts/deletion.r $param{'del_chrom'} $param{'del_start'} $param{'del_end'} $param{'sim_rest'}\n\n";
if ($param{'sim_rest'} eq "yes"){
	map {
		chomp $_;
		print out "$oligoMatch ../scripts/$re_cut Normal.$_.genome.fasta Normal.$_.genome.$re_cut\_cut.bed\n";
		print out "Rscript ../scripts/calculate.distance.r Normal.$_.genome.$re_cut\_cut.bed Normal $chrom_id{$_}\n";
		if ($_ ne $param{'del_chrom'}){
			print out "shuf -n $param{'reads'} IMR90.$_.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$_.Normal.re.shuffled.interaction.txt\n";	
		}
		else {
			print out "shuf -n $normal_reads IMR90.$_.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$_.Normal.re.shuffled.interaction.txt\n";
		}
		print out "perl ../scripts/generatePairReads.pl IMR90.$_.Normal.re.shuffled.interaction.txt IMR90.$_.Normal.re.shuffled.bed Normal.$_.genome.fasta Normal $_ $bedtools\n";
	}@chr_list;
	print out "$oligoMatch ../scripts/$re_cut Deleted.$_.genome.fasta Deleted.$_.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Deleted.$_.genome.$re_cut\_cut.bed Deleted $chrom_id{$_}\n";
	print out "shuf -n $deletion_reads IMR90.$_.Deleted.re.interaction.txt|sort -g -k 2 > IMR90.$_.Deleted.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$_.Deleted.re.shuffled.interaction.txt IMR90.$_.Deleted.re.shuffled.bed Deleted.$_.genome.fasta Deleted $_ $bedtools\n\n";
}
else {
	print out "$oligoMatch ../scripts/$re_cut Normal.$param{'del_chrom'}.genome.fasta Normal.$param{'del_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Normal.$param{'del_chrom'}.genome.$re_cut\_cut.bed Normal $chrom_id{$param{'del_chrom'}}\n";
	print out "shuf -n $normal_reads IMR90.$param{'del_chrom'}.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$param{'del_chrom'}.Normal.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'del_chrom'}.Normal.re.shuffled.interaction.txt IMR90.$param{'del_chrom'}.Normal.re.shuffled.bed Normal.$param{'del_chrom'}.genome.fasta Normal $param{'del_chrom'} $bedtools\n";
	print out "$oligoMatch ../scripts/$re_cut Deleted.$param{'del_chrom'}.genome.fasta Deleted.$param{'del_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Deleted.$param{'del_chrom'}.genome.$re_cut\_cut.bed Deleted $chrom_id{$param{'del_chrom'}}\n";
	print out "shuf -n $deletion_reads IMR90.$param{'del_chrom'}.Deleted.re.interaction.txt|sort -g -k 2 > IMR90.$param{'del_chrom'}.Deleted.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'del_chrom'}.Deleted.re.shuffled.interaction.txt IMR90.$param{'del_chrom'}.Deleted.re.shuffled.bed Deleted.$param{'del_chrom'}.genome.fasta Deleted $param{'del_chrom'} $bedtools\n\n";
}
print out "cat *.re.shuffled.bed.re_R1.fastq > Deletion_Simulation_$param{'del_chrom'}\_$param{'del_start'}\_$param{'del_end'}\_$stamp.$job.re_R1.fastq\n";
print out "cat *.re.shuffled.bed.re_R2.fastq > Deletion_Simulation_$param{'del_chrom'}\_$param{'del_start'}\_$param{'del_end'}\_$stamp.$job.re_R2.fastq\n";
print out "mkdir FinalFastqs\n";
print out "mkdir IntermediateFiles\n";
print out "mv Deletion_Simulation_$param{'del_chrom'}\_$param{'del_start'}\_$param{'del_end'}\_$stamp.$job.re_*.fastq FinalFastqs/\n";
print out "mv *.fastq IntermediateFiles/\n";
print out "mv *.bed IntermediateFiles/\n";
print out "mv *.interaction.txt IntermediateFiles/\n";
print out "mv *.fasta IntermediateFiles/\n";
print out "mv *.fai IntermediateFiles/\n";
close out;

`chmod 755 Deletion_Simulation_$param{'del_chrom'}\_$param{'del_start'}\_$param{'del_end'}\_$stamp.$job.sh`;
undef %chrom_id;
undef %param;
undef @param_file;
