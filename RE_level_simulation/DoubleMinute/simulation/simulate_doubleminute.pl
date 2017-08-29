$stamp = `date|awk '{print \$3_\$2_\$6}'`;
$job   = $$;
chomp ($stamp,$job);

$re_cut     = "HindIII.fa";
$oligoMatch = "/mnt/BioApps/UCSC/oligoMatch";
$bedtools   = "/mnt/BioApps/bedtools/bin/bedtools";
@chr_list   = `cat ../scripts/chr.list`;
$chrom_num  = 0; map{chomp $_; $chrom_num++; $chrom_id{$_} = $chrom_num;}@chr_list;
@param_file = `cat dm.simulation.config.txt`;
$i = 0;
while ($i <= $#param_file){
	chomp $param_file[$i];
	$param_file[$i] =~ s/\s+//g;
	@paramData = split(/=/,$param_file[$i]);
	$param{@paramData[0]} = @paramData[1];
	undef @paramData;
	$i++;
}

open (out,">DoubleMinute_Simulation_$param{'dm_chrom'}\_$param{'dm_start'}\_$param{'dm_end'}\_$stamp.$job.sh");
print out "#HiCnv script to call CNVs from Hi-C data\n";
print out "#Job performed on $stamp; Job id: $job\n\n";
print out "Rscript ../scripts/normal.r\n\n";
map {
	chomp $_;
	print out "$oligoMatch ../scripts/$re_cut Normal.$_.genome.fasta Normal.$_.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Normal.$_.genome.$re_cut\_cut.bed Normal $chrom_id{$_}\n";
	print out "shuf -n $param{'reads'} IMR90.$_.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$_.Normal.re.shuffled.interaction.txt\n";
	if ($_ eq $param{'dm_chrom'}){
		print out "awk \'{if(\$2 >= $param{'dm_start'} && \$3 <= $param{'dm_end'}){print}}\' IMR90.$_.Normal.re.interaction.txt|shuf -n $param{'reads'}|sort -g -k 2 > DM.re.shuffled.interaction.txt\n";
		print out "var=\$(cat DM.re.shuffled.interaction.txt|wc -l)\n";
		print out "perl ../scripts/generatePairReads.pl DM.re.shuffled.interaction.txt DM.Normal.re.shuffled.bed Normal.$param{'dm_chrom'}.genome.fasta Normal $param{'dm_chrom'} $bedtools\n";
	}	
	print out "perl ../scripts/generatePairReads.pl IMR90.$_.Normal.re.shuffled.interaction.txt IMR90.$_.Normal.re.shuffled.bed Normal.$_.genome.fasta Normal $_ $bedtools\n";
}@chr_list;

map {
	chomp $_;
	print out "cat IMR90.$_.Normal.re.shuffled.bed.with_seq.bed |grep -ie \"AAGCT\$\"|shuf -n \$var|awk '{print \$1\"\\t\"\$2\"\\t\"\$3}' > $_\_Normal_$param{'dm_chrom'}_DM.seq_A.bed\n"; 
	print out "cat DM.Normal.re.shuffled.bed.with_seq.bed |grep -ie \"\\bAGCTT\" |shuf -n \$var|awk '{print \$1\"\\t\"\$2\"\\t\"\$3}' > $_\_Normal_$param{'dm_chrom'}_DM.seq_B.bed\n";
}@chr_list;
print out "perl ../scripts/DM_pairEndFileGenerate.pl $param{'dm_chrom'} $bedtools\n";
print out "cat IMR90.*.Normal.re.shuffled.bed.re_R1.fastq *_Normal_$param{'dm_chrom'}_DM.re_R1.fastq > Normal_DM.re_R1.fastq\n";
print out "cat IMR90.*.Normal.re.shuffled.bed.re_R2.fastq *_Normal_$param{'dm_chrom'}_DM.re_R2.fastq > Normal_DM.re_R2.fastq\n";
print out "mkdir FinalFastqs\n";
print out "mkdir IntermediateFiles\n";
print out "mv Normal_DM.re_*.fastq FinalFastqs/\n";
print out "mv *.fastq IntermediateFiles/\n";
print out "mv *.bed IntermediateFiles/\n";
print out "mv *.interaction.txt IntermediateFiles/\n";
print out "mv *.fasta IntermediateFiles/\n";
print out "mv *.fai IntermediateFiles/\n";
close out;
`chmod 755 DoubleMinute_Simulation_$param{'dm_chrom'}\_$param{'dm_start'}\_$param{'dm_end'}\_$stamp.$job.sh`;
undef %chrom_id;
undef %param;
undef @param_file;
