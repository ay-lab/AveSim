$stamp = `date|awk '{print \$3_\$2_\$6}'`;
$job   = $$;
chomp ($stamp,$job);

$re_cut     = "HindIII.fa";
$oligoMatch = "/mnt/BioApps/UCSC/oligoMatch";
$bedtools   = "/mnt/BioApps/bedtools/bin/bedtools";
@chr_list   = `cat ../scripts/chr.list`;
$chrom_num  = 0; map{chomp $_; $chrom_num++; $chrom_id{$_} = $chrom_num;}@chr_list;
@param_file = `cat trans.simulation.config.txt`;
$i = 0;
while ($i <= $#param_file){
	chomp $param_file[$i];
	$param_file[$i] =~ s/\s+//g;
	@paramData = split(/=/,$param_file[$i]);
	$param{@paramData[0]} = @paramData[1];
	undef @paramData;
	$i++;
}
$trans_reads = int($param{'reads'}*$param{'hetero'});
$normal_reads  = int($param{'reads'}*(1-$param{'hetero'}));

open (out,">Translocation_Simulation_$param{'don_chrom'}\_$param{'don_start'}\_$param{'don_end'}\_$stamp.$job.sh");
print out "#HiCnv script to call CNVs from Hi-C data\n";
print out "#Job performed on $stamp; Job id: $job\n\n";
print out "Rscript ../scripts/translocation.r $param{'don_chrom'} $param{'rec_chrom'} $param{'don_start'} $param{'don_end'} $param{'rec_put'} $param{'sim_rest'}\n\n";
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
	print out "$oligoMatch ../scripts/$re_cut Translocation_Donor.$param{'don_chrom'}.genome.fasta Translocation_Donor.$param{'don_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "$oligoMatch ../scripts/$re_cut Translocation_Receiver.$param{'rec_chrom'}.genome.fasta Translocation_Receiver.$param{'rec_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Translocation_Donor.$param{'don_chrom'}.genome.$re_cut\_cut.bed Translocation_Donor $chrom_id{$param{'don_chrom'}}\n";
	print out "Rscript ../scripts/calculate.distance.r Translocation_Receiver.$param{'rec_chrom'}.genome.$re_cut\_cut.bed Translocation_Receiver $chrom_id{$param{'rec_chrom'}}\n";
	print out "shuf -n $trans_reads IMR90.$param{'don_chrom'}.Translocation_Donor.re.interaction.txt|sort -g -k 2 > IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.interaction.txt\n";
	print out "shuf -n $trans_reads IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.interaction.txt|sort -g -k 2 > IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.interaction.txt IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.bed Translocation_Donor.$param{'don_chrom'}.genome.fasta Translocation_Donor $param{'don_chrom'} $bedtools\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.interaction.txt IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.bed Translocation_Receiver.$param{'rec_chrom'}.genome.fasta Translocation_Receiver $param{'rec_chrom'} $bedtools\n\n";

}
else {
	print out "$oligoMatch ../scripts/$re_cut Normal.$param{'don_chrom'}.genome.fasta Normal.$param{'don_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "$oligoMatch ../scripts/$re_cut Normal.$param{'rec_chrom'}.genome.fasta Normal.$param{'rec_chrom'}.genome.$re_cut\_cut.bed\n";
	print out "Rscript ../scripts/calculate.distance.r Normal.$param{'don_chrom'}.genome.$re_cut\_cut.bed Normal $chrom_id{$param{'don_chrom'}}\n";
	print out "Rscript ../scripts/calculate.distance.r Normal.$param{'rec_chrom'}.genome.$re_cut\_cut.bed Normal $chrom_id{$param{'rec_chrom'}}\n";
	print out "shuf -n $normal_reads IMR90.$param{'don_chrom'}.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$param{'don_chrom'}.Normal.re.shuffled.interaction.txt\n";
	print out "shuf -n $normal_reads IMR90.$param{'rec_chrom'}.Normal.re.interaction.txt|sort -g -k 2 > IMR90.$param{'rec_chrom'}.Normal.re.shuffled.interaction.txt\n";
	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'don_chrom'}.Normal.re.shuffled.interaction.txt IMR90.$param{'don_chrom'}.Normal.re.shuffled.bed Normal.$param{'don_chrom'}.genome.fasta Normal $param{'don_chrom'} $bedtools\n";
 	print out "perl ../scripts/generatePairReads.pl IMR90.$param{'rec_chrom'}.Normal.re.shuffled.interaction.txt IMR90.$param{'rec_chrom'}.Normal.re.shuffled.bed Normal.$param{'rec_chrom'}.genome.fasta Normal $param{'rec_chrom'} $bedtools\n";
	print out "$oligoMatch ../scripts/$re_cut Translocation_Donor.$param{'don_chrom'}.genome.fasta Translocation_Donor.$param{'don_chrom'}.genome.$re_cut\_cut.bed\n";
        print out "$oligoMatch ../scripts/$re_cut Translocation_Receiver.$param{'rec_chrom'}.genome.fasta Translocation_Receiver.$param{'rec_chrom'}.genome.$re_cut\_cut.bed\n";
        print out "Rscript ../scripts/calculate.distance.r Translocation_Donor.$param{'don_chrom'}.genome.$re_cut\_cut.bed Translocation_Donor $chrom_id{$param{'don_chrom'}}\n";
        print out "Rscript ../scripts/calculate.distance.r Translocation_Receiver.$param{'rec_chrom'}.genome.$re_cut\_cut.bed Translocation_Receiver $chrom_id{$param{'rec_chrom'}}\n";
        print out "shuf -n $trans_reads IMR90.$param{'don_chrom'}.Translocation_Donor.re.interaction.txt|sort -g -k 2 > IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.interaction.txt\n";
        print out "shuf -n $trans_reads IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.interaction.txt|sort -g -k 2 > IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.interaction.txt\n";
        print out "perl ../scripts/generatePairReads.pl IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.interaction.txt IMR90.$param{'don_chrom'}.Translocation_Donor.re.shuffled.bed Translocation_Donor.$param{'don_chrom'}.genome.fasta Translocation_Donor $param{'don_chrom'} $bedtools\n";
        print out "perl ../scripts/generatePairReads.pl IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.interaction.txt IMR90.$param{'rec_chrom'}.Translocation_Receiver.re.shuffled.bed Translocation_Receiver.$param{'rec_chrom'}.genome.fasta Translocation_Receiver $param{'rec_chrom'} $bedtools\n\n";
}

print out "cat *.re.shuffled.bed.re_R1.fastq > Translocation_Simulation_$param{'don_chrom'}\_$param{'don_start'}\_$param{'don_end'}\_$stamp.$job.re_R1.fastq\n";
print out "cat *.re.shuffled.bed.re_R2.fastq > Translocation_Simulation_$param{'don_chrom'}\_$param{'don_start'}\_$param{'don_end'}\_$stamp.$job.re_R2.fastq\n";
print out "mkdir FinalFastqs\n";
print out "mkdir IntermediateFiles\n";
print out "mv Translocation_Simulation_$param{'don_chrom'}\_$param{'don_start'}\_$param{'don_end'}\_$stamp.$job.re_*.fastq FinalFastqs/\n";
print out "mv *.fastq IntermediateFiles/\n";
print out "mv *.bed IntermediateFiles/\n";
print out "mv *.interaction.txt IntermediateFiles/\n";
print out "mv *.fasta IntermediateFiles/\n";
print out "mv *.fai IntermediateFiles/\n";
close out;

`chmod 755 Translocation_Simulation_$param{'don_chrom'}\_$param{'don_start'}\_$param{'don_end'}\_$stamp.$job.sh`;
undef %chrom_id;
undef %param;
undef @param_file;
