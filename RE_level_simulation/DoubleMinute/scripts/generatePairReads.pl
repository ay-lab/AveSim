@array = (1..15);
@choose = qw(1 0);
$input_a = @ARGV[0]; #IMR90.$chr[$i].Normal.re.shuffled.interaction.txt
$output_a = @ARGV[1]; #IMR90.$chr[$i].Normal.re.shuffled.bed
$genome_name = @ARGV[2]; # chr13.fasta
$fastq_id = @ARGV[3]; #Any number or string to recognize fastq files
$chrom_name = @ARGV[4];
$bedtools_path = @ARGV[5];
chomp ($input_a,$output_a,$genome_name,$fastq_id,$chrom_name,$bedtools_path);
open (in,"$input_a");
open (out,">$output_a");
while (<in>){
	chomp $_;
	$randomAdd = $array[rand @array];
	$randomChoose = $choose[rand @choose];
	@data = split(/\s+/,$_);
	if (@data[1] > 200){
	if ($randomChoose == 1){
		print out @data[0],"\t",@data[1]-(23+$randomAdd),"\t",@data[1]+5,"\t",$.,"\n";
		print out @data[0],"\t",@data[2]+1,"\t",@data[2]+(23+(@arra[$#array]-$randomAdd)),"\t",$.,"\n";
		push @number,"$.\n";
	}else {
		print out @data[0],"\t",@data[1]-(23+(@arra[$#array]-$randomAdd)),"\t",@data[1]+5,"\t",$.,"\n";
             	print out @data[0],"\t",@data[2]+1,"\t",@data[2]+(23+$randomAdd),"\t",$.,"\n";
		push @number,"$.\n";
	}}
	undef @data;
}
close in;
close out;
`$bedtools_path nuc -fi $genome_name -bed $output_a -seq|grep -v "#" > $output_a.with_seq.bed`;
open (seq_in,"$output_a.with_seq.bed");
$c = 0;
while (<seq_in>){
	chomp $_;
	@data = split(/\s+/,$_);
	chomp @data[13];
	@data[13] =~ s/\s+//g;
	$pair{@data[3]} += 1;
	$seq{@data[3]}{$pair{@data[3]}} = @data[13];
	$c++; 	
	undef @data;
}
close seq_in;
open (pairA,">$output_a.re_R1.fastq");
open (pairB,">$output_a.re_R2.fastq");
$j = 1;
while ($j <= $#number){
	chomp $number[$j];	
	$job_id = int(rand($$*rand(10)))+int(rand($$*rand(10)))+int(rand($$*rand(10)));
	$s = $seq{$number[$j]}{1}.$seq{$number[$j]}{2};
	$k = 0;
	while ($k < length($s)){
       		$array .="2";
       		$k++;
	}
	print pairA "\@seq_$job_id\_$number[$j]:$fastq_id:0_7:0:0_0/1\n";
	print pairA "$s\n";
	print pairA "+\n";
	print pairA "$array\n";
	$s =~ tr /atcgATCG/tagcTAGC/;
	$s = reverse($s);	

	print pairB "\@seq_$job_id\_$number[$j]:$fastq_id:0_7:0:0_0/2\n";
        print pairB "$s\n";
        print pairB "+\n";
        print pairB "$array\n";	
	$array = undef;
	$j++;	
}
close pairA;
close pairB;
undef @number;
undef %seq;
undef %pair;
undef @array;
undef @choose;
close list;
