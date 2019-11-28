@chr = `cat ../scripts/chr.list`;
$fastq_id = "DM_MIXED";
$dm_chrom = @ARGV[0];
$bedtools_path = @ARGV[1];
chomp ($dm_chrom, $bedtools_path);
$i = 0;
while ($i <= $#chr){
	chomp $chr[$i];	
	`paste $chr[$i]_Normal_$dm_chrom\_DM.seq_A.bed $chr[$i]_Normal_$dm_chrom\_DM.seq_B.bed > $chr[$i]_Normal_$dm_chrom\_DM.seq.temp`;
        open (temp_out_a,">$chr[$i]_Normal_$dm_chrom\_DM.seq_temp.A.bed");
        open (temp_out_b,">$chr[$i]_Normal_$dm_chrom\_DM.seq_temp.B.bed");
        open (temp_in,"$chr[$i]_Normal_$dm_chrom\_DM.seq.temp");
        while (<temp_in>){
                chomp $_;
                @data = split(/\s+/,$_);
                if ((@data[2]-@data[1]) > 25){
                        $subtract = (@data[2]-@data[1])-25;
                        print temp_out_a @data[0],"\t",@data[1]+$subtract,"\t",@data[2],"\t$.\n";
                }
                elsif ((@data[2]-@data[1]) < 25){
                        $add = 25-(@data[2]-@data[1]);
                        if ((@data[1]-$add) > 0){
                                print temp_out_a @data[0],"\t",@data[1]-$add,"\t",@data[2],"\t$.\n";
                        }
                        else {
                                print temp_out_a @data[0],"\t",@data[1],"\t",@data[2],"\t$.\n";
                        }
                }
                elsif ((@data[2]-@data[1]) == 25){
                        print temp_out_a @data[0],"\t",@data[1],"\t",@data[2],"\t$.\n";
                }
                if ((@data[5]-@data[4]) > 25){
                        $subtract = (@data[5]-@data[4])-25;
                        print temp_out_b @data[3],"\t",@data[4],"\t",@data[5]-$subtract,"\t$.\n";
                }
                elsif ((@data[5]-@data[4]) < 25){
                        $add = 25-(@data[5]-@data[4]);
                        print temp_out_b @data[3],"\t",@data[4],"\t",@data[5]+$add,"\t$.\n";
                }
                elsif ((@data[5]-@data[4]) == 25){
                        print temp_out_b @data[3],"\t",@data[4],"\t",@data[5],"\t$.\n";
                }
                undef @data;
        }
        close temp_in;
        close temp_out_a;
        close temp_out_b;
        `$bedtools_path nuc -fi Normal.$chr[$i].genome.fasta -bed $chr[$i]_Normal_$dm_chrom\_DM.seq_temp.A.bed -seq|grep -v "#" > $chr[$i]_Normal_$dm_chrom\_DM.seq_A.bed`;
        `$bedtools_path nuc -fi Normal.$dm_chrom.genome.fasta -bed $chr[$i]_Normal_$dm_chrom\_DM.seq_temp.B.bed -seq|grep -v "#"|sed 's/$dm_chrom\\b/$chr[$i]/g' > $chr[$i]_Normal_$dm_chrom\_DM.seq_B.bed`;
        `cat $chr[$i]_Normal_$dm_chrom\_DM.seq_A.bed $chr[$i]_Normal_$dm_chrom\_DM.seq_B.bed > $chr[$i]_Normal_$dm_chrom\_DM.seq.bed`;
	open (seq_in,"$chr[$i]_Normal_$dm_chrom\_DM.seq.bed");
        $c = 0;
        while (<seq_in>){
                chomp $_;
                @data = split(/\s+/,$_);
                chomp @data[13];
                @data[13] =~ s/\s+//g;
                $pair{@data[3]} += 1;
                $seq{@data[3]}{$pair{@data[3]}} = @data[13];
		if($check{@data[3]} eq ""){
			push @number,"@data[3]\n";
			$check{@data[3]} = "yes";
		}
                $c++;
                undef @data;
        }
        close seq_in;
        open (pairA,">$chr[$i]_Normal_$dm_chrom\_DM.re_R1.fastq");
        open (pairB,">$chr[$i]_Normal_$dm_chrom\_DM.re_R2.fastq");
        $j = 0;
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
	undef %pair;
	undef %seq;
	undef %check;
	undef @number;
	$i++;
}
