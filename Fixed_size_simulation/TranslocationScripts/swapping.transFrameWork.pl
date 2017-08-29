#author abhijit 
#created on 09/24/16
#This file will generate the frameWork files for translocation simulation
#
@fWork = `cat @ARGV[0]`;
open(a,">@ARGV[0].A");
open(b,">@ARGV[0].B");
$i = 0;
while ($i <= $#fWork)
{
	chomp $fWork[$i];
	if ($fWork[$i] =~ /^#/){
		print a "$fWork[$i]\n";
		print b "$fWork[$i]\n";
	}
	elsif ($fWork[$i] =~ /^bedFile/	|| $fWork[$i] =~ /^matFile/ || $fWork[$i] =~ /^chromosomes/){
		if ($fWork[$i] =~ /^chromosomes/){		
			@chrData = split(/\s+/,$fWork[$i]);
			$chrA_name = @chrData[1];
			$chrB_name = @chrData[2];
			undef @chrData;
		}
		if ($fWork[$i] =~ /^bedFile/){	
			
			@bedData = split(/\s+/,$fWork[$i]);
			`awk '{if(\$1==\"$chrA_name\" || \$1==\"$chrB_name\"){c++;print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"c}}' @bedData[1] > chr.bed`;
			undef @bedData;
		}
		if ($fWork[$i] =~ /^matFile/){	
			@data = split(/\s+/,$fWork[$i]);
			open (bed,"chr.bed");
			open (bed_out,">chr.mod.bed");
			while (<bed>){
				chomp $_;
				@bedData = split(/\s+/,$_);
				$bed_old_to_new{@bedData[3]} = @bedData[4];
				$bed_chr_name{@bedData[3]} = @bedData[0];
				print bed_out "@bedData[0]\t@bedData[1]\t@bedData[2]\t@bedData[4]\n";
				undef @bedData;
			}
			close bed;
			close bed_out;
			print a "bedFile\tchr.mod.bed\n";
			print a "matFile\tmat.mod.matrix\n";
			print b "bedFile\tchr.mod.bed\n";
                        print b "matFile\tmat.mod.matrix\n";
			undef @data;
		}
	}
	elsif ($fWork[$i] =~ /^chrA/){
		@data = split(/\s+/,$fWork[$i]);
		$chrA_name = @data[1];
		$Cell  = @data[2];
		$chrA_start = @data[3];		
		$chrA_end = @data[4];
		undef @data;
	}
	elsif ($fWork[$i] =~ /^chrB/){
                @data = split(/\s+/,$fWork[$i]);
                $chrB_name = @data[1];
                $chrB_start = @data[3];
                $chrB_end = @data[4];
                undef @data;
        }
	elsif ($fWork[$i] =~ /^New_Chromosome1/){
		@data  = split(/\s+/,$fWork[$i]);
		@mergingData = split(/\)/,@data[1]);
		$new_Chr_A_length = 0;
		print b "chr\t$chrB_name\n";
	        print b "output\t$Cell\n";
                print b ">TR:RMV:Start\n";
		$countB1 = 0;
		$k = 0;
		while ($k <= $#mergingData){
			chomp $mergingData[$k];	
			@partnerData = split(/~/,$mergingData[$k]);
			if ($#partnerData == 1){
				chomp (@partnerData[0],@partnerData[1]);
				@cisChr = split(/\:/,@partnerData[0]);
				@transChr = split(/\:/,@partnerData[1]);
				if (@cisChr[0] < @cisChr[1]){
					$m = @cisChr[0];
					while ($m <= @cisChr[1]){
						$new_Chr_A_length += 1;
						$new_Chr_A_bin_name{$new_Chr_A_length} = $chrA_name;
						$new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
						push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
							push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
						}
						if ($bed_chr_name{$m} eq $chrB_name){
							$countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
						$m++;
					}
				}	
				else {
					$m = @cisChr[0];
	                                while ($m >= @cisChr[1]){
	                                        $new_Chr_A_length += 1;
	                                        $new_Chr_A_bin_name{$new_Chr_A_length} = $chrA_name;
	                                        $new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
						push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        $countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
	                                        $m--;
	                                }
				}
				if (@transChr[0] < @transChr[1]){
					print b "@transChr[0]:@transChr[1]\n";
					$m = @transChr[0];
	                                while ($m <= @transChr[1]){
	                                        $new_Chr_A_length += 1;
        	                                $new_Chr_A_bin_name{$new_Chr_A_length} = $chrB_name;
        	                                $new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
        	                                push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        $countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
        	                                $m++;
        	                        }
        	                }
        	                else {
					print b "@transChr[1]:@transChr[0]\n";
					$m = @transChr[0];
        	                        while ($m >= @transChr[1]){
        	                                $new_Chr_A_length += 1;
        	                                $new_Chr_A_bin_name{$new_Chr_A_length} = $chrB_name;
        	                                $new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
        	                                push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        $countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
        	                                $m--;
        	                        }
        	                }
				undef @cisChr;
	                        undef @transChr;
			}
			else {
				@cisChr = split(/\:/,@partnerData[0]);
				if (@cisChr[0] < @cisChr[1]){
                                        $m = @cisChr[0];
                                        while ($m <= @cisChr[1]){
                                                $new_Chr_A_length += 1;
                                                $new_Chr_A_bin_name{$new_Chr_A_length} = $chrA_name;
                                                $new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
                                                push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        $countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
                                                $m++;
                                        }
                                }
                                else {
                                        $m = @cisChr[0];
                                        while ($m >= @cisChr[1]){
                                                $new_Chr_A_length += 1;
                                                $new_Chr_A_bin_name{$new_Chr_A_length} = $chrA_name;
                                                $new_Chr_A_new_To_old_bin_id{$new_Chr_A_length} = $m;
                                                push @new_chromosome_A,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$new_Chr_A_length\t+\t-\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        $countB1++;
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$countB1\t-\t+\n";
                                                }
                                                $m--;
                                        }
                                }			
				undef @cisChr;
			}
			undef @partnerData;
			$k++;
		}
		print b ">TR:RMV:End";
		undef @mergingData;
		undef @data;
	}
	elsif ($fWork[$i] =~ /^New_Chromosome2/){
		@data  = split(/\s+/,$fWork[$i]);
                @mergingData = split(/\)/,@data[1]);
                $new_Chr_B_length = 0;
                print a "chr\t$chrA_name\n";
                print a "output\t$Cell\n";
                print a ">TR:RMV:Start\n";
		$countA2 = 0;
                $k = 0;
                while ($k <= $#mergingData){
                        chomp $mergingData[$k];
                        @partnerData = split(/~/,$mergingData[$k]);
			if ($#partnerData == 1){
	                        chomp (@partnerData[0],@partnerData[1]);
	                        @cisChr = split(/\:/,@partnerData[0]);
	                        @transChr = split(/\:/,@partnerData[1]);
				$difference = abs(@transChr[0]-@transChr[1]);
	                        if (@cisChr[0] < @cisChr[1]){
	                                $m = @cisChr[0];
	                                while ($m <= @cisChr[1]){
	                                        $new_Chr_B_length += 1;
	                                        $new_Chr_B_bin_name{$new_Chr_B_length} = $chrB_name;
	                                        $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;
	                                        push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
	                                        $m++;
	                                }
	                        }
	                        else {
	                                $m = @cisChr[0];
	                                while ($m >= @cisChr[1]){
	                                        $new_Chr_B_length += 1;
	                                        $new_Chr_B_bin_name{$new_Chr_B_length} = $chrB_name;
	                                        $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;	
	                                        push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
	                                        $m--;
	                                }
	                        }
	                        if (@transChr[0] < @transChr[1]){
	                                print a "@transChr[0]:@transChr[1]\n";
	                                $m = @transChr[0];
	                                while ($m <= @transChr[1]){
	                                        $new_Chr_B_length += 1;
	                                        $new_Chr_B_bin_name{$new_Chr_B_length} = $chrA_name;
	                                        $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;
	                                        push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
	                                        $m++;
	                                }
	                        }
	                        else {
	                                print a "@transChr[1]:@transChr[0]\n";
	                                $m = @transChr[0];
	                                while ($m >= @transChr[1]){
	                                        $new_Chr_B_length += 1;
	                                        $new_Chr_B_bin_name{$new_Chr_B_length} = $chrA_name;
	                                        $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;
        	                                push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
        	                                $m--;
        	                        }
        	                }
                	        undef @cisChr;
                	        undef @transChr;
			}
			else {
				@cisChr = split(/\:/,@partnerData[0]);
				if (@cisChr[0] < @cisChr[1]){
                                        $m = @cisChr[0];
                                        while ($m <= @cisChr[1]){
                                                $new_Chr_B_length += 1;
                                                $new_Chr_B_bin_name{$new_Chr_B_length} = $chrB_name;
                                                $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;
                                                push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
                                                $m++;
                                        }
                                }
                                else {
                                        $m = @cisChr[0];
                                        while ($m >= @cisChr[1]){
                                                $new_Chr_B_length += 1;
                                                $new_Chr_B_bin_name{$new_Chr_B_length} = $chrB_name;
                                                $new_Chr_B_new_To_old_bin_id{$new_Chr_B_length} = $m;
                                                push @new_chromosome_B,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\tNoraml\tNormal\n";
						if ($bed_chr_name{$m} eq $chrA_name){
                                                        $countA2++;
                                                        push @chrA_index,"$bed_chr_name{$m}\t$m\t$countA2\t-\t+\n";
                                                }
                                                if ($bed_chr_name{$m} eq $chrB_name){
                                                        push @chrB_index,"$bed_chr_name{$m}\t$m\t$new_Chr_B_length\t+\t-\n";
                                                }
                                                $m--;
                                        }
                                }
				undef @cisChr;
			}
                        undef @partnerData;
                        $k++;
                }
                print a ">TR:RMV:End";
                undef @mergingData;
                undef @data;
        }
	$i++;
}
close a;
close b;

open (chrA_out,">$chrA_name.modified.index");
open (chrB_out,">$chrB_name.modified.index");
print chrA_out @chrA_index;
print chrB_out @chrB_index;
close chrA_out;
close chrB_out;
undef @chrA_index;
undef @chrB_index;

print @new_chromosome_A;

open(chrom,">$Cell.$chrA_name-$chrB_name.Simulation.matrix");
$i = 0;
while ($i <= $#new_chromosome_A){
	chomp $new_chromosome_A[$i];
	@chrData1 = split(/\s+/,$new_chromosome_A[$i]);
	$j = $i+1;
	while ($j <= $#new_chromosome_A){
		chomp $new_chromosome_A[$j];
		@chrData2 = split(/\s+/,$new_chromosome_A[$j]);
		if (@chrData1[0] eq $chrA_name && @chrData2[0] eq $chrB_name){
			print chrom "@chrData1[1]\t@chrData2[1]\t@chrData1[2]\t@chrData2[2]\t1\tNRM\tNRM\t-\t-\n";
		}
		if (@chrData1[0] eq $chrB_name && @chrData2[0] eq $chrA_name){
                        print chrom "@chrData2[1]\t@chrData1[1]\t@chrData2[2]\t@chrData1[2]\t1\tNRM\tNRM\t-\t-\n";
                }
		undef @chrData2;
		$j++;
	}
	undef @chrData1;
	$i++;
}

$i = 0;
while ($i <= $#new_chromosome_B){
        chomp $new_chromosome_B[$i];
        @chrData1 = split(/\s+/,$new_chromosome_B[$i]);
        $j = $i+1;
        while ($j <= $#new_chromosome_B){
                chomp $new_chromosome_B[$j];
                @chrData2 = split(/\s+/,$new_chromosome_B[$j]);
                if (@chrData1[0] eq $chrA_name && @chrData2[0] eq $chrB_name){
                        print chrom "@chrData1[1]\t@chrData2[1]\t@chrData1[2]\t@chrData2[2]\t1\tNRM\tNRM\t-\t-\n";
                }
                if (@chrData1[0] eq $chrB_name && @chrData2[0] eq $chrA_name){
                        print chrom "@chrData2[1]\t@chrData1[1]\t@chrData2[2]\t@chrData1[2]\t1\tNRM\tNRM\t-\t-\n";
                }
                undef @chrData2;
                $j++;
        }
        undef @chrData1;
        $i++;
}

$i = 0;
while ($i <= $#new_chromosome_A){
        chomp $new_chromosome_A[$i];
        @chrData1 = split(/\s+/,$new_chromosome_A[$i]);
        $j = $i+1;
        while ($j <= $#new_chromosome_B){
                chomp $new_chromosome_B[$j];
                @chrData2 = split(/\s+/,$new_chromosome_B[$j]);
                if (@chrData1[0] eq $chrA_name && @chrData2[0] eq $chrB_name){
                        #print chrom "@chrData1[1]\t@chrData2[1]\t1\t1000\t1\tNRM\tNRM\t-\t-\n";
                }
                if (@chrData1[0] eq $chrB_name && @chrData2[0] eq $chrA_name){
                        #print chrom "@chrData2[1]\t@chrData1[1]\t1\t1000\t1\tNRM\tNRM\t-\t-\n";
                }
                undef @chrData2;
                $j++;
        }
        undef @chrData1;
        $i++;
}
close chrom
