# author abhijit
# created August 8th, 2016
#
#@frameWork = `cat frameWork.txt`;
@frameWork = `cat @ARGV[0]`;
$simNum = 0;
$check_inv = 0;
$check_trAdd = 0;
$check_trRmv = 0;
$i = 0;
while ($i <= $#frameWork)
{
	chomp $frameWork[$i];
	if ($frameWork[$i] =~ /^#/)
	{
		$simNum++;
		@frameDataSim = split(/\s+/,$frameWork[$i]);
		@frameDataBed = split(/\s+/,$frameWork[$i+1]);
		@frameDataMat = split(/\s+/,$frameWork[$i+2]);
		@frameDataChr = split(/\s+/,$frameWork[$i+4]);
		@frameDataOut = split(/\s+/,$frameWork[$i+3]);
		$bed{@frameDataSim[1]} = @frameDataBed[1];
		$mat{@frameDataSim[1]} = @frameDataMat[1];
		$chr{@frameDataSim[1]} = @frameDataChr[1];
		$out{@frameDataSim[1]} = @frameDataOut[1];
		undef @frameDataSim;
		undef @frameDataBed;
		undef @frameDataMat;
		undef @frameDataChr;
		undef @frameDataOut;
	}
	elsif ($frameWork[$i] =~ /^>AMP:Start/)
	{
		$j = $i+1;
		if ($frameWork[$j] !~ /^>AMP:End/)
                {
			do
			{
				chomp $frameWork[$j];
				@bins = split(/:/,$frameWork[$j]);
				$k = @bins[0];	
				while ($k <= @bins[1])
				{
					$cnv{$k} = "AMP";
					$k++;
				}
				undef @bins;
				$j++;
			}
			while ($frameWork[$j] !~ /^>AMP:End/);
		}
	}
	elsif ($frameWork[$i] =~ /^>HTD:Start/)
        {
                $j = $i+1;
		if ($frameWork[$j] !~ /^>HTD:End/)
		{
                	do
                	{
                	        chomp $frameWork[$j];
                	        @bins = split(/:/,$frameWork[$j]);
                	        $k = @bins[0];
                	        while ($k <= @bins[1])
                	        {
                	                $cnv{$k} = "HTD";
					$k++;
                	        }
                	        undef @bins;
                	        $j++;
                	}
                	while ($frameWork[$j] !~ /^>HTD:End/);
		}
        }
	elsif ($frameWork[$i] =~ /^>HMD:Start/)
        {
                $j = $i+1;
		if ($frameWork[$j] !~ /^>HMD:End/)
                {
                	do
                	{
                	        chomp $frameWork[$j];
                	        @bins = split(/:/,$frameWork[$j]);
                	        $k = @bins[0];
                	        while ($k <= @bins[1])
                	        {
                	                $cnv{$k} = "HMD";
					$k++;
                	        }
                	        undef @bins;
                	        $j++;
                	}
                	while ($frameWork[$j] !~ /^>HMD:End/);
		}
        }
	elsif ($frameWork[$i] =~ /^>INV:Start/)
        {
                $j = $i+1;
		if ($frameWork[$j] !~ /^>INV:End/)
                {
			$check_inv = 1;
                	do
                	{
                	        chomp $frameWork[$j];
                	        @bins = split(/:/,$frameWork[$j]);
                	        $k = @bins[0];
				$v = @bins[1];
                	        while ($k <= @bins[1])
                	        {
                	                $inv{$k} = $v; 
                	                $k++;
					$v--;
                	        }
                	        undef @bins;
                	        $j++;
                	}
                	while ($frameWork[$j] !~ /^>INV:End/);
        	}
	}
	elsif ($frameWork[$i] =~ /^>TR:ADD:Start/)
        {
		$add = 0;
                $j = $i+1;
                if ($frameWork[$j] !~ /^>TR:ADD:End/)
                {
			$check_trAdd = 1;
                        do
                        {
                                chomp $frameWork[$j];
                                @bins = split(/:/,$frameWork[$j]);
				$add += @bins[1]; 
				push @start,"@bins[0]\n";
				push @add,"$add\n";	
				push @end,"@bins[2]\n";
				undef @bins;
                                $j++;
                        }
                        while ($frameWork[$j] !~ /^>TR:ADD:End/);
			
			$j1 = 0;
			while ($j1 <= $#start)
			{
				chomp $start[$j1];
				chomp $add[$j1];
				chomp $end[$j1];
				if ($j1 != $#start)
				{
					$k = $start[$j1] + 1;
					$end = $start[$j1+1];
					while ($k <= $end)
					{
						$trAdd{$k} = $add[$j1] + $k;
						$k++;
					}
				}
				elsif ($j1 == $#start)
				{
					$k = $start[$j1] + 1;
                                        $end = 7000;
                                        while ($k <= $end)
                                        {
                                                $trAdd{$k} = $add[$j1] + $k;
                                                $k++;
                                        }
				}
				$j1++;
			}
			undef @start;
			undef @add;
			undef @end;
                }
        }
	elsif ($frameWork[$i] =~ /^>TR:RMV:Start/)
        {
                $add = 0;
                $j = $i+1;
                if ($frameWork[$j] !~ /^>TR:RMV:End/)
                {
                        $check_trRmv = 1;
                        do
                        {
                                chomp $frameWork[$j];
                                @bins = split(/:/,$frameWork[$j]);
                                $add -= abs(@bins[1] - @bins[0]) - 1;
                                
				push @start,"@bins[0]\n";
                                push @add,"$add\n";
                                push @end,"@bins[1]\n";
				
				$k = @bins[0]+1;
				while ($k < @bins[1])
				{
					$trName{$k} = "TRRMV";
					$k++;
				}
				
                                undef @bins;
                                $j++;
                        }
                        while ($frameWork[$j] !~ /^>TR:RMV:End/);

                        $j1 = 0;
                        while ($j1 <= $#start)
                        {
                                chomp $start[$j1];
                                chomp $add[$j1];
                                chomp $end[$j1];
                                if ($j1 != $#start)
                                {
                                        $k = $end[$j1];
                                        $end = $start[$j1+1];
                                        while ($k <= $end)
                                        {
                                                $trRmv{$k} = $add[$j1] + $k;
                                                $k++;
                                        }
                                }
                                elsif ($j1 == $#start)
                                {
                                        $k = $end[$j1];
                                        $end = 7000;
                                        while ($k <= $end)
                                        {
                                                $trRmv{$k} = $add[$j1] + $k;
                                                $k++;
                                        }
                                }
                                $j1++;
                        }
                        undef @start;
                        undef @add;
                        undef @end;
                }
        }		
	$i++;
}

$i = 1;
while ($i <= $simNum)
{
	open (bed,"$bed{$i}");
	open (out,">$out{$i}.$chr{$i}.bed");
	while (<bed>)
	{
	        chomp $_;
	        @bedData = split(/\s+/,$_);
	        $bedInfo_chr{$i}{@bedData[3]} = @bedData[0];
	        $bedInfo_pos{$i}{@bedData[3]}{'s'} = @bedData[1];
	        $bedInfo_pos{$i}{@bedData[3]}{'e'} = @bedData[2];
	        $resolution{$i} = abs(@bedData[2] - @bedData[1]);
	        #$bedInfo_index_count{$i}{@bedData[0]} += 1;
	        $bedInfo_index_count{$i}{@bedData[0]} = @bedData[3];
	        $bedInfo_index_span{$i}{@bedData[0]}{$bedInfo_index_count{@bedData[0]}} = @bedData[3];
		#$bedInfo_index_conv_id{$i}{@bedData[0]}{@bedData[3]} = $bedInfo_index_count{$i}{@bedData[0]};
		$bedInfo_index_conv_id{$i}{@bedData[0]}{@bedData[3]} = @bedData[3];
		if ($chr{$i} eq $bedInfo_chr{$i}{@bedData[3]})
		{
			print out "$bedInfo_chr{$i}{@bedData[3]}\t$bedInfo_pos{$i}{@bedData[3]}{'s'}\t$bedInfo_pos{$i}{@bedData[3]}{'e'}\t$bedInfo_index_count{$i}{@bedData[0]}\n";
		}
	        undef @bedData;
	}
	close out;
	close bed;
	
	open (bed,"$mat{$i}");
        open (out,">$out{$i}.$chr{$i}.PreSimulation.matrix");
        while (<bed>)
        {
                chomp $_;
                @matData = split(/\s+/,$_);
		if ($bedInfo_chr{$i}{@matData[0]} eq $chr{$i} && $bedInfo_chr{$i}{@matData[1]} eq $chr{$i})	
		{	
			if ($cnv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} eq ""){ $cnv1 = "NRM"; }
			else {$cnv1 = $cnv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}}; }
			if ($cnv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} eq ""){ $cnv2 = "NRM"; }
                        else {$cnv2 = $cnv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}}; }
		
			if ($check_inv == 1 && $check_trAdd == 0 && $check_trRmv == 0)
			{	
				if ($inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} ne ""){$index1 = $inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}}; $tr1 = "-";}
				else {$index1 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}; $tr1 = "-";}
				if ($inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} ne ""){$index2 = $inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}}; $tr2 = "-";}
				else {$index2 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}; $tr2 = "-";}
			}
			elsif ($check_inv == 0 && $check_trAdd == 1 && $check_trRmv == 0)
			{
				if ($trAdd{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} ne ""){$index1 = $trAdd{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}}; $tr1 = "-";}
                                else {$index1 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}; $tr1 = "-";}
                                if ($trAdd{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} ne ""){$index2 = $trAdd{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}}; $tr2 = "-";}
                                else {$index2 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}; $tr2 = "-";}
			}
			elsif ($check_inv == 0 && $check_trAdd == 0 && $check_trRmv == 1)
                        {
                                if ($trRmv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} ne ""){$index1 = $trRmv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}};}
                                else {$index1 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]};}
                                if ($trRmv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} ne ""){$index2 = $trRmv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}};}
                                else {$index2 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]};}
				if ($trName{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} ne "") {$tr1 = "TRRMV";}
				else {$tr1 = "-";}
				if ($trName{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} ne "") {$tr2 = "TRRMV";}
				else {$tr2 = "-";}
                        }
			elsif ($check_inv == 0 && $check_trAdd == 0 && $check_trRmv == 0)	
			{
				if ($inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}} ne ""){$index1 = $inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}}; $tr1 = "-";}
                                else {$index1 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}; $tr1 = "-";}
                                if ($inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}} ne ""){$index2 = $inv{$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}}; $tr2 = "-";}
                                else {$index2 = $bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}; $tr2 = "-";}
			}
			else
			{
				print "Can't handle both inversion and translocation right now. Still under development\n";
				exit;
			}
			print out "$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[0]}\t$bedInfo_index_conv_id{$i}{$chr{$i}}{@matData[1]}\t$index1\t$index2\t@matData[2]\t$cnv1\t$cnv2\t$tr1\t$tr2\n";
		}
                undef @matData;
        }
        close out;
        close bed;
	$i++;
} 
