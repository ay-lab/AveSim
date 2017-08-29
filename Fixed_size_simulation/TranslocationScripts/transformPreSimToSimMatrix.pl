#author abhijit
#created on April 16, 2017
#This script will transform the PreSimulation intr matrices to Simulation matrices. This is do the neccessary insertion and deletion from the Intra Hi-C matrix
#
#@chrA_index  = `cat chr19.modified.index`; 
#@chrA_preSim = `cat IMR90.chrA.Intra.chr19.PreSimulation.matrix`;
#$chrA_contact_type = "HTD"; 
#@chrB_index  = `cat chr20.modified.index`;
#@chrB_preSim = `cat IMR90.chrB.Intra.chr20.PreSimulation.matrix`;
#$chrB_contact_type = "HTD";

@chrA_index  = `cat @ARGV[0]`;
@chrA_preSim = `cat @ARGV[1]`;
$chrA_contact_type = @ARGV[2];
@chrB_index  = `cat @ARGV[3]`;
@chrB_preSim = `cat @ARGV[4]`;
$chrB_contact_type = @ARGV[5];

$i = 0;
while ($i <= $#chrA_index){
	chomp $chrA_index[$i];
	@indexData = split(/\s+/,$chrA_index[$i]);
	$chromosomeInfo{@indexData[1]} = $chrA_index[$i];
	undef @indexData;
	$i++;
}
$i = 0;
while ($i <= $#chrB_index){
        chomp $chrB_index[$i];
        @indexData = split(/\s+/,$chrB_index[$i]);
        $chromosomeInfo{@indexData[1]} = $chrB_index[$i];
        undef @indexData;
        $i++;
}

open (preSim_outA, ">@ARGV[6]");
$i = 0;
while ($i <= $#chrA_preSim){
	chomp $chrA_preSim[$i];
	@preSimData = split(/\s+/,$chrA_preSim[$i]);
	@dataA = split(/\s+/,$chromosomeInfo{@preSimData[0]});
	@dataB = split(/\s+/,$chromosomeInfo{@preSimData[1]});
	if ($chromosomeInfo{@preSimData[0]} ne "" && $chromosomeInfo{@preSimData[1]} ne "" ){
		if ($chrA_contact_type eq "HTD"){
			if ((@dataA[3] eq "+" && @dataB[3] eq "+") || (@dataA[3] eq "-" && @dataB[3] eq "-")){
				$diffA = abs(@dataA[1]-@dataB[1]);
				$diffB = abs(@dataA[2]-@dataB[2]);
				if ($diffA == $diffB){
					print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
				else {
					print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
			}
			else {
					print preSim_outA "@preSimData[0]\t@preSimData[1]\t1\t1000\t@preSimData[4]\tNRM\tHTD\t-\t-\n";
			}
		}
		if ($chrA_contact_type eq "HMD"){
	                if ((@dataA[3] eq "+" && @dataB[3] eq "+") || (@dataA[3] eq "-" && @dataB[3] eq "-")){
				$diffA = abs(@dataA[1]-@dataB[1]);
                                $diffB = abs(@dataA[2]-@dataB[2]);
                                if ($diffA == $diffB){
					print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
				else {
		                        print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
	                }
	        }
		if ($chrA_contact_type eq "NRM"){
			$diffA = abs(@dataA[1]-@dataB[1]);
                      	$diffB = abs(@dataA[2]-@dataB[2]);
			if ($diffA == $diffB){
				print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
			}
			else {
				print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
			}
		}
	}
	$difference = abs(@preSimData[0]-@preSimData[1]);
	$distanceWiseCount{$difference} += @preSimData[4];
	$frequency{$difference} += 1;
	if ($check_diff{$difference} eq ""){
		push @diffList,"$difference\n";
		$check_diff{$difference} = "yes";
	}
	undef @dataA;
	undef @dataB;
	undef @preSimData;
	$i++;
}
close preSim_outA;

open (preSim_outB,">@ARGV[7]");
$i = 0;
while ($i <= $#chrB_preSim){
        chomp $chrB_preSim[$i];
        @preSimData = split(/\s+/,$chrB_preSim[$i]);
        @dataA = split(/\s+/,$chromosomeInfo{@preSimData[0]});
        @dataB = split(/\s+/,$chromosomeInfo{@preSimData[1]});
	if ($chromosomeInfo{@preSimData[0]} ne "" && $chromosomeInfo{@preSimData[1]} ne "" ){
	        if ($chrB_contact_type eq "HTD"){
	                if ((@dataA[3] eq "+" && @dataB[3] eq "+") || (@dataA[3] eq "-" && @dataB[3] eq "-")){
				$diffA = abs(@dataA[1]-@dataB[1]);
	                        $diffB = abs(@dataA[2]-@dataB[2]);
        	                if ($diffA == $diffB){
	        	                print preSim_outB "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
				else {
					print preSim_outB "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
	                }
	                else {
	                        print preSim_outB "@preSimData[0]\t@preSimData[1]\t1\t1000\t@preSimData[4]\tNRM\tHTD\t-\t-\n";
	                }
	        }
	        if ($chrB_contact_type eq "HMD"){
	                if ((@dataA[3] eq "+" && @dataB[3] eq "+") || (@dataA[3] eq "-" && @dataB[3] eq "-")){
				$diffA = abs(@dataA[1]-@dataB[1]);
                                $diffB = abs(@dataA[2]-@dataB[2]);
                                if ($diffA == $diffB){
		                        print preSim_outB "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
				else {
					print preSim_outB "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
	                }
	        }
		if ($chrB_contact_type eq "NRM"){
			if ((@dataA[3] eq "+" && @dataB[3] eq "+") || (@dataA[3] eq "-" && @dataB[3] eq "-")){
				$diffA = abs(@dataA[1]-@dataB[1]);
        	          	$diffB = abs(@dataA[2]-@dataB[2]);
                	      	if ($diffA == $diffB){
			                print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
				else {
					print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[2]\t@dataB[2]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
				}
			}
			else {
				print preSim_outA "@preSimData[0]\t@preSimData[1]\t@dataA[1]\t@dataB[1]\t@preSimData[4]\tNRM\tNRM\t-\t-\n";
			}
	        }
	}
	$difference = abs(@preSimData[0]-@preSimData[1]);
        $distanceWiseCount{$difference} += @preSimData[4];
        $frequency{$difference} += 1;
	if ($check_diff{$difference} eq ""){
                push @diffList,"$difference\n";
                $check_diff{$difference} = "yes";
        }
        undef @dataA;
        undef @dataB;
        undef @preSimData;
        $i++;
}
close preSim_outB;

open (mean_out,">Mean.DistanceWise.contact");
@sorted_diffList = sort{$a<=>$b}@diffList;
$i = 0;
while ($i <= $#sorted_diffList){
	chomp $sorted_diffList[$i];
	$mean = int($distanceWiseCount{$sorted_diffList[$i]}/$frequency{$sorted_diffList[$i]});
	print mean_out "$sorted_diffList[$i]\t$mean\t$frequency{$sorted_diffList[$i]}\n";
	$i++;
}
close mean_out;

undef @chrA_index;
undef @chrB_index;
undef @chrA_preSim;
undef @chrB_preSim;
undef @diffList;
undef @sorted_diffList;
undef %chromosomeInfo;
undef %distanceWiseCount;
undef %frequency;
undef %check_diff;
