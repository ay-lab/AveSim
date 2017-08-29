#author abhijit
#created on April 16, 2017
#This script will generate the intrachromosomal Hi-C matrices adjusted with insertion and deletion events
#

@comb = qw(
NRM_NRM
NRM_AMP
NRM_HMD
NRM_HTD
AMP_AMP
AMP_HMD
AMP_HTD
HMD_HMD
HTD_HMD
HTD_HTD
);

$i = 0;
while ($i <= $#comb)
{
	chomp $comb[$i];
	@cnvData = split(/_/,$comb[$i]);
	@file = `cat ../data/ENCODE.CancerCell.matrix.$comb[$i].param`;
	$j = 0;
	while ($j <= $#file)
	{
		chomp $file[$j];
		@fileData = split(/\s+/,$file[$j]);
		
		$dist     = @fileData[1];
		$nbTheta  = @fileData[2];
		$nbMu     = @fileData[3];
		$nbBIC    = @fileData[4];
		$psMu     = @fileData[5];
		$psBIC    = @fileData[6];
		$n        = @fileData[7];
		
		$paramInfo_nbTheta{$dist}{@cnvData[0]}{@cnvData[1]} = $nbTheta;
		$paramInfo_nbTheta{$dist}{@cnvData[1]}{@cnvData[0]} = $nbTheta;
		$paramInfo_nbMu{$dist}{@cnvData[0]}{@cnvData[1]}    = $nbMu;
		$paramInfo_nbMu{$dist}{@cnvData[1]}{@cnvData[0]}    = $nbMu;
		$paramInfo_nbBIC{$dist}{@cnvData[0]}{@cnvData[1]}   = $nbBIC;
		$paramInfo_nbBIC{$dist}{@cnvData[1]}{@cnvData[0]}   = $nbBIC;
		$paramInfo_psMu{$dist}{@cnvData[0]}{@cnvData[1]}    = $psMu;
		$paramInfo_psMu{$dist}{@cnvData[1]}{@cnvData[0]}    = $psMu;
		$paramInfo_psBIC{$dist}{@cnvData[0]}{@cnvData[1]}   = $psBIC;
		$paramInfo_psBIC{$dist}{@cnvData[1]}{@cnvData[0]}   = $psBIC;
		$paramInfo_pop{$dist}{@cnvData[0]}{@cnvData[1]}     = $n;
		$paramInfo_pop{$dist}{@cnvData[1]}{@cnvData[0]}     = $n;
				
		if (@cnvData[0] eq "HMD" || @cnvData[1] eq "HMD")
                {
                        $paramInfo_nbTheta{$dist}{@cnvData[0]}{@cnvData[1]} = 0;
                        $paramInfo_nbTheta{$dist}{@cnvData[1]}{@cnvData[0]} = 0;
                        $paramInfo_nbMu{$dist}{@cnvData[0]}{@cnvData[1]}    = 0;
                        $paramInfo_nbMu{$dist}{@cnvData[1]}{@cnvData[0]}    = 0;
                        $paramInfo_psMu{$dist}{@cnvData[0]}{@cnvData[1]}    = 0;
                        $paramInfo_psMu{$dist}{@cnvData[1]}{@cnvData[0]}    = 0;
                }

		if ($j == $#file)
		{
			$k = $dist;
			while ($k <= 7000)
			{
				$paramInfo_nbTheta{$k}{@cnvData[0]}{@cnvData[1]} = $nbTheta;
		                $paramInfo_nbTheta{$k}{@cnvData[1]}{@cnvData[0]} = $nbTheta;
        		        $paramInfo_nbMu{$k}{@cnvData[0]}{@cnvData[1]}    = $nbMu;
                		$paramInfo_nbMu{$k}{@cnvData[1]}{@cnvData[0]}    = $nbMu;
                		$paramInfo_nbBIC{$k}{@cnvData[0]}{@cnvData[1]}   = $nbBIC;
                		$paramInfo_nbBIC{$k}{@cnvData[1]}{@cnvData[0]}   = $nbBIC;
                		$paramInfo_psMu{$k}{@cnvData[0]}{@cnvData[1]}    = $psMu;
                		$paramInfo_psMu{$k}{@cnvData[1]}{@cnvData[0]}    = $psMu;
                		$paramInfo_psBIC{$k}{@cnvData[0]}{@cnvData[1]}   = $psBIC;
                		$paramInfo_psBIC{$k}{@cnvData[1]}{@cnvData[0]}   = $psBIC;
                		$paramInfo_pop{$k}{@cnvData[0]}{@cnvData[1]}     = $n;
                		$paramInfo_pop{$k}{@cnvData[1]}{@cnvData[0]}     = $n;
				if (@cnvData[0] eq "HMD" || @cnvData[1] eq "HMD")
                		{
                        		$paramInfo_nbTheta{$k}{@cnvData[0]}{@cnvData[1]} = 0;
                        		$paramInfo_nbTheta{$k}{@cnvData[1]}{@cnvData[0]} = 0;
                        		$paramInfo_nbMu{$k}{@cnvData[0]}{@cnvData[1]}    = 0;
                        		$paramInfo_nbMu{$k}{@cnvData[1]}{@cnvData[0]}    = 0;
                        		$paramInfo_psMu{$k}{@cnvData[0]}{@cnvData[1]}    = 0;
                        		$paramInfo_psMu{$k}{@cnvData[1]}{@cnvData[0]}    = 0;
                		}
				$k++;
			}
		}
		undef @fileData;
		$j++;
	}
	undef @file; 
	undef @cnvData;
	$i++;
}

open (mat,"@ARGV[0]");
while (<mat>)
{
	chomp $_;
	@matData  = split(/\s+/,$_);
	$distance_original = abs(@matData[0]-@matData[1]);
	$distance = abs(@matData[2]-@matData[3]);
	$obs_count= @matData[4];
	$cnv1     = @matData[5];
	$cnv2     = @matData[6];
	$trRmv1   = @matData[7];
	$trRmv2   = @matData[8];
	$function     = "mean";
	$sample       = 1;
	if ($paramInfo_nbBIC{$distance}{"NRM"}{"NRM"} < $paramInfo_psBIC{$distance}{"NRM"}{"NRM"})
        {
		if ($paramInfo_nbMu{$distance}{"NRM"}{"NRM"} == 0) {$nrm = 1;}
		else {$nrm = $paramInfo_nbMu{$distance}{"NRM"}{"NRM"};}	
		if ($paramInfo_nbMu{$distance_original}{"NRM"}{"NRM"} == 0) {$nrm_ori = 1;}
                else {$nrm_ori = $paramInfo_nbMu{$distance_original}{"NRM"}{"NRM"};}
	}	
	else
	{
		if ($paramInfo_psMu{$distance}{"NRM"}{"NRM"} == 0) {$nrm = 1;}
                else {$nrm = $paramInfo_psMu{$distance}{"NRM"}{"NRM"};}
		if ($paramInfo_psMu{$distance_original}{"NRM"}{"NRM"} == 0) {$nrm_ori = 1;}
                else {$nrm_ori = $paramInfo_psMu{$distance_original}{"NRM"}{"NRM"};}
	}
	if ($paramInfo_nbBIC{$distance}{$cnv1}{$cnv2} < $paramInfo_psBIC{$distance}{$cnv1}{$cnv2})
        {
		if ($paramInfo_nbMu{$distance}{$cnv1}{$cnv2} == 0) {$cnv = 1;}
		else {$cnv = $paramInfo_nbMu{$distance}{$cnv1}{$cnv2};}
		$ratio_to_scale_down_count = $nrm/$nrm_ori;
		$current_obs_count = $obs_count*$ratio_to_scale_down_count;
		$ratio = $cnv/$nrm;
             	$random_count_cnv = sprintf("%0.0f",$current_obs_count * $ratio);
		print "@matData[0]\t@matData[1]\t$random_count_cnv\n";
        }
	else
	{
                if ($paramInfo_psMu{$distance}{$cnv1}{$cnv2} == 0) {$cnv = 1;}
		else {$cnv = $paramInfo_psMu{$distance}{$cnv1}{$cnv2};}
		$ratio_to_scale_down_count = $nrm/$nrm_ori;
                $current_obs_count = $obs_count*$ratio_to_scale_down_count;
		$ratio =$cnv/$nrm;
                $random_count_cnv = sprintf("%0.0f",$current_obs_count * $ratio);
                print "@matData[0]\t@matData[1]\t$random_count_cnv\n";
	}
	undef @matData; 
}
close out;
