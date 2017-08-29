#Author abhijit
#Created on April 15, 2017
#This script will handle the primary configuration file to simulate/generate both intra and inter-chromosomal Hi-C matrices
#

sub missing_file{
	my $fileMissing = $_[0];
	print "Didn't found $fileMissing! Exiting.\n";
	exit;
}

my $filename = @ARGV[0];
open(fh,$filename);
@file = <fh>;
$i = 0;
while ($i <= $#file){
	chomp $file[$i];
	if ($file[$i] eq "#Input Files"){
		if ($file[$i+1] =~ /^bedFile/){$bedFileName = $file[$i+1];} else{missing_file("bedFile")}
		if ($file[$i+2] =~ /^matFile/){$matFileName = $file[$i+2];} else{missing_file("matFile")}
		@bedName = split(/\s+/,$bedFileName);
		@matName = split(/\s+/,$matFileName);
		$bedFileName = @bedName[1];
		$matFileName = @matName[1];
		undef @bedName;
		undef @matName;
	}
	if ($file[$i] eq "#Output File Prefix"){
		if ($file[$i+1] =~ /^output/){$outFileName = $file[$i+1];} else{missing_file("output File Name")}
		@outName = split(/\s+/,$outFileName);
		$outFileName = @outName[1];
		undef @outName;
	}
	if ($file[$i] =~ /^#1st Chromosome/){
		open (chrA_out,">chrA.Intra.conf");
		print chrA_out "#Simulation 1\n";
		print chrA_out "bedFile\t$outFileName.New.bed\n";
		print chrA_out "matFile\t$outFileName.New.matrix\n";
		print chrA_out "output\t$outFileName.chrA.Intra\n";
		$j = $i+1;
		do{
			chomp $file[$j];
			if ($file[$j] =~ /^chr/){
				@chrA = split(/\s+/,$file[$j]);
				$chrA_name = @chrA[1];
				undef @chrA;
			} 	
			print chrA_out "$file[$j]\n";
			$j++;
		}while ($file[$j] !~ /^#/);
		close chrA_out;
	}
	if ($file[$i] =~ /^#2nd Chromosome/){
                open (chrB_out,">chrB.Intra.conf");
                print chrB_out "#Simulation 1\n";
                print chrB_out "bedFile\t$outFileName.New.bed\n";
                print chrB_out "matFile\t$outFileName.New.matrix\n";
                print chrB_out "output\t$outFileName.chrB.Intra\n";
                $j = $i+1;
                do{
                        chomp $file[$j];
			if ($file[$j] =~ /^chr/){
                                @chrB = split(/\s+/,$file[$j]);
                                $chrB_name = @chrB[1];
                                undef @chrB;
                        }
                        print chrB_out "$file[$j]\n";
                        $j++;
                }while ($file[$j] !~ /^#/);
                close chrB_out;
        }
	if ($file[$i] =~ /^#Translocation regions/){
		open (chrAB_out,">chrAB.Inter.conf");
		print chrAB_out "#Simulation 1\n";
		print chrAB_out "chromosomes\t$chrA_name\t$chrB_name\n";
		print chrAB_out "bedFile\t$outFileName.New.bed\n";
		print chrAB_out "matFile\t$outFileName.New.matrix\n";
		print chrAB_out "chrA\t$chrA_name\t$outFileName.chrAB.Inter\t$chrA_start\t$chrA_end\n";
		print chrAB_out "chrB\t$chrB_name\t$outFileName.chrAB.Inter\t$chrB_start\t$chrB_end\n";
		$j = $i+1;
                do{
                        chomp $file[$j];
                        print chrAB_out "$file[$j]\n";
                        $j++;
                }while ($file[$j] !~ /^#/);
		close chrAB_out;
	}
	if ($file[$i] eq "#Contact count treatment"){	
		@chrA_treat = split(/\s+/,$file[$i+1]);
		@chrB_treat = split(/\s+/,$file[$i+2]);
		$chrA_count_treatment = @chrA_treat[1];
		$chrB_count_treatment = @chrB_treat[1]; 
		undef @chrA_treat;
		undef @chrB_treat;
	}
	if ($file[$i] eq "#Contact counts"){
		@intra = split(/\s+/,$file[$i+1]);
		@inter = split(/\s+/,$file[$i+2]);
		if (@intra[1] eq ""){
			print "Please prodive the mode (scaled/random/self/self_random) to generate intra contact counts! Exiting\n";
			exit;
		}
		if (@inter[1] eq ""){
                        print "Please prode the mode (scaled/random/self/self_random) to generate inter contact counts! Exiting\n";
                        exit;
                }
		$intra_count = @intra[1];
		$inter_count = @inter[1];
		undef @intra;
		undef @inter;
		
	}
	if ($file[$i] eq "#Chromosome sizes"){
		@chrA = split(/\s+/,$file[$i+1]);
		@chrB = split(/\s+/,$file[$i+2]);
		$chrA_start = @chrA[1];
		$chrA_end   = @chrA[2];
		$chrB_start = @chrB[1];
                $chrB_end   = @chrB[2];
		undef @chrA;
		undef @chrB;
	}	
	$i++;
}


#print "$bedFileName\n";
#print "$matFileName\n";
#print "$chrA_name\t$chrA_start\t$chrA_end\n";
#print "$chrB_name\t$chrB_start\t$chrB_end\n";
#print "$outFileName\n";
#print "$translocationType\n";
#print "$balanceType\n";

$c = 0;
open (bed,"$bedFileName");
while (<bed>){
	chomp $_;
	@bedData = split(/\s+/,$_);
	if (@bedData[0] eq $chrA_name){
		$c++;
		$bedIndexChr{@bedData[3]} = @bedData[0];
		$bedNewIndex{@bedData[0]} += 1;	
		$bedConvertedIndex{@bedData[3]} = $bedNewIndex{@bedData[0]};
		push @newBed,"@bedData[0]\t@bedData[1]\t@bedData[2]\t$bedNewIndex{@bedData[0]}\n";
	}
	undef @bedData;
}
close bed;

$c = 0;
$bedNewIndex{$chrB_name} = $bedNewIndex{$chrA_name};
open (bed,"$bedFileName");
while (<bed>){
        chomp $_;
        @bedData = split(/\s+/,$_);
        if (@bedData[0] eq $chrB_name){
		$c++;
		$bedIndexChr{@bedData[3]} = @bedData[0];
                $bedNewIndex{@bedData[0]} += 1;
                $bedConvertedIndex{@bedData[3]} = $bedNewIndex{@bedData[0]};
                push @newBed,"@bedData[0]\t@bedData[1]\t@bedData[2]\t$bedNewIndex{@bedData[0]}\n";
        }
        undef @bedData;
}
close bed;

open (bed_out,">$outFileName.New.bed");
print bed_out @newBed;
close bed_out;

open (mat,"$matFileName");
open (mat_out,">$outFileName.New.matrix");
while (<mat>){
	chomp $_;
	@matData = split(/\s+/,$_);
	if ($bedIndexChr{@matData[0]} eq $chrA_name && $bedIndexChr{@matData[1]} eq $chrA_name){
		print mat_out "$bedConvertedIndex{@matData[0]}\t$bedConvertedIndex{@matData[1]}\t@matData[2]\n";
	}
	if ($bedIndexChr{@matData[0]} eq $chrB_name && $bedIndexChr{@matData[1]} eq $chrB_name){
                print mat_out "$bedConvertedIndex{@matData[0]}\t$bedConvertedIndex{@matData[1]}\t@matData[2]\n";
        }
	if (($bedIndexChr{@matData[0]} eq $chrA_name && $bedIndexChr{@matData[1]} eq $chrB_name) || ($bedIndexChr{@matData[0]} eq $chrB_name && $bedIndexChr{@matData[1]} eq $chrA_name)){
		print mat_out "$bedConvertedIndex{@matData[0]}\t$bedConvertedIndex{@matData[1]}\t@matData[2]\n";
	}
	undef @matData;
}
close mat;
close mat_out;

`perl generateFrameWork.Intra.pl chrA.Intra.conf`;
`perl generateFrameWork.Intra.pl chrB.Intra.conf`;
`perl swapping.transFrameWork.pl chrAB.Inter.conf`;
`perl transformPreSimToSimMatrix.pl $chrA_name.modified.index $outFileName.chrA.Intra.$chrA_name.PreSimulation.matrix $chrA_count_treatment $chrB_name.modified.index $outFileName.chrB.Intra.$chrB_name.PreSimulation.matrix $chrB_count_treatment $outFileName.chrA.Intra.$chrA_name.Simulation.matrix $outFileName.chrB.Intra.$chrB_name.Simulation.matrix`;
`perl generateInterMatrix.pl $outFileName.chrAB.Inter.$chrA_name-$chrB_name.Simulation.matrix`;
`perl generateIntraMatrix.pl $outFileName.chrA.Intra.$chrA_name.Simulation.matrix > $chrA_name.Intra.sparse.matrix`;
`perl generateIntraMatrix.pl $outFileName.chrB.Intra.$chrB_name.Simulation.matrix > $chrB_name.Intra.sparse.matrix`;
`rm *.Intra-Inter.sparse.matrix`;
`cat *.sparse.matrix > $chrA_name.$chrB_name.Intra-Inter.sparse.matrix`;
`perl createSparseToMatrix.pl $outFileName.New.bed $chrA_name.$chrB_name.Intra-Inter.sparse.matrix $outFileName.$chrA_name.$chrB_name.Figure.axis $outFileName.$chrA_name.$chrB_name.Figure.matrix`;
`Rscript figure.translocation.r $outFileName.$chrA_name.$chrB_name.Figure.matrix $outFileName.$chrA_name.$chrB_name.Figure.axis $outFileName.$chrA_name.$chrB_name.Figure.png`;
`mv $outFileName.New.bed $outFileName.bed`;
`mv $outFileName.New.matrix $outFileName.original.matrix`;
`rm chrA.Intra.conf chrB.Intra.conf chrAB.Inter.conf`;
`rm *.modified.index`;
`rm *.PreSimulation.matrix`;
`rm *.Simulation.matrix`;
`rm *.Intra.sparse.matrix`;
`rm *.Figure.axis`;
`rm *.Figure.matrix`;
`rm InterChromsomal.sparse.matrix`;
`rm chrAB.Inter.conf.*`;
`rm valid.counts.r`;
`rm Mean.DistanceWise.contact`;
`rm chr.bed chr.mod.bed`;
`rm *.*.Intra.*.bed`; 
