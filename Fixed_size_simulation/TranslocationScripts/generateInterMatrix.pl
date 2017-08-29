#author abhijit
#created on April 16, 2017
#This script will create the inter-chromosomal sparse matrix 
#
#$dist_cutoff = 100;
#$num  = `cat IMR90.chrAB.Inter.chr19-chr20.Simulation.matrix|wc -l`;
#$top  = int($num*0.15);
#@valid  = `awk '{d = \$4-\$3; if(d >= (-1*$dist_cutoff) && d <= $dist_cutoff){print}}' IMR90.chrAB.Inter.chr19-chr20.Simulation.matrix`;
#@random = `awk '{d = \$4-\$3; if(d < (-1*$dist_cutoff) || d > $dist_cutoff){print}}' IMR90.chrAB.Inter.chr19-chr20.Simulation.matrix |shuf|head -$top`;
@valid  = `cat @ARGV[0]`;

@mean = `cat Mean.DistanceWise.contact`;
$i = 1;
while ($i <= $#mean){
	chomp $mean[$i];
	@meanData = split(/\s+/,$mean[$i]);
	$meanCount{@meanData[0]} = @meanData[1]; 
	$freqCount{@meanData[0]} = @meanData[2];
	undef @meanData;
	$i++;
}

open (out,">InterChromsomal.sparse.matrix");

open (valid_file,">valid.counts.r");
$meanCount{0} = $meanCount{1};
$i = 0;
while ($i <= $#valid){
        chomp $valid[$i];
        @validData = split(/\s+/,$valid[$i]);
        $diff = (abs(@validData[2]-@validData[3]));
        if ($meanCount{$diff} eq "" ){
		
        }
        else {
                $p = $freqCount{$diff}/$freqCount{1};
		if ($p > 1){$p = 1;}
                $q = 1-($p)**3;
                print valid_file "print (paste(as.character(\"@validData[0]-@validData[1]-$meanCount{$diff}\"),(sample(c(0,1), size=1,replace=TRUE, prob=c($q,$p))),sep=\"#\"))\n";;
        }
        undef @validData;
        $i++;
}
close valid_file;
@validFinal = `Rscript valid.counts.r |grep "#1"|awk '{print \$2}'|sed 's/"//g'|sed 's/-/ /g'|sed 's/#/ /g'|awk '{print \$1\"\\t\"\$2\"\\t\"\$3}'`;

print out @validFinal;

=cut
$i = 0;
while ($i <= $#random){
        chomp $random[$i];
        @randomData = split(/\s+/,$random[$i]);
        $diff = abs(@randomData[2]-@randomData[3]);
        print out @randomData[0],"\t",@randomData[1],"\t",1,"\n";
        undef @randomData;
        $i++;
}
=cut

close out;

undef @valid;
undef @validFinal;
undef @random;
#`rm valid.counts.r`;
