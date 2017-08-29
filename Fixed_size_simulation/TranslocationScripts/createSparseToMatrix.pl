#author abhijit
#created on April 16, 2017
#This script will will convert the sparse format to matrix  format for figure generation
#
@bed = `cat @ARGV[0]`;
open (bed_out,">@ARGV[2]");
$i = 0;
while ($i <= $#bed){
	chomp $bed[$i];
	@bedData = split(/\s+/,$bed[$i]);
	$bedChrom{@bedData[3]} = @bedData[0];
	$bedMid{@bedData[3]} = @bedData[1]+2e4;	
	#print bed_out "$bedMid{@bedData[3]}\n";
	print bed_out "@bedData[3]\n";
	push @index,"@bedData[3]\n";
	undef @bedData;
	$i++;
}
close bed_out;

open (in,"@ARGV[1]");
while (<in>){	
	chomp $_;
	@matData = split(/\s+/,$_);
	$count{@matData[0]}{@matData[1]} = @matData[2];
	$count{@matData[1]}{@matData[0]} = @matData[2];
	undef @matData;
}
close in;

open (mat_out,">@ARGV[3]");
$i = 0;
while ($i <= $#index){
	chomp $index[$i];
	$j = 0;
	while ($j <= $#index){
		chomp $index[$j];
		if ($count{$index[$i]}{$index[$j]} eq ""){$c = 0;}
		else {$c = $count{$index[$i]}{$index[$j]};}
		$seq .= "$c\t";
		$j++;
	}
	chop $seq;
	print mat_out "$seq\n";
	$seq = undef;
	$i++;
}
close mat_out;
