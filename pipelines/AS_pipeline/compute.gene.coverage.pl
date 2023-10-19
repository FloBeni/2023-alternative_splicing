use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $exons=$_[1];
    my $geneexons=$_[2];
    
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $geneid=$s[0];
	my $exonid=$s[1];
	my $chr=$s[2];
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];
	
	if(exists $geneexons->{$geneid}){
	    push(@{$geneexons->{$geneid}}, $exonid);
	} else{
	    $geneexons->{$geneid}=[$exonid];
	}
	
	$exons->{$exonid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand};
	    	
	$line=<$input>;
    }
    
    close($input);
}


##############################################################

sub orderExons{
    my $exons=$_[0];
    my $refordered=$_[1];

    my %hashstart;
    
    foreach my $exid (keys %{$exons}){
	my $chr=$exons->{$exid}{"chr"};
	my $b=$exons->{$exid}{"start"};
	my $e=$exons->{$exid}{"end"};
	my $s=$exons->{$exid}{"strand"};
	
	if(exists $hashstart{$chr}){
	    if(exists $hashstart{$chr}{$b}){
		push(@{$hashstart{$chr}{$b}{"end"}},$e);
		push(@{$hashstart{$chr}{$b}{"strand"}},$s);
		push(@{$hashstart{$chr}{$b}{"id"}},$exid);
	    }
	    else{
		$hashstart{$chr}{$b}={"end"=>[$e],"strand"=>[$s],"id"=>[$exid]};
	    }
	}
	else{
	    $hashstart{$chr}={$b=>{"end"=>[$e],"strand"=>[$s],"id"=>[$exid]}};
	}
    }
    
    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={};
	
	my @uniquestart=keys %{$hashstart{$chr}};
	
	my @sortedstart = sort {$a <=> $b} @uniquestart;
	
	foreach my $b (@sortedstart){
	    
	    my $nbblocks=@{$hashstart{$chr}{$b}{"end"}};
	    
	    for(my $i=0;$i<$nbblocks;$i++){
		my $strand=${$hashstart{$chr}{$b}{"strand"}}[$i];
		
		if(!(exists $refordered->{$chr}{$strand})){
		    $refordered->{$chr}{$strand}={"start"=>[], "end"=>[], "id"=>[]};
		}
		
		push(@{$refordered->{$chr}{$strand}{"start"}},$b);
		push(@{$refordered->{$chr}{$strand}{"end"}},${$hashstart{$chr}{$b}{"end"}}[$i]);
		push(@{$refordered->{$chr}{$strand}{"id"}},${$hashstart{$chr}{$b}{"id"}}[$i]);
	    }
	}
    }
}

##############################################################

sub readCoverage{
    my $pathin=$_[0];
    my $coverage=$_[1];

    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    }
    else{
	open($input, $pathin);
    }
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t" , $line);
	
	my $chr=$s[0];
	my $start=$s[1]+1; ## originally 0-based, included, now 1-based, included
	my $end=$s[2]; ## originally 0-based, not included, now 1-based, included
	my $score=$s[3]+0.0;
	
	if(exists $coverage->{$chr}){
	    my $laststart=${$coverage->{$chr}{"start"}}[-1];
	    
	    if($start<=$laststart){
		print "Data are not ordered! ".$start." ".$laststart."\n";
		exit(1);
	    }
	    
	    push(@{$coverage->{$chr}{"start"}}, $start);
	    push(@{$coverage->{$chr}{"end"}}, $end);
	    push(@{$coverage->{$chr}{"score"}}, $score);
	}
	else{
	    $coverage->{$chr}={"start"=>[$start], "end"=>[$end], "score"=>[$score]};
	}
	
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub computeExonCoverage{
    my $orderedexons=$_[0];
    my $coverage=$_[1];
    my $covexons=$_[2];

    foreach my $chr (keys %{$orderedexons}){
	foreach my $strand (keys %{$orderedexons->{$chr}}){
	    print $chr."\t".$strand."\n";

	    if(exists $coverage->{$chr}){
		my $nbex=@{$orderedexons->{$chr}{$strand}{"start"}};
		my $nbreg=@{$coverage->{$chr}{"start"}};
		
		my $firstreg=0;
		
		for(my $i=0; $i<$nbex; $i++){
		    my $startex=${$orderedexons->{$chr}{$strand}{"start"}}[$i];
		    my $endex=${$orderedexons->{$chr}{$strand}{"end"}}[$i];
		    my $idex=${$orderedexons->{$chr}{$strand}{"id"}}[$i];
		    
		    my $lenex=($endex-$startex+1.0);
		    my $sumcov=0;
		    
		    my $j=$firstreg;
		    
		    while($j<$nbreg && ${$coverage->{$chr}{"end"}}[$j] < $startex){
			$j++;
		    }

		    $firstreg=$j;
		    
		    while($j<$nbreg && ${$coverage->{$chr}{"start"}}[$j] <= $endex){
			my $startreg=${$coverage->{$chr}{"start"}}[$j];
			my $endreg=${$coverage->{$chr}{"end"}}[$j];
			
			my $M=max($startex, $startreg);
			my $m=min($endex, $endreg);
			
			my $lenov=($m-$M+1);
			
			if($lenov>=1){
			    my $score=${$coverage->{$chr}{"score"}}[$j];
			    
			    $sumcov+=($lenov+0.0)*($score+0.0);
			}
			
			$j++;
		    }

		    my $meancov=($sumcov+0.0)/($lenex);
		    
		    $covexons->{$idex}=$meancov;
		}
	    }
	    
	}
	
    }
}

##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes exon and gene RNA-seq coverage.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }

    print "\n";
}

##############################################################
##############################################################

my %parameters;
$parameters{"pathExonBlocks"}="NA";
$parameters{"pathCoverage"}="NA";
$parameters{"pathOutputExons"}="NA";
$parameters{"pathOutputGenes"}="NA";

my @defaultpars=("pathExonBlocks", "pathCoverage", "pathOutputExons", "pathOutputGenes");

my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}

## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##############################################################
##############################################################

print "Reading annotations...\n";

my %exonblocks;
my %geneexons;

readExonBlocks($parameters{"pathExonBlocks"}, \%exonblocks, \%geneexons);

my $nbex=keys %exonblocks;
my $nbg=keys %geneexons;

print "Found ".$nbex." exon blocks and ".$nbg." genes.\n";

print "Done.\n";

print "Ordering exons...\n";

my %orderedexons;

orderExons(\%exonblocks, \%orderedexons);

print "Done.\n";

##############################################################

print "Reading RNA-seq coverage...\n";

my %coverage;
readCoverage($parameters{"pathCoverage"}, \%coverage);

print "Done.\n";

##############################################################

print "Computing exon coverage...\n";

my %exoncoverage;
computeExonCoverage(\%orderedexons, \%coverage, \%exoncoverage);

print "Done.\n";

##############################################################

print "Writing output for exons...\n";

open(my $outputex, ">".$parameters{"pathOutputExons"});
open(my $outputgenes, ">".$parameters{"pathOutputGenes"});

print $outputex "ExonID\tChr\tStart\tEnd\tStrand\tGeneID\tLength\tCoverage\n";
print $outputgenes "GeneID\tTotalExonicLength\tCoverage\n";

foreach my $geneid (keys %geneexons){
    my $totlen=0; 
    my $sumcov=0;
    
    foreach my $exonid (@{$geneexons{$geneid}}){
	my $chr=$exonblocks{$exonid}{"chr"};
	my $start=$exonblocks{$exonid}{"start"};
	my $end=$exonblocks{$exonid}{"end"};
	my $strand=$exonblocks{$exonid}{"strand"};
	
	my $exlen=$end-$start+1;
	my $cov=$exoncoverage{$exonid};

	print $outputex $exonid."\t".$chr."\t".$start."\t".$end."\t".$strand."\t".$geneid."\t".$exlen."\t".$cov."\n";

	$totlen+=$exlen;
	$sumcov+=($cov*$exlen);
    }

    my $meancov=$sumcov/($totlen+0.0);

    print $outputgenes $geneid."\t".$totlen."\t".$meancov."\n";
}

close($outputex);
close($outputgenes);

print "Done.\n";

##############################################################
