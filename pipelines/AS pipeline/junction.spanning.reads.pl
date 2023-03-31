#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readIntrons{
    my $path=$_[0];
    my $refintrons=$_[1];

    open(my $input,$path);
    
    my $line=<$input>; # header
    $line=<$input>;

    my $nbnogene=0;
    
    while($line){
	chomp $line;
	
	my @info=split("\t",$line);

	my $chr=$info[0];
	my $start=$info[1];
	my $end=$info[2];
	my $strand=$info[3];
	my $gene=$info[5];
    
	if($gene ne "NA"){
	    if(exists $refintrons->{$gene}){
		push(@{$refintrons->{$gene}{"start"}},$start);
		push(@{$refintrons->{$gene}{"end"}},$end);
	    }
	    else{
		$refintrons->{$gene}={"chr"=>$chr,"strand"=>$strand,"start"=>[$start],"end"=>[$end]};
	    }
	} else{
	    $nbnogene++;
	}
	
	$line=<$input>;
    }

    close($input);

    print "Removed ".$nbnogene." introns that were not assigned to annotated genes.\n";
}

##############################################################

sub extractIntronBorders{
    my $refintrons=$_[0];
    my $halfsize=$_[1];
    my $refborders=$_[2];

    foreach my $gene (keys %{$refintrons}){
	
	my $nbintrons=@{$refintrons->{$gene}{"start"}};
	my $chr=$refintrons->{$gene}{"chr"};
	my $strand=$refintrons->{$gene}{"strand"};
	
	for(my $i=0;$i<$nbintrons;$i++){
	    
	    my $start=${$refintrons->{$gene}{"start"}}[$i];
	    my $end=${$refintrons->{$gene}{"end"}}[$i];
	    
	    my $key=$gene.",".$chr.",".$strand.",".$start.",".$end;

	    my $leftstart=$start-$halfsize;
	    my $leftend=$start+$halfsize-1;
	    my $rightstart=$end-$halfsize+1;
	    my $rightend=$end+$halfsize;

	    if(exists $refborders->{$key}){
		print "error !! i've seen this border (".$key.") before\n";
	    } else{
		if($strand eq "1" || $strand eq "+"){
		    $refborders->{$key}={"start"=>[$leftstart,$rightstart], "end"=>[$leftend,$rightend], "type"=>["5prime","3prime"]};
		} else{
		    if($strand eq "-1" || $strand eq "-"){
			$refborders->{$key}={"start"=>[$leftstart,$rightstart], "end"=>[$leftend,$rightend], "type"=>["3prime","5prime"]};
		    }
		    else{
			print "unknown strand: ".$strand."\n";
			exit(1);
		    }
		}
	    }
	}
	
    }
    
}

##############################################################

sub orderBorders{

    my $refunordered=$_[0];
    my $refordered=$_[1];
         
    my %hashstart;

    foreach my $key (keys %{$refunordered}){
	my @s=split(",",$key);
	my $chr=$s[1];
	
	for(my $i=0;$i<2;$i++){
	   
	    my $b=${$refunordered->{$key}{"start"}}[$i];
	    my $e=${$refunordered->{$key}{"end"}}[$i];
	    my $t=${$refunordered->{$key}{"type"}}[$i];

	    if(exists $hashstart{$chr}){
		if(exists $hashstart{$chr}{$b}){
		    push(@{$hashstart{$chr}{$b}{"end"}},$e);
		    push(@{$hashstart{$chr}{$b}{"type"}},$t);
		    push(@{$hashstart{$chr}{$b}{"key"}},$key);
		}
		else{
		    $hashstart{$chr}{$b}={"end"=>[$e],"type"=>[$t],"key"=>[$key]};
		}
	    }
	    else{
		$hashstart{$chr}={$b=>{"end"=>[$e],"type"=>[$t],"key"=>[$key]}};
	    }
	    
	}
    }
    
    foreach my $chr (keys %hashstart){
	$refordered->{$chr}={"start"=>[],"end"=>[],"type"=>[],"key"=>[]};

	my @uniquestart=keys %{$hashstart{$chr}}; 
	my @sortedstart=sort {$a <=> $b} @uniquestart;
	

	foreach my $start (@sortedstart){
	    my $nb=@{$hashstart{$chr}{$start}{"end"}};
	    
	    for(my $i=0;$i<$nb;$i++){
		push(@{$refordered->{$chr}{"start"}},$start);
		push(@{$refordered->{$chr}{"end"}},${$hashstart{$chr}{$start}{"end"}}[$i]);
		push(@{$refordered->{$chr}{"type"}},${$hashstart{$chr}{$start}{"type"}}[$i]);
		push(@{$refordered->{$chr}{"key"}},${$hashstart{$chr}{$start}{"key"}}[$i]);
	    }
	}
    }
}


##############################################################

sub findNbHits{
    my $s=$_[0];

    my $nh="NA";

    my @grepres=grep(/NH:i/,@{$s});
    
    if(@grepres==1){
	my @t=split(":",$grepres[0]);
	$nh=$t[2]+0;
    }

    return $nh;
}

##############################################################

sub extractReads{
    my $pathaln=$_[0]; ## path for sorted alignments
    my $refborders=$_[1]; ## ordered exon-intron borders
    my $onlyunique=$_[2]; ## do we keep only unique reads
    my $refmapped=$_[3]; ## results with the number of reads for each border
      
    my @ext=split("\\.",$pathaln);
    my $nb=@ext;
    my $extension=$ext[$nb-1];
    
    my $input;

    if($extension eq "gz"){
	open($input,"zcat $pathaln |");
    }
    else{
	if($extension eq "bam"){
	    open($input,"samtools view $pathaln |");
	} else{
	    open($input,$pathaln);
	}
    }


    my $totreadsseen=0;

    my $currentchr="NA";
    my $firstindex=0; ## where we should start looking for intron borders

    my $line=<$input>;
      
    while($line){

	my $firstchar = substr $line,0,1;
	if($firstchar eq "@"){
	    $line=<$input>;
	    next;
	}
	my @s=split("\t",$line);
	my $chr=$s[2];
	
	if(!($chr eq $currentchr)){
	    $currentchr=$chr;
	    $firstindex=0; ##we restart junction index counting

	    # print "Looking at chromosome ".$chr."\n";
	}
	

	my $alnstart=$s[3]+0; ## position starts at 1
	
	## check if read is unique 
	my $okuniqueness=1;

	if($onlyunique!=0){
	    my $nbhits=findNbHits(\@s);
	    if($nbhits!=1){
		$okuniqueness=0;
	    }
	}

	## check alignment type - only keep unspliced reads
	
	my $cigar=$s[5];
	my $indexN=index($cigar,'N'); ## $indexN has to be equal to -1
		
	if($okuniqueness==1 && $indexN==-1){

	    ## get genome start and end positions
	    
	    my @cigarsplit1=split(/[a-zA-Z]/, $cigar);
	    my @cigarsplit2=grep { /\S/ } split(/[0-9]/, $cigar);
	    	 
	    my $gstart=$alnstart;
	    my $gend=$alnstart;
	    
	    my @allgenomestart;
	    my @allgenomeend;
	    
	    my $alnlen=0;
	    
	    for(my $i=0; $i<@cigarsplit1; $i++){
		if($cigarsplit2[$i] eq "N"){
		    print "Weird! we are not supposed to have introns in this read.\n";
		    print $line."\n";
		    exit(1);
		}
		
		if($cigarsplit2[$i] eq "M"){
		    $alnlen+=$cigarsplit1[$i];
		    $gend+=$cigarsplit1[$i];
		} 
		else{
		    if($cigarsplit2[$i] eq "D"){
			$gend+=$cigarsplit1[$i];
		    }
		    ## other CIGARS: S, I do not affect genome position
		    ## there are no Ns
		}
	    }
	    
	    ## don't forget the last segment
	    
	    push(@allgenomestart, $gstart);
	    push(@allgenomeend, ($gend-1));

	    ## final position for the read on the genome
	    
	    my $start=min @allgenomestart;
	    my $end=max @allgenomeend;

	    ## check if the read overlaps with intron boundaries
	    
	    if(exists $refborders->{$currentchr}){
		my $nb=@{$refborders->{$currentchr}{"start"}};
		
		my $j=$firstindex;
		
		while($j<$nb && ${$refborders->{$currentchr}{"end"}}[$j]<$start){
		    $j++;
		}
		$firstindex=$j;

		while($j<$nb && ${$refborders->{$currentchr}{"start"}}[$j]<=$end){
		    if(${$refborders->{$currentchr}{"end"}}[$j]<=$end && ${$refborders->{$currentchr}{"start"}}[$j]>=$start){
			
			## we found a read which includes the border

			my $key=${$refborders->{$currentchr}{"key"}}[$j];

			my @tt=split(",",$key);
			
			my $type=${$refborders->{$currentchr}{"type"}}[$j];
			
			if(exists $refmapped->{$key}){
			    if(exists $refmapped->{$key}{$type}){
				$refmapped->{$key}{$type}{"nbreads"}++;
			    }
			    else{
				$refmapped->{$key}{$type}={"nbreads"=>1};
			    }
			}
			else{
			    $refmapped->{$key}={$type=>{"nbreads"=>1}};
			}
		    }
		    
		    $j++;
		}
	    }
	}
	
	$line=<$input>;
	$totreadsseen++;

	if($totreadsseen%1000000==0){
	    print "Seen ".$totreadsseen." reads so far.\n";
	}
    }

    close($input);

}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script counts the number of unspliced reads that span intron-exon boundaries.\n";
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
$parameters{"pathIntronLibrary"}="NA";
$parameters{"pathAlignment"}="NA";
$parameters{"sizeOverlap"}=8;
$parameters{"onlyUnique"}=1;
$parameters{"pathOutput"}="junctions_spanning_reads.txt";
 
my @defaultpars=("pathIntronLibrary", "pathAlignment", "sizeOverlap", "onlyUnique", "pathOutput");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

my %numericpars;

$numericpars{"sizeOverlap"}=1;

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
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0;
	}
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

print "Reading introns...\n";
my %introns;

readIntrons($parameters{"pathIntronLibrary"},\%introns);

my $nbg=keys %introns;

print "Found introns for ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Extracting borders....\n";

my %unorderedBorders;

my $size=$parameters{"sizeOverlap"}+0;

print "Size of exon-intron overlap: 2 * ".$size."\n";

extractIntronBorders(\%introns,$size,\%unorderedBorders);

print "Done.\n";

print "Ordering borders...\n";

my %borders;

orderBorders(\%unorderedBorders,\%borders);

print "Done.\n";

##############################################################

print "Extracting reads that overlap exon-intron boundaries...\n";

my %mapped;

my $onlyunique=$parameters{"onlyUnique"}+0;

if($onlyunique==0){
    print "We keep all reads.\n";
} else{
    print "We keep only unique reads.\n";
}

extractReads($parameters{"pathAlignment"},\%borders, $onlyunique, \%mapped);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output,">".$parameters{"pathOutput"});
print $output "GeneID\tChr\tIntronStart\tIntronEnd\tStrand\tNbReads5Prime\tNbReads3Prime\n";

open(my $input, $parameters{"pathIntronLibrary"});
my $lineinput=<$input>; ## header
$lineinput=<$input>;

while($lineinput){
    chomp $lineinput;
    my @info=split("\t", $lineinput);
    my $chr=$info[0];
    my $start=$info[1];
    my $end=$info[2];
    my $strand=$info[3];
    my $gene=$info[5];
    
    my $line=$gene."\t".$chr."\t".$start."\t".$end."\t".$strand;

    my $key=$gene.",".$chr.",".$strand.",".$start.",".$end;

    if(exists $mapped{$key}){
	if(exists $mapped{$key}{"5prime"}){
	    $line.="\t".$mapped{$key}{"5prime"}{"nbreads"};
	} else{
	    $line.="\t0";
	}

	if(exists $mapped{$key}{"3prime"}){
	    $line.="\t".$mapped{$key}{"3prime"}{"nbreads"};
	} else{
	    $line.="\t0";
	}
    }
    else{
	$line.="\t0\t0";
    }
       
    $line.="\n";

    print $output $line;
    $lineinput=<$input>;
}

close($output);
close($input);
print "Done.\n";

##############################################################
    
