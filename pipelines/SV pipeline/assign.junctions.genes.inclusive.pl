#!/usr/bin/perl

use strict;

##############################################################

sub readExons{
    my $pathin=$_[0];
    my $exons=$_[1];
    my $chrlist=$_[2];
   
    open(my $input, $pathin);
    
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon" || $type eq "CDS"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $ss=$s[6];
	    my $strand="NA";
	    
	    if($ss eq "+"){
		$strand="1";
	    } else{
		if($ss eq "-"){
		    $strand="-1";
		} else{
		    print "Weird strand!\n";
		    print $line."\n";
		    exit(1);
		} 
	    }

	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $gene=findInfo("gene_id", \@infoarray);

	    if($gene eq "NA"){
		$gene=findInfo("geneID", \@infoarray);
	    }

	    if($gene ne "NA"){
		if(!(exists $exons->{$gene})){
		    if(exists $chrlist->{$chr}){
			push(@{$chrlist->{$chr}}, $gene);
		    } else{
			$chrlist->{$chr}=[$gene];
		    }
		    
		    $exons->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[$start], "end"=>[$end]};
		    
		} else{
		    push(@{$exons->{$gene}{"start"}}, $start);
		    push(@{$exons->{$gene}{"end"}}, $end);
		}
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

############################################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];
    
    my $res="NA";
    my $nbfound=0;
    
    my @grepres=grep(/${pattern}/,@{$array});
    
    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
	$nbfound++;
    }
    else{
	foreach my $g (@grepres){
	    my @u=split(" ",$g);
	    
	    if($u[0] eq $pattern){
		$nbfound++;
		my @t=split("\"",$g);
		$res=$t[1];
	    }
	}
    }
    
    if($nbfound==1){
	return $res;
    } else{
	print "Cannot find gene_id in this line: ".join(" ", @{$array})."\n";
	
	return "NA";
    }
}

############################################################################

sub makeExonBlocks{
    my $exons=$_[0];
    my $exonblocks=$_[1];

    foreach my $gene (keys %{$exons}){
	my $chr=$exons->{$gene}{"chr"};
	my $strand=$exons->{$gene}{"strand"};
	my $nbex=@{$exons->{$gene}{"start"}};

	my %hashcoords;
	    
	for(my $i=0; $i<$nbex; $i++){
	    my $b=${$exons->{$gene}{"start"}}[$i];
	    my $e=${$exons->{$gene}{"end"}}[$i];
	    
	    if(exists $hashcoords{$b}){
		if($e>$hashcoords{$b}){
		    $hashcoords{$b}=$e;   
		}
	    }
	    else{
		$hashcoords{$b}=$e;
	    }
	}	    
	
	$exonblocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[], "end"=>[]};
	
	my @uniquestart = keys %hashcoords;
	my @sortedstart = sort {$a <=> $b} @uniquestart;

	my $nb=@sortedstart;
	
	my $currentstart=$sortedstart[0];
	my $currentend=$hashcoords{$currentstart};
	
	for(my $u=1; $u<$nb; $u++){
	    my $thisstart=$sortedstart[$u];
	    my $thisend=$hashcoords{$thisstart};
	    
	    ## cluster blocks if they overlap
	    
	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){  
		
		## we only change the end if it's larger than the current position
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    }
	    else{
		push(@{$exonblocks->{$gene}{"start"}},$currentstart);
		push(@{$exonblocks->{$gene}{"end"}},$currentend);
		
		$currentstart=$thisstart;
		$currentend=$thisend;
		
	    }
	}
	
	## don't forget the last block
	
	push(@{$exonblocks->{$gene}{"start"}},$currentstart);
	push(@{$exonblocks->{$gene}{"end"}},$currentend);

    }
}

############################################################################

sub makeHashCoordinates{
    my $exons=$_[0];
    my $genelist=$_[1];
    my $margin=$_[2];
    my $hashcoords=$_[3];

    foreach my $gene (@{$genelist}){
	my $strand=$exons->{$gene}{"strand"};
	my $nbexons=@{$exons->{$gene}{"start"}};
	
	for(my $i=0; $i<$nbexons; $i++){
	    my $start=${$exons->{$gene}{"start"}}[$i]-$margin;
	    my $end=${$exons->{$gene}{"end"}}[$i]+$margin;
	    
	    for(my $k=$start; $k<=$end; $k++){
		
		if(exists $hashcoords->{$strand}){
		    if(exists $hashcoords->{$strand}{$k}){
			$hashcoords->{$strand}{$k}{$gene}=1;
		    }
		    else{
			$hashcoords->{$strand}{$k}={$gene=>1};
		    }
		}
		else{
		    $hashcoords->{$strand}={$k=>{$gene=>1}};
		}
	    }
	}
    }
}

##############################################################

sub readSpliceJunctions{
    my $pathin=$_[0];
    my $junctions=$_[1];

    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }
    
    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    my $nbjunc=0;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $chrname="Chr";
	my $strname="Strand";
	
	if(exists $header{"Chromosome"}){
	    $chrname="Chromosome";
	}
	
	if(exists $header{"ProbableStrand"}){
	    $strname="ProbableStrand";
	}
	
	my $chr=$s[$header{$chrname}];
	my $strand=$s[$header{$strname}];
	my $start=$s[$header{"Start"}]+0; ## intron start
	my $end=$s[$header{"End"}]+0; ## intron end
	my $nbreads=$s[$header{"NbReads"}]+0;
	my $splicesite=$s[$header{"SpliceSite"}];
 
	if(exists $junctions->{$chr}){
	    if(exists $junctions->{$chr}{$strand}){
		push(@{$junctions->{$chr}{$strand}{"start"}}, $start);
		push(@{$junctions->{$chr}{$strand}{"end"}}, $end);
	    } else{
		$junctions->{$chr}{$strand}={"start"=>[$start], "end"=>[$end]};
	    }
	} else{
	    $junctions->{$chr}={$strand=>{"start"=>[$start], "end"=>[$end]}};
	}

	$nbjunc++;
	
	$line=<$input>;
    }
    
    close($input);

    print "Found ".$nbjunc." junctions.\n";

}

##############################################################

sub identifyGenes{
    my $junctions=$_[0];
    my $exons=$_[1];
    my $margin=$_[2];
    my $chrlist=$_[3];
    my $juncgene=$_[4];

    foreach my $chr (keys %{$junctions}){

	if(exists $chrlist->{$chr}){
	    my %coordinates;
	    makeHashCoordinates($exons, $chrlist->{$chr}, $margin, \%coordinates);
		
	    foreach my $strand (keys %{$junctions->{$chr}}){
		
		if(exists $coordinates{$strand}){
		    my $nbjunc=@{$junctions->{$chr}{$strand}{"start"}};
		    
		    for(my $i=0; $i<$nbjunc; $i++){
			my $start=${$junctions->{$chr}{$strand}{"start"}}[$i];
			my $end=${$junctions->{$chr}{$strand}{"end"}}[$i];
			
			my $id=$chr.",".$start.",".$end.",".$strand;
			
			my $ex1=$start-1; ## exon coordinates
			my $ex2=$end+1;
			
			if(exists $coordinates{$strand}{$ex1} && exists $coordinates{$strand}{$ex2}){
			    my @genes1 = keys %{$coordinates{$strand}{$ex1}};
			    my @genes2 = keys %{$coordinates{$strand}{$ex2}};
			    
			    my %commongenes;
			
			    foreach my $gene1 (@genes1){
				foreach my $gene2 (@genes2){
				    if($gene1 eq $gene2){
					$commongenes{$gene1}=1;
					last;
				    }
				}
			    }
			    
			    my @common=keys %commongenes;
			    my $nbcommon=@common;
			    
			    if($nbcommon!=0){
				$juncgene->{$id}={"genes"=>[], "type"=>"OK"};
				push(@{$juncgene->{$id}{"genes"}}, @common);
			    }
			    else{
				$juncgene->{$id}={"genes"=>["NA"], "type"=>"Joins2Genes"};
			    }
			}
			else{
			    if(exists $coordinates{$strand}{$ex1}){
				my @genes1 = keys %{$coordinates{$strand}{$ex1}};
				$juncgene->{$id}={"genes"=>[], "type"=>"LeftMatchOnly"};
				push(@{$juncgene->{$id}{"genes"}}, @genes1);
			    }
			    else{
				if(exists $coordinates{$strand}{$ex2}){
				    my @genes2 = keys %{$coordinates{$strand}{$ex2}};
				    $juncgene->{$id}={"genes"=>[], "type"=>"RightMatchOnly"};
				    push(@{$juncgene->{$id}{"genes"}}, @genes2);
				    
				}
			    }
			}
		    }
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
    print "This script assigns splice junctions to genes.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

## parameters 

my %parameters;
$parameters{"pathGTF"}="NA";
$parameters{"pathSpliceJunctions"}="NA";
$parameters{"margin"}=0;
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGTF","pathSpliceJunctions", "margin","pathOutput");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

my %numericpars;

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

print "Reading exon coordinates...\n";

my %exons;
my %allchromo;
readExons($parameters{"pathGTF"}, \%exons, \%allchromo);

my $nbg=keys %exons;
my $nbchr=keys %allchromo;

print "There are ".$nbg." genes on ".$nbchr." chromosomes.\n";
print "Done.\n";

##############################################################

print "Making exon blocks...\n";

my %blocks;
makeExonBlocks(\%exons, \%blocks);
   
print "Done.\n";

##############################################################

print "Reading splice junctions...\n";

my %junctions;

readSpliceJunctions($parameters{"pathSpliceJunctions"}, \%junctions);
 
print "Done.\n";

##############################################################

print "Identifying genes for junctions...\n";

my $margin=$parameters{"margin"}+0;

print "Adding a margin of ".$margin."bp on each side of exons.\n";

my %idgene;

identifyGenes(\%junctions, \%blocks, $margin, \%allchromo, \%idgene);
    
print "Done.\n";

##############################################################

print "Reading junctions and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tStrand\tSpliceSite\tNbReads\tGene\tComment\n";

my $pathin=$parameters{"pathSpliceJunctions"};

my @s=split("\\.",$pathin);
my $ext=$s[-1];

my $input;

if($ext eq "gz"){
    open($input, "zcat $pathin |");
} else{
    open($input, $pathin);
}

my $line=<$input>; ## header
chomp $line;
my @s=split("\t", $line);
my %header;

for(my $i=0; $i<@s; $i++){
    $header{$s[$i]}=$i;
}

$line=<$input>;

my $nbjunc=0;

while($line){
    chomp $line;
    my @s=split("\t", $line);
    
    my $chrname="Chr";
    my $strname="Strand";
    
    if(exists $header{"Chromosome"}){
	$chrname="Chromosome";
    }
    
    if(exists $header{"ProbableStrand"}){
	$strname="ProbableStrand";
    }
    
    my $chr=$s[$header{$chrname}];
    my $strand=$s[$header{$strname}];
    my $start=$s[$header{"Start"}]+0; ## intron start
    my $end=$s[$header{"End"}]+0; ## intron end
    my $nbreads=$s[$header{"NbReads"}]+0;
    my $splicesite=$s[$header{"SpliceSite"}];
    
    my $id=$chr.",".$start.",".$end.",".$strand;
    
    if(exists $idgene{$id}){
	print $output $chr."\t".$start."\t".$end."\t".$strand."\t".$splicesite."\t".$nbreads."\t".join(",",@{$idgene{$id}{"genes"}})."\t".$idgene{$id}{"type"}."\n";
    } else{
	print $output $chr."\t".$start."\t".$end."\t".$strand."\t".$splicesite."\t".$nbreads."\tNA\tNoMatch\n";
    }
 
    $line=<$input>;
      
}

close($input);

close($output);

print "Done.\n";

##############################################################
