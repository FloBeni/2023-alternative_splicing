#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readIntronLibrary{
    my $pathin=$_[0];
    my $intronlib=$_[1];
    my $allids=$_[2];

    open(my $input, $pathin);

    my $line=<$input>;## header
    my %header;
    chomp $line;
    my @s=split("\t", $line);
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    my $nbkept=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	my $genelist=$s[$header{"Gene"}];
	
	if($genelist ne "NA"){
	    my $splice5="NA";
	    my $splice3="NA";
	    
	    if($strand eq "1"){
		$splice5=$start;
		$splice3=$end;
	    } else{
		if($strand eq "-1"){
		    $splice5=$end;
		    $splice3=$start;
		}
		else{
		    print "Weird strand! ".$strand."".$line."\n";
		    exit(1);
		}
	    }
	    
	    my @genes=split(",", $genelist);
	    
	    foreach my $gene (@genes){
		
		my $id=$gene.",".$splice5.",".$splice3;
		$allids->{$id}=1;

		if(exists $intronlib->{$gene}){
		    if(exists $intronlib->{$gene}{"splice5"}{$splice5}){
			$intronlib->{$gene}{"splice5"}{$splice5}{$splice3}=0;
		    } else{
			$intronlib->{$gene}{"splice5"}{$splice5}={$splice3=>0};
		    }
		    
		    if(exists $intronlib->{$gene}{"splice3"}{$splice3}){
			$intronlib->{$gene}{"splice3"}{$splice3}{$splice5}=0;
		    } else{
			$intronlib->{$gene}{"splice3"}{$splice3}={$splice5=>0};
		    }
		} else{
		    $intronlib->{$gene}={"splice5"=>{}, "splice3"=>{}, "chr"=>$chr, "strand"=>$strand};
		    $intronlib->{$gene}{"splice5"}{$splice5}={$splice3=>0};
		    $intronlib->{$gene}{"splice3"}{$splice3}={$splice5=>0};
		}
	    }
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readSpliceJunctions{
    my $pathin=$_[0];
    my $juncbygene=$_[1];
 
    open(my $input, $pathin);

    my $line=<$input>;## header
    my %header;
    chomp $line;
    my @s=split("\t", $line);
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $chr=$s[$header{"Chr"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];
	
	my $splicesite=$s[$header{"SpliceSite"}];
	my $genes=$s[$header{"Gene"}];
	my $nbuniquereads=$s[$header{"NbReads"}]+0;
	
	my $splice5="NA";
	my $splice3="NA";
	
	if($strand eq "1"){
	    $splice5=$start;
	    $splice3=$end;
	} else{
	    if($strand eq "-1"){
		$splice5=$end;
		$splice3=$start;
	    }
	    else{
		print "Weird strand! ".$strand."".$line."\n";
		exit(1);
	    }
	}
	
	if($genes ne "NA"){
	    my @genelist=split(",",$genes);
	    
	    foreach my $gene (@genelist){
		if(exists $juncbygene->{$gene}){
		    if(exists $juncbygene->{$gene}{"splice5"}{$splice5}{$splice3} && exists $juncbygene->{$gene}{"splice3"}{$splice3}{$splice5}){
			$juncbygene->{$gene}{"splice5"}{$splice5}{$splice3}=$nbuniquereads;
			
			$juncbygene->{$gene}{"splice3"}{$splice3}{$splice5}=$nbuniquereads;
		    }
		}
	    }
	}
        
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub computeAlternativeSplicingFrequency{
    my $allids=$_[0]; 
    my $juncbygene=$_[1];
    my $pathoutput=$_[2];
   
    open(my $output, ">".$pathoutput);
   
    print $output "Gene\tChr\tStrand\tSplice5\tSplice3\tNbReadsIntron\tOtherSplice5Coords\tOtherSplice3Coords\tNbReadsOtherSplice5\tNbReadsOtherSplice3\n";
   
    foreach my $idsplice (@{$allids}){
	my @s=split(",", $idsplice);
	my $gene=$s[0];
	my $splice5=$s[1]+0;
	my $splice3=$s[2]+0;
	
	my $chr=$juncbygene->{$gene}{"chr"};
	my $strand=$juncbygene->{$gene}{"strand"};
	
	
	if(!(exists $juncbygene->{$gene}{"splice5"}{$splice5}{$splice3}) || !(exists $juncbygene->{$gene}{"splice3"}{$splice3}{$splice5})){
	    print "Weird! cannot find intron ".$splice5." ".$splice3." for ".$gene."\n";
	    exit(1);
	}
	    
	my %hashothersplice5uniq;
	my %hashothersplice3uniq;
	
	my $thisuniq=$juncbygene->{$gene}{"splice5"}{$splice5}{$splice3};
	
	my $nbuniqueother3=0;
	my $nbuniqueother5=0;
	
	foreach my $othersplice3  (keys %{$juncbygene->{$gene}{"splice5"}{$splice5}}){
	    if($othersplice3 != $splice3){
		my $otheruniq=$juncbygene->{$gene}{"splice5"}{$splice5}{$othersplice3};
		
		if($otheruniq!=0){
		    $hashothersplice3uniq{$othersplice3}=1;
		}
		
		$nbuniqueother3+=$otheruniq;
	    }
	}
	
	foreach my $othersplice5  (keys %{$juncbygene->{$gene}{"splice3"}{$splice3}}){
	    if($othersplice5 != $splice5){
		my $otheruniq=$juncbygene->{$gene}{"splice3"}{$splice3}{$othersplice5};
		
		if($otheruniq!=0){
		    $hashothersplice5uniq{$othersplice5}=1;
		}
		
		$nbuniqueother5+=$otheruniq;
		
	    }
	}
	
	
	my @sorted5uniq=sort (keys %hashothersplice5uniq);
	my $keyother5uniq=join(",",@sorted5uniq);
	if($keyother5uniq eq ""){
	    $keyother5uniq="NA";
	}
	
	my @sorted3uniq=sort (keys %hashothersplice3uniq);
	my $keyother3uniq=join(",",@sorted3uniq);
	if($keyother3uniq eq ""){
	    $keyother3uniq="NA";
	}
	
	print $output $gene."\t".$chr."\t".$strand."\t".$splice5."\t".$splice3."\t".$thisuniq."\t".$keyother5uniq."\t".$keyother3uniq."\t".$nbuniqueother5."\t".$nbuniqueother3."\n";
	
    }
    
    close($output);
    
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script computes alternative splicing frequencies.\n";
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
$parameters{"pathIntronLibrary"}="NA";
$parameters{"pathSpliceJunctions"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathIntronLibrary",  "pathSpliceJunctions","pathOutput");
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

print "Reading intron library...\n";

my %intronlib;
my %allids;

readIntronLibrary($parameters{"pathIntronLibrary"}, \%intronlib, \%allids);

my $nbintrons=keys %allids;

print "We found ".$nbintrons." introns (that are assigned to genes) in the library. All analyses will be restricted to these introns.\n";

print "Done.\n";

##############################################################

print "Reading splice junctions...\n";

readSpliceJunctions($parameters{"pathSpliceJunctions"},  \%intronlib);

print "Done.\n";

##############################################################

print "Computing AS frequency and writing output...\n";

my @uniqids=keys %allids;
my @sortedids=sort @uniqids;

computeAlternativeSplicingFrequency(\@sortedids, \%intronlib, $parameters{"pathOutput"});

print "Done.\n";

##############################################################
