#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

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

sub readFasta{

    my $path=$_[0];
    my $reffasta=$_[1];
   
    open(my $input,$path);
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];
	    
	    $reffasta->{$id}="";

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

##############################################################

sub reverseComplement{
    my $sequence=$_[0];
    
    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}

##############################################################

sub extractJunctions{
    my $pathhits=$_[0];
    my $maxmismatch=$_[1];
    my $refjunctions=$_[2];

    my $input;
    
    ## see if the file is zipped or not
    
    my @ext=split("\\.",$pathhits);
    my $nb=@ext;
    my $extension=$ext[$nb-1];
    
    if($extension eq "gz"){
	open($input,"zcat $pathhits |");
    }
    else{
	if($extension eq "bam"){
	    open($input,"samtools view $pathhits |");
	}
	else{
	    open($input,$pathhits);
	}
    }
        
    my $line=<$input>;
      
    my $nbreads=0;
    my $nbsplicedreads=0;
    my $nbsplicednonss=0;

     while($line){
	
	 chomp $line;
	 
	 if($line eq ""){
	     $line=<$input>;
	     next;
	 }

	 my $firstchar=substr $line,0,1;
	 
	 if($firstchar eq "@"){
	     $line=<$input>;
	     next;
	 }
	 
	 my @s=split("\t",$line);
	 my $cigar=$s[5];

	 my $nbhits=findNbHits(\@s);

	 if($nbhits==1){
	     
	     my $countN= ($cigar =~ tr/N//);
	     
	     if($countN>0){
	     
		 my $countD= ($cigar =~ tr/D//);
		 my $countI= ($cigar =~ tr/I//);
		 my $countM= ($cigar =~ tr/M//);
		 
		 my $nbmismatch=findMismatches(\@s);
		 my $freqmismatch=($nbmismatch+0.0)/($countM+0.0);
		 
		 if($freqmismatch<=$maxmismatch){
		     
		     my $strand=findStrand(\@s);
		     
		     if($strand eq "NA"){
			 print "Warning: no XS tag for spliced read.\n";
			 $nbsplicednonss++;
		     }
		     else{
			 my @cigarsplit1=split(/[a-zA-Z]/, $cigar);
			 my @cigarsplit2=grep { /\S/ } split(/[0-9]/, $cigar);
			 
			 my $chr=$s[2];
			 my $start=$s[3]+0;
			 
			 my $gstart=$start;
			 my $gend=$start;
			 
			 my @allgenomestart;
			 my @allgenomeend;
			 my @positionsintrons;
			 
			 my $alnlen=0;
			 
			 for(my $i=0; $i<@cigarsplit1; $i++){
			     if($cigarsplit2[$i] eq "M"){
				 $alnlen+=$cigarsplit1[$i];
				 $gend+=$cigarsplit1[$i];
			     } 
			     else{
				 if($cigarsplit2[$i] eq "D"){
				     $gend+=$cigarsplit1[$i];
				 }
				 else{
				     if($cigarsplit2[$i] eq "N"){
					 ## first part
					 
					 push(@allgenomestart, $gstart);
					 push(@allgenomeend, ($gend-1));
					 
					 my $nbcurrent=@allgenomestart;
					 push(@positionsintrons, ($nbcurrent-1)); 
					 
					 ## second part
					 
					 $gstart=$gend+$cigarsplit1[$i];
					 $gend=$gstart;
				     }
				     ## other CIGARS: S, I do not affect genome position
				 }
			     }
			 }
			 
			 ## don't forget the last segment
			 
			 push(@allgenomestart, $gstart);
			 push(@allgenomeend, ($gend-1));
			 			 
			 my $id=$s[0]; ## read id
			 
			 foreach my $k (@positionsintrons){
			     ## exon left is at $allgenomestart[$k], $allgenomeend[$k], 
			     ## exon right is at  $allgenomestart[$k+1], $allgenomeend[$k+1],
			     
			     my $startexonleft=$allgenomestart[$k];
			     my $endexonleft=$allgenomeend[$k];
			     my $startexonright=$allgenomestart[$k+1];
			     my $endexonright=$allgenomeend[$k+1];
			     
			     my $startintron=$endexonleft+1;
			     my $endintron=$startexonright-1;
			     
			     my $anchorleft=($endexonleft-$startexonleft+1);
			     my $anchorright=($endexonright-$startexonright+1);
			     
			     my $key=$startintron."_".$endintron."_".$strand; 
			     
			     if(exists $refjunctions->{$chr}){
				 if(exists $refjunctions->{$chr}{$key}){
				     push(@{$refjunctions->{$chr}{$key}{"readids"}},$id);
				     push(@{$refjunctions->{$chr}{$key}{"anchorleft"}},$anchorleft);
				     push(@{$refjunctions->{$chr}{$key}{"anchorright"}},$anchorright);
				     push(@{$refjunctions->{$chr}{$key}{"nbmismatch"}},$nbmismatch);
				 }
				 else{
				     $refjunctions->{$chr}{$key}={"readids"=>[$id],"anchorleft"=>[$anchorleft],"anchorright"=>[$anchorright],"nbmismatch"=>[$nbmismatch]};
				 }
			     }
			     else{
				 $refjunctions->{$chr}={$key=>{"readids"=>[$id],"anchorleft"=>[$anchorleft],"anchorright"=>[$anchorright],"nbmismatch"=>[$nbmismatch]}};	
			     }
			 }
		     }
		 }
	     }
	 }
	 
	 $nbreads++;

	 if($nbreads%1000000==0){
	     print $nbreads." alignments read\n";
	 }
	
	 $line=<$input>;
     }
    
    close($input);
    
    print "found ".$nbsplicedreads." spliced reads with at most ".$maxmismatch." mismatches.\n";
    print "found ".$nbsplicednonss." spliced reads without strand information.\n";
    
}

##############################################################

sub determineSpliceSite{
    my $refjunctions=$_[0];
    my $fastagenome=$_[1];
    
    foreach my $chr (keys %{$refjunctions}){
	
	my $sequence = $fastagenome->{$chr};
	    
	foreach my $keyjunc (keys %{$refjunctions->{$chr}}){
	    
	    my @s=split("_",$keyjunc);
	    
	    my $startint=$s[0]+0; 
	    my $endint=$s[1]+0;
	    
	    ## coordinates start at 1
	    
	    my $splicesite1=substr $sequence, ($startint-1),2;
	    my $splicesite2=substr $sequence, ($endint-2),2;
	    
	    my $site1=$splicesite1.$splicesite2;
	    $site1=uc $site1; 
	
	    my $site2=reverseComplement($site1);
	    
	    my $strand="NA";
	    
	    my $site=$site1."_".$site2;

	    
	    $refjunctions->{$chr}{$keyjunc}{"splicesite"}=$site;
	    	    
	    if($site eq "GTAG_CTAC"){
		$strand="1";
	    }
	    if($site eq "GCAG_CTGC"){
		$strand="1";
	    }
	    if($site eq "ATAC_GTAT"){
		$strand="1";
	    }
	    if($site eq "CTAC_GTAG"){
		$strand="-1";
	    }
	    if($site eq "CTGC_GCAG"){
		$strand="-1";
	    }
	    if($site eq "GTAT_ATAC"){
		$strand="-1";
	    }
	    $refjunctions->{$chr}{$keyjunc}{"splicestrand"}=$strand;
	    
	}
	
    }
    
}

##############################################################

sub writeJunctions{

    my $refjunctions=$_[0];
    my $pathout=$_[1];

    open(my $output,">".$pathout);

    print $output "Chromosome\tStart\tEnd\tSpliceSite\tProbableStrand\tAnchorLeft\tAnchorRight\tNbReads\n";

    foreach my $chr (keys %{$refjunctions}){
	foreach my $keyjunc (keys %{$refjunctions->{$chr}}){
	    my @s=split("_",$keyjunc);
		
	    my $startint=$s[0]+0; 
	    my $endint=$s[1]+0;
	    
	    my $nbreads=$refjunctions->{$chr}{$keyjunc}{"nbreads"};
	    my $splicesite=$refjunctions->{$chr}{$keyjunc}{"splicesite"};
	    my $strand=$refjunctions->{$chr}{$keyjunc}{"splicestrand"};
	    my $anchorleft=max (@{$refjunctions->{$chr}{$keyjunc}{"anchorleft"}});
	    my $anchorright=max (@{$refjunctions->{$chr}{$keyjunc}{"anchorright"}});
	    
	    print $output $chr."\t".$startint."\t".$endint."\t".$splicesite."\t".$strand."\t".$anchorleft."\t".$anchorright."\t".$nbreads."\n";
		
	}
    }

    close($output);

}

##############################################################

sub findMismatches{
    my $array=$_[0];
    
    my $nb=@{$array};
     
    for(my $i=11; $i<$nb; $i++){
	my @t=split(":",$array->[$i]);
	
	if($t[0] eq "NM"){
	    my $mis=$t[2]+0;
	    
	    return $mis;
	} 
    }

    ## this is in fact total edit distance, includes I & D
    print "couldn't find number of mismatches!!!\n";
    print join("\t",@{$array});
    exit(1);
}

##############################################################

sub findStrand{
    my $s=$_[0];

    my $strand="NA";
    my @grepres=grep(/XS:A/,@{$s});
    if(@grepres==1){
	my @t=split(":",$grepres[0]);
	if($t[2] eq "+"){
	    $strand="1";
	}
	else{
	    if($t[2] eq "-"){
		$strand="-1";
	    }
	}
    }

    return $strand;
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts splice junction coordinates from an alignment in SAM/BAM format.\n";
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
$parameters{"pathAln"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"anchordetection"}=8;
$parameters{"anchorquantification"}=8;
$parameters{"maxmismatch"}=0.01;
$parameters{"pathOutput"}="junctions.txt";
$parameters{"pathOutputWrongStrand"}="junctions.txt";

my @defaultpars=("pathAln","pathGenomeSequence","anchordetection","anchorquantification","maxmismatch","pathOutput","pathOutputWrongStrand");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

my %numericpars;

$numericpars{"anchordetection"}=1;
$numericpars{"anchorquantification"}=1;
$numericpars{"maxmismatch"}=1;

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

print "Extracting junctions from the accepted hits...\n";

my %junctions;

my $mis=($parameters{"maxmismatch"}+0);
print "we keep reads with at most ".$mis." mismatches/indels frequency.\n";

if($mis>=1){
    print "maxmismatch should be between 0 and 1\n";
    exit(1);
}

extractJunctions($parameters{"pathAln"},$mis,\%junctions);

my $nbtot=0;

foreach my $chr (keys %junctions){
    my $thisnb=keys %{$junctions{$chr}};
    $nbtot+=$thisnb;
}

print "There are ".$nbtot." junctions in total after the first detection step.\n";

print "Done.\n";

print "Reading genome sequence from ".$parameters{"pathGenomeSequence"}."\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"},\%genome);

my $nbchr=keys %genome;
print "There are ".$nbchr." chromosomes.\n";

print "Done.\n";

print "Extracting splice site sequence...\n";
determineSpliceSite(\%junctions,\%genome);

print "Done.\n";

print "Filtering junctions...\n";
open(my $outstrand,">".$parameters{"pathOutputWrongStrand"});

my $anchordetection=$parameters{"anchordetection"}+0;
my $anchorquantification=$parameters{"anchorquantification"}+0;

print "minimum anchor for detection: ".$anchordetection."\n";
print "minimum anchor for quantification: ".$anchorquantification."\n";

my $nbsmallanchor=0;
my $nbwrongstrand=0;
my $nbnoreads=0;

foreach my $chr (keys %junctions){
    my @keysjunc=keys %{$junctions{$chr}};
    
    foreach my $key (@keysjunc){
	my $maxleft=max @{$junctions{$chr}{$key}{"anchorleft"}};
	my $maxright=max @{$junctions{$chr}{$key}{"anchorright"}};

	if($maxleft>=$anchordetection && $maxright>=$anchordetection){
	    my @s=split("_",$key);
	    my $strandread=$s[2];
	    my $splicestrand=$junctions{$chr}{$key}{"splicestrand"};
	    my $splicesite=$junctions{$chr}{$key}{"splicesite"};

	    if($strandread eq $splicestrand){
		my $nbreads=0;
		my $nbtot=@{$junctions{$chr}{$key}{"readids"}};
				
		for(my $i=0; $i<$nbtot; $i++){
		    my $thisleft=$junctions{$chr}{$key}{"anchorleft"}[$i];
		    my $thisright=$junctions{$chr}{$key}{"anchorright"}[$i];

		    if($thisleft>=$anchorquantification && $thisright>=$anchorquantification){
			$nbreads++;
		    }
		}

		if($nbreads>0){
		    $junctions{$chr}{$key}{"nbreads"}=$nbreads;
		}
		else{
		    delete $junctions{$chr}{$key};
		    $nbnoreads++;
		}
	    }
	    else{
		$nbwrongstrand++;
		print $outstrand $chr." ".$key." strandread ".$strandread." splicestrand ".$splicestrand." splicesite ".$splicesite."\n";
		delete $junctions{$chr}{$key};
	    }
	}
	else{
	    $nbsmallanchor++;
	    delete $junctions{$chr}{$key};
	}
    }
}

print "There are ".$nbwrongstrand. " junctions with weird strands.\n";
print "There are ".$nbsmallanchor. " junctions with small anchor.\n";
print "There are ".$nbnoreads. " junctions with no reads.\n";

my $nbkept=0;

foreach my $chr (keys %junctions){
    my $thisnb=keys %{$junctions{$chr}};
    $nbkept+=$thisnb;
}

print "There are ".$nbkept." junctions after the filtering step.\n";
close($outstrand);
print "Done.\n";


print "Writing output for junctions...\n";
writeJunctions(\%junctions, $parameters{"pathOutput"});
print "Done.\n";

##############################################################
