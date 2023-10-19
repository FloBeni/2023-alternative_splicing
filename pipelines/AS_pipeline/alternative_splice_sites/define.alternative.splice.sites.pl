#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(uniq);
use strict;

##############################################################

sub readExonBlocks{
    my $pathin=$_[0];
    my $exonblocks=$_[1];

    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $gene=$s[0];
	my $chr=$s[2];
	my $start=$s[3]+0;
	my $end=$s[4]+0;
	my $strand=$s[5];

	if($strand eq "1"){
	    $strand="+";
	}

	if($strand eq "-1"){
	    $strand="-";
	}   
	
	if(exists $exonblocks->{$gene}){
	    push(@{$exonblocks->{$gene}{"start"}}, $start);
	    push(@{$exonblocks->{$gene}{"end"}}, $end);
	} else{
	    $exonblocks->{$gene}={"chr"=>$chr, "strand"=>$strand, "start"=>[$start], "end"=>[$end]};
	}
	
	$line=<$input>;	
    }
    
    
    close($input);
}

##############################################################

sub readSpliceJunctions{
    my $pathin=$_[0];
    my $junctions=$_[1];
    
    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    ## define chromosome and strand fields

    my $chrfield;

    if(exists $header{"Chromosome"}){
	$chrfield="Chromosome";
    } else{
	if(exists $header{"Chr"}){
	    $chrfield="Chr";
	} else{
	    print "Weird! cannot find chromosome field in ".$pathin."\n";
	    exit(1);
	}
    }

    ## define chromosome and strand fields

    my $strandfield;

    if(exists $header{"Strand"}){
	$strandfield="Strand";
    } else{
	if(exists $header{"ProbableStrand"}){
	    $strandfield="ProbableStrand";
	} else{
	    print "Weird! cannot find strand field in ".$pathin."\n";
	    exit(1);
	}
    }
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $chr=$s[$header{$chrfield}];
	my $strand=$s[$header{$strandfield}];

	if($strand eq "1"){
	    $strand="+";
	}

	if($strand eq "-1"){
	    $strand="-";
	} 

##	print "chr ".$chr." strand ".$strand."\n";
	
	my $start=$s[$header{"Start"}]+0; ## intron start, 1-based
	my $end=$s[$header{"End"}]+0; ## intron end, 1-based

	if($start>$end){
	    print "Weird junction coordinates: start>end in ".$line."\n";
	    exit(1);
	}
	
	my $nbreads=0;
	
	if(exists $header{"NbReads"}){
	    $nbreads=$s[$header{"NbReads"}]+0;
	}

	if(exists $junctions->{$chr}){
	    if(exists $junctions->{$chr}{$strand}){
		if(exists $junctions->{$chr}{$strand}{$start}){
		    $junctions->{$chr}{$strand}{$start}{$end}=$nbreads;
		} else{
		    $junctions->{$chr}{$strand}{$start}={$end=>$nbreads};
		}
	    } else{
		$junctions->{$chr}{$strand}={$start=>{$end=>$nbreads}};
	    }
	} else{
	    $junctions->{$chr}={$strand=>{$start=>{$end=>$nbreads}}};
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub extractAllSpliceSites{
    my $exonblocks=$_[0];
    my $splicejunctions=$_[1];
    my $splicesites=$_[2];


    my $nbgenes=0;
    
    foreach my $gene (keys %{$exonblocks}){
	my $chr=$exonblocks->{$gene}{"chr"};
	my $strand=$exonblocks->{$gene}{"strand"};

	my $nbexons=@{$exonblocks->{$gene}{"start"}};

	if($nbexons>=1){
	    if(exists $splicejunctions->{$chr}){
		if(exists $splicejunctions->{$chr}{$strand}){

	##	    print "found junctions on chr ".$chr." and strand ".$strand."\n";
		    
		    for(my $i=0; $i<$nbexons; $i++){
			my $thisexonstart=${$exonblocks->{$gene}{"start"}}[$i];
			my $thisexonend=${$exonblocks->{$gene}{"end"}}[$i];
			my $thiskeyexon=$thisexonstart.",".$thisexonend;

			for(my $junctionstart=($thisexonstart-1); $junctionstart<=($thisexonend+1); $junctionstart++){
			    if(exists $splicejunctions->{$chr}{$strand}{$junctionstart}){
				## check if the other splice site is in a different exon of the same gene
				
				foreach my $junctionend (keys %{$splicejunctions->{$chr}{$strand}{$junctionstart}}){
				    for(my $j=$i; $j<$nbexons; $j++){
					my $nextexonstart=${$exonblocks->{$gene}{"start"}}[$j];
					my $nextexonend=${$exonblocks->{$gene}{"end"}}[$j];

					if($junctionend>=($nextexonstart-1) && $junctionend<=($nextexonend+1)){
					    my $nextkeyexon=$nextexonstart.",".$nextexonend;

					    ## we attribute this junction to the exons

					    my $type1="NA";
					    my $type2="NA";
					    
					    if($strand eq "+" || $strand eq "1"){
						$type1="donor";
						$type2="acceptor";
					    } else{
						if($strand eq "-" || $strand eq "-1"){
						    $type1="acceptor";
						    $type2="donor";
						} else{
						    print "Weird strand: ".$strand."\n";
						    exit(1);
						}
					    }
						
					    if(exists $splicesites->{$gene}){
						## fill in info for exon $i
						
						if(exists $splicesites->{$gene}{$thiskeyexon}){
						    if(exists $splicesites->{$gene}{$thiskeyexon}{$junctionstart}){
							push(@{$splicesites->{$gene}{$thiskeyexon}{$junctionstart}{"otherend"}}, $junctionend);
							push(@{$splicesites->{$gene}{$thiskeyexon}{$junctionstart}{"type"}}, $type1);
						    } else{
							$splicesites->{$gene}{$thiskeyexon}{$junctionstart}={"otherend"=>[$junctionend], "type"=>[$type1]};
						    }
						} else{
						    $splicesites->{$gene}{$thiskeyexon}={$junctionstart=>{"otherend"=>[$junctionend], "type"=>[$type1]}};
						}

						## fill in info for exon $j

						if(exists $splicesites->{$gene}{$nextkeyexon}){
						    if(exists $splicesites->{$gene}{$nextkeyexon}{$junctionend}){
							push(@{$splicesites->{$gene}{$nextkeyexon}{$junctionend}{"otherend"}}, $junctionstart);
							push(@{$splicesites->{$gene}{$nextkeyexon}{$junctionend}{"type"}}, $type2);
						    } else{
							$splicesites->{$gene}{$nextkeyexon}{$junctionend}={"otherend"=>[$junctionstart], "type"=>[$type2]};
						    }
						} else{
						    $splicesites->{$gene}{$nextkeyexon}={$junctionend=>{"otherend"=>[$junctionstart], "type"=>[$type2]}};
						}
					    } else{
						## fill info for both exons
						$splicesites->{$gene}={$thiskeyexon=>{$junctionstart=>{"otherend"=>[$junctionend], "type"=>[$type1]}}, $nextkeyexon=>{$junctionend=>{"otherend"=>[$junctionstart], "type"=>[$type2]}}};
					    }

					    last;
					}
				    }
				}
			    }
			}			
		    }
		    
		}
	    }
	}

	$nbgenes++;

	if($nbgenes%1000==0){
	    print "Checked ".$nbgenes."\n";
	}
    }

}

##############################################################

sub extractAlternativeSpliceSites{
    my $exonblocks=$_[0];
    my $splicesites=$_[1]; 
    my $junctionquantification=$_[2]; ##  has number of reads information
    my $alternativesites=$_[3];

    foreach my $gene (keys %{$exonblocks}){
	my $chr=$exonblocks->{$gene}{"chr"};
	my $strand=$exonblocks->{$gene}{"strand"};

	my $nbexons=@{$exonblocks->{$gene}{"start"}};

	my $typestart="acceptor";
	my $typeend="donor";

	if($strand eq "-" || $strand eq "-1"){
	    $typestart="donor";
	    $typeend="acceptor";
	}

	if($nbexons>=1){
	    for(my $i=0; $i<$nbexons; $i++){
		my $exonstart=${$exonblocks->{$gene}{"start"}}[$i];
		my $exonend=${$exonblocks->{$gene}{"end"}}[$i];
		
		my $keyexon=$exonstart.",".$exonend;

		if(exists $splicesites->{$gene}{$keyexon}){
		    my @allpositions=keys %{$splicesites->{$gene}{$keyexon}};
		    my $nbpos=@allpositions;

		    if($nbpos>=2){
			my @sortedpositions=sort {$a<=>$b} @allpositions;

			my %acceptedsites;

			## first direction, check type start
			my $nbtypestart=0;
			
			for(my $k=0; $k<$nbpos; $k++){
			    my $thispos=$sortedpositions[$k];
			    my @thistypes=uniq @{$splicesites->{$gene}{$keyexon}{$thispos}{"type"}};

			    my $nbtypes=@thistypes;

			    if($nbtypes==1){
				my $type=${$splicesites->{$gene}{$keyexon}{$thispos}{"type"}}[0];

				if($type eq $typeend){
				    last;
				} else{
				    ## we accept this site, we compute its number of reads

				    my $totreads=0;
				    
				    foreach my $otherend (@{$splicesites->{$gene}{$keyexon}{$thispos}{"otherend"}}){
					if($thispos<$otherend && exists $junctionquantification->{$chr}{$strand}{$thispos}{$otherend}){
					    $totreads+=$junctionquantification->{$chr}{$strand}{$thispos}{$otherend};
					} else{ 
					    if($thispos>$otherend && exists $junctionquantification->{$chr}{$strand}{$otherend}{$thispos}){
						$totreads+=$junctionquantification->{$chr}{$strand}{$otherend}{$thispos};
					    }
					}
				    }
				    
				    $acceptedsites{$thispos}={"type"=>$typestart, "nbreads"=>$totreads};
				    $nbtypestart++;
				}
			    } else{
				last;
			    }
			}

			## opposite direction, check type end
			my $nbtypeend=0;

			for(my $k=($nbpos-1); $k>=0; $k--){
			    my $thispos=$sortedpositions[$k];
			    my @thistypes=uniq @{$splicesites->{$gene}{$keyexon}{$thispos}{"type"}};

			    my $nbtypes=@thistypes;

			    if($nbtypes==1){
				my $type=${$splicesites->{$gene}{$keyexon}{$thispos}{"type"}}[0];

				if($type eq $typestart){
				    last;
				} else{
				    ## we accept this site, we compute its number of reads

				    my $totreads=0;
				    
				    foreach my $otherend (@{$splicesites->{$gene}{$keyexon}{$thispos}{"otherend"}}){
					if($thispos<$otherend && exists $junctionquantification->{$chr}{$strand}{$thispos}{$otherend}){
					    $totreads+=$junctionquantification->{$chr}{$strand}{$thispos}{$otherend};
					} else{ 
					    if($thispos>$otherend && exists $junctionquantification->{$chr}{$strand}{$otherend}{$thispos}){
						$totreads+=$junctionquantification->{$chr}{$strand}{$otherend}{$thispos};
					    }
					}
				    }
				    
				    $acceptedsites{$thispos}={"type"=>$typeend, "nbreads"=>$totreads};
				    $nbtypeend++;
				}
			    } else{
				last;
			    }
			}

			## are there at least 2 sites
			
			if($nbtypestart>1 || $nbtypeend>1){
			    foreach my $site (keys %acceptedsites){
				my $thistype=$acceptedsites{$site}{"type"};
				my $nbreads=$acceptedsites{$site}{"nbreads"};
				
				if(($thistype eq $typestart && $nbtypestart>1) || ($thistype eq $typeend && $nbtypeend>1)){
				    if(exists $alternativesites->{$gene}){
					if(exists $alternativesites->{$gene}{$keyexon}){
					    if(exists $alternativesites->{$gene}{$keyexon}{$thistype}){
						$alternativesites->{$gene}{$keyexon}{$thistype}{$site}=$nbreads;
					    }  else{
						$alternativesites->{$gene}{$keyexon}{$thistype}={$site=>$nbreads};
					    }
					} else{
					    $alternativesites->{$gene}{$keyexon}={$thistype=>{$site=>$nbreads}};
					}
				    } else{
					$alternativesites->{$gene}={$keyexon=>{$thistype=>{$site=>$nbreads}}};
				    }
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
    print "This script defines alternative splice sites for genes.\n";
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
$parameters{"pathExonBlocks"}="NA";
$parameters{"pathIntronLibrary"}="NA";
$parameters{"pathIntronQuantification"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathExonBlocks", "pathIntronLibrary", "pathIntronQuantification", "pathOutput");
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

print "Reading exon blocks...\n";

my %exonblocks;

readExonBlocks($parameters{"pathExonBlocks"}, \%exonblocks);

my $nbg=keys %exonblocks;

print "Found ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Reading intron library...\n";

my %intronlibrary;

readSpliceJunctions($parameters{"pathIntronLibrary"}, \%intronlibrary);

print "Done.\n";

##############################################################

print "Reading splice junction quantification for this sample...\n";

my %splicejunctions;

readSpliceJunctions($parameters{"pathIntronQuantification"}, \%splicejunctions);
print "Done.\n";

##############################################################

print "Extracting all potential splice sites...\n";

my %allsplicesites;
extractAllSpliceSites(\%exonblocks, \%intronlibrary, \%allsplicesites);
  
print "Done.\n";

##############################################################

print "Checking and extracting alternative splice sites...\n";

my %alternativesites;

extractAlternativeSpliceSites(\%exonblocks, \%allsplicesites, \%splicejunctions, \%alternativesites);
  
print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "GeneID\tChr\tStrand\tExonBlockID\tType\tNbSites\tPositions\tNbReads\n";

foreach my $gene (keys %alternativesites){
    my $chr=$exonblocks{$gene}{"chr"};
    my $strand=$exonblocks{$gene}{"strand"};
    
    foreach my $exon (keys %{$alternativesites{$gene}}){
	foreach my $type (keys %{$alternativesites{$gene}{$exon}}){
	    my $nbsites=keys %{$alternativesites{$gene}{$exon}{$type}};

	    my @positions=sort (keys %{$alternativesites{$gene}{$exon}{$type}});
	    my @nbreads;

	    foreach my $pos (@positions){
		push(@nbreads, $alternativesites{$gene}{$exon}{$type}{$pos});
	    }

	    print $output $gene."\t".$chr."\t".$strand."\t".$exon."\t".$type."\t".$nbsites."\t".join(",", @positions)."\t".join(",",@nbreads)."\n";
	}
    }
}

close($output);

print "Done.\n";

##############################################################
