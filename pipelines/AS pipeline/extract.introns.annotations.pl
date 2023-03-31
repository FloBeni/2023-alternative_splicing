use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];
    my $txgenes=$_[2];
    
    open(my $input, $pathin);
    
    my $line=<$input>;
    
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "exon" || $type eq "CDS"){
	
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	  
	    my $newstrand=$strand;

	    if($strand eq "+"){
		$newstrand="1";
	    }
	    
	    if($strand eq "-"){
		$newstrand="-1";
	    }
	    	    
	    my $exonid=$chr.",".$start.",".$end.",".$newstrand;

	    my $info=$s[8];
	    my @t=split(";", $info);
	    my $txid=findInfo("transcript_id", \@t);

	    my $geneid=findInfo("gene_id", \@t);

	    if($geneid eq "NA"){
		$geneid=findInfo("geneID", \@t);
	    }

	    if($geneid ne "NA" && $txid ne "NA"){
		my $newid=$txid.":".$type;
		
		if(exists $txgenes->{$newid}){
		    $txgenes->{$newid}{$geneid}=1;	
		} else{
		    $txgenes->{$newid}={$geneid=>1};	
		}
		
		if(exists $transcripts->{$newid}){
		    push(@{$transcripts->{$newid}{"start"}}, $start);
		    push(@{$transcripts->{$newid}{"end"}}, $end);
		}
		else{
		    $transcripts->{$newid}={"chr"=>$chr, "strand"=>$newstrand, "start"=>[$start], "end"=>[$end]};
		}
	    }
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

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
	return "NA";
    }
}

##############################################################

sub extractIntrons{
    my $transcripts=$_[0];
    my $introns=$_[1];

    foreach my $txid (keys %{$transcripts}){
	my $nbex=@{$transcripts->{$txid}{"start"}};
	
	if($nbex>=2){
	    my $chr=$transcripts->{$txid}{"chr"};
	    my $strand=$transcripts->{$txid}{"strand"};

	    my %hashexons;

	    for(my $i=0; $i<$nbex; $i++){
		my $start=${$transcripts->{$txid}{"start"}}[$i];
		my $end=${$transcripts->{$txid}{"end"}}[$i];
		
		if(exists $hashexons{$start}){
		    print "Weird! already saw exon start ".$start." for ".$txid."\n";
		    exit(1);
		}

		$hashexons{$start}=$end;
	    }

	    my @startpos=keys %hashexons;
	    my @sortedstart=sort {$a<=>$b} @startpos;

	    for(my $i=0; $i<($nbex-1); $i++){
		my $start1=$sortedstart[$i];
		my $end1=$hashexons{$start1};
		
		my $start2=$sortedstart[$i+1];
			
		if($start2<=$end1){
		    print "Weird! overlapping exons for ".$txid."\n";
		}
		else{
		    my $startint=$end1+1;
		    my $endint=$start2-1;

		    my $keyint=$chr.",".$startint.",".$endint.",".$strand;

		    if(exists $introns->{$keyint}){
			push(@{$introns->{$keyint}}, $txid);
		    }
		    else{
			$introns->{$keyint}=[$txid];
		    }
		}
	    }
	}
    }
}

##############################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
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

################################################################

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

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts intron coordinates from annotations.\n";
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

$parameters{"pathGTF"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGTF", "pathGenomeSequence", "pathOutput");

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

print "Reading exon coordinates...\n";

my %transcripts;
my %txgenes;

readGTF($parameters{"pathGTF"}, \%transcripts, \%txgenes);

my $nbtx=keys %transcripts;

print "Found ".$nbtx." transcripts.\n";
print "Done.\n";

##############################################################

print "Extracting introns...\n";

my %introns;
extractIntrons(\%transcripts, \%introns);
my $nbint=keys %introns;

print "Found ".$nbint." distinct introns.\n";

print "Done.\n";

##############################################################

print "Reading genome sequence...\n";

my %genome;
readFasta($parameters{"pathGenomeSequence"}, \%genome);

my $nbchr=keys %genome;

print "Found ".$nbchr." chromosomes/contigs.\n";

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tStrand\tTranscripts\tGenes\t5SpliceSignal\t3SpliceSignal\n";

foreach my $int (keys %introns){
    my @s=split(",",$int);
    my $chr=$s[0];

    my $start=$s[1]+0;
    my $end=$s[2]+0;
    my $strand=$s[3];
    
    my $splice5="NA";
    my $splice3="NA";
    
    if(exists $genome{$chr}){
	if($strand eq "1"){
	    $splice5=uc(substr $genome{$chr}, ($start-1), 2);
	    $splice3=uc(substr $genome{$chr}, ($end-2), 2);
	} else{
	    if($strand eq "-1"){
		my $s3=uc(substr $genome{$chr}, ($start-1), 2);
		$splice3=reverseComplement($s3);
		
		my $s5=uc(substr $genome{$chr}, ($end-2), 2);
		$splice5=reverseComplement($s5);
	    } else{
		print "Weird strand ".$strand." for ".$int."\n";
		exit(1);
	    }
	}
	
	
	my %genes;
	foreach my $tx (@{$introns{$int}}){
	    if(!(exists $txgenes{$tx})){
		print "Weird! cannot find gene info for ".$tx."\n";
	    }
	    
	    foreach my $gene (keys %{$txgenes{$tx}}){
		$genes{$gene}=1;
	    }
	}
	
	print $output $chr."\t".$start."\t".$end."\t".$strand."\t".join(",",@{$introns{$int}})."\t".join(",",keys %genes)."\t".$splice5."\t".$splice3."\n";
	
    }
     else{
	 print "Could not find sequence for ".$chr."\n";
    }
}

close($output);

print "Done.\n";

##############################################################
