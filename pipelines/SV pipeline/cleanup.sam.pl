use strict;

##############################################################
##############################################################

sub printHelp{
    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script makes SAM format files lighter.\n";
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
$parameters{"pathInput"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathInput", "pathOutput");

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

my $path=$parameters{"pathInput"};

if(-e $path){
    print "Reading alignments from ".$path."\n";

    open(my $output, ">".$parameters{"pathOutput"});

    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "sam"){
	open($input, $path);
    } else{
	if($ext eq "bam"){
	    open($input, "samtools view -h $path |");
	} else{
	    if($ext eq "gz"){
		open($input, "zcat $path |");
	    } else{
		print "Unknown file format\n";
		exit(1);
	    }
	}
    }

    my $line=<$input>;
  
    while($line){
	if($line eq ""){
	    $line=<$input>;
	    next;
	}
	
	my $firstchar=substr $line,0,1;
	
	if($firstchar eq "@"){
	    print $output $line; ## not chomped up yet
	    $line=<$input>;
	    next;
	}
	
	chomp $line;
	my @s=split("\t",$line);

	my $id=$s[0];
	
	## replace sequence and quality with "*", to make files lighter
	$s[9]="*";
	$s[10]="*"; 
	print $output join("\t", @s)."\n";
		
	$line=<$input>;
    }
    
    close($output);
    close($input);
}
else{
    print "Cannot find file: ".$path."\n";
}
print "Done.\n";

##############################################################
