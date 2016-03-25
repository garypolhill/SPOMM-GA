#!/usr/bin/perl
#
# spomm-MGA.pl
#
# Script to run a GA over a SPOMM to find parameters that fit to an occupancy
# pattern. We use ter Braak's (2004) repopulation method, and a cost function
# based on the area underneath the 
#
# The genome consists of a vector of B parameters that are multiplied by
# a biophysical characteristics vector on each cell to give the habitat
# available for a species on that cell. Also on the genome are A and M for
# the species parameter. The GA must find B, A and M to minimise occupancy
# error against a dataset.
#
# The inputs to this program are:
#
# 1. A SPOMM parameter file and all associated files, including a species file.
#
# 2. A file containing the occupancy for each cell to compare the solution 
#    with.
#
# 3. A file or files containing the biophysical characteristics for each cell.
#
# This program uses a multicriteria fitness function.
#
# Bugs to fix 2014-05-16:
#
# (i)  The full length of the gene is not saved when using Gaussian betas
# (ii) There can be a filename clash in the Hall of Fame if lots of options
#      are saved, which causes the program to crash. PID probably can't be
#      relied on as an identifier.

use strict;
use Errno qw(EAGAIN);
use Math::Trig;
use Cwd;

if(scalar(@ARGV) < 4 && $ARGV[0] ne 'continue') {
  die "Usage: $0 [-maxproc <max concurrent processes>] ",
  "<GA info> <SPOMM file> <occupancy> <biophysical files...>\n",
  "\n",
  "\tThe GA information includes the population size, rule and rule\n",
  "\tparameters, as well as prior distributions and constraints for the\n",
  "\tsearch space. This should be stored in a GA file.\n",
  "\n",
  "\tAll files named in SPOMM parameter file are expected to exist in the\n",
  "\tsame directory, except those named *dummy*, any reporting files, and\n",
  "\tthe patchFile.\n",
  "\n",
  "\tThe occupancy and biophysical files can be in a number of formats,\n",
  "\twhich the program will try to guess.\n",
  "\n",
  "\tThe number of biophysical files determines the length of the beta\n",
  "\tvector.\n";
}

my @contpop;
my @contfit;
my @contid;
my $continuing = 0;
my @contargv;
my $contfile;
my $prevwd;
my $prevdata = 0;
if($ARGV[0] eq 'continue') {
  $continuing = 1;
  shift(@ARGV);
  $contfile = shift(@ARGV);
  my @loadargv;
  my ($cwd, $wd) = &loadga($contfile, \@contpop, \@contfit, \@contid,
			   \@loadargv);
  $prevwd = $wd;
  chdir $cwd or die "Cannot chdir $cwd: $!\n";
  push(@contargv, @ARGV);
  @ARGV = @loadargv;

  my $prevdata = 1;
  for(my $i = 0; $i <= $#contid; $i++) {
    if(!-d "$prevwd/$contid[$i]") {
      $prevdata = 0;
      last;
    }
  }
}  

# TODO Add an option to use the initial occupation probability as a searched
# parameter (it is otherwise assumed to be 1.0).

# TODO Add an option to use condor

my $maxproc = 1;
my $wd = $continuing ? $prevwd : "/var/tmp/spomm-GA/$$";
my $delete = 0;			# 1 => Delete runs
my $keephof = 1;		# 1 => Keep only hall of fame
my $cmd = "/home/gp/swarm/fearlus/model1-1-12/spom-2.3";
my $blob = "/home/gp/swarm/fearlus/expt/spomm-climate/blob.R";
my $doblob = 1;
my $blobkey = 20;
my $fitnessnooccup = 0;
my $habmodifier = 'none';
my $conflatefitness = 'no';
my $useextinct = 0;
my $usenoequi = 0;
my $nohalf = 1;
my $betaMode = 'linear';

my $log = "LOG";

my @saveargv;
push(@saveargv, @ARGV);

while($ARGV[0] =~ /^-/ || scalar(@contargv) > 0) {
  my $option;
  my $shiftarr;

  if($ARGV[0] =~ /^-/) {
    $option = shift(@ARGV);
    $shiftarr = \@ARGV;
  }
  elsif($contargv[0] =~ /^-/) {
    $option = shift(@contargv);
    $shiftarr = \@contargv;
  }
  else {
    die "Invalid argument to continuation run: $contargv[0]\n";
  }

  if($option eq '-maxproc') {
    $maxproc = shift(@$shiftarr);
  }
  elsif($option eq '-keepall') {
    $keephof = 0;
    $delete = 0;
  }
  elsif($option eq '-deleteall') {
    $delete = 1;
    $keephof = 1;
  }
  elsif($option eq '-keephof') {
    $delete = 0;
    $keephof = 1;
  }
  elsif($option eq '-habmod') {
    $habmodifier = shift(@$shiftarr);
  }
  elsif($option eq '-fitnessnooccup') {
    $fitnessnooccup = 1;
  }
  elsif($option eq '-fitnesshalfpenalty') {
    $nohalf = 0;
  }
  elsif($option eq '-conflate') {
    $conflatefitness = shift(@$shiftarr);
  }
  elsif($option eq '-useextinct') {
    $useextinct = 1;
  }
  elsif($option eq '-usenoequilibrium') {
    $usenoequi = 1;
  }
  elsif($option eq '-blob') {
    $blob = shift(@$shiftarr);
  }
  elsif($option eq '-noblob') {
    $doblob = 0;
  }
  elsif($option eq '-blobkey') {
    $blobkey = shift(@$shiftarr);
  }
  elsif($option eq '-wd') {
    $wd = shift(@$shiftarr);
  }
  elsif($option eq '-log') {
    $log = shift(@$shiftarr);
  }
  elsif($option eq '-cmd') {
    $cmd = shift(@$shiftarr);
  }
  else {
    die "Unrecognised command-line option: $option\n";
  }
}

if(!-e $wd) {
  mkdir $wd or die "Cannot create working directory $wd: $!\n";
}

$log = "$wd/$log" if $log !~ /\//;

if(!-e "$log") {
  open(LOG, ">>$log") or die "Cannot create log file $log: $!\n";
  close(LOG);
}

my $gaFile = shift(@ARGV);
my $spommFile = shift(@ARGV);
my $occupFile = shift(@ARGV);
my @biophysFiles = @ARGV;

&log("Started. GA file $gaFile, SPOMM file $spommFile, observed occupancy",
     "file $occupFile");

if($continuing) {
  &log("N.B. Continuing from file $contfile");
}

my ($npop, $ngen, $rule, $param, $search, $priors, $constraints)
  = &readga($gaFile, \@biophysFiles);

# Initialise population

my @population;
if($continuing) {		# Continuing: use loaded population
  @population = @contpop;
  if(scalar(@population) != $npop) {
    die "Cannot continue: failed to load exactly $npop instances from file ",
    "$contfile. Number loaded: ", scalar(@population), "\n";
  }
}
else {				# Not continuing: build population
  @population = &buildpop($npop, $priors, $constraints);
}

my @identifiers;
my @fitness;
if($prevdata) {			# Continuing and data from previous runs found
  @identifiers = @contid;
  @fitness = @contfit;
  if($wd ne $prevwd) {
    for(my $i = 0; $i <= $#identifiers; $i++) {
      &cp("$prevwd/$identifiers[$i]", "$wd/$identifiers[$i]");
    }
  }
}
else {				# (Re)calculate the fitness
  @fitness = &calcfitness($spommFile, $occupFile, $search, \@biophysFiles,
			  \@population, \@identifiers);
  if($continuing) {
    for(my $i = 0; $i <= $#fitness; $i++) {
      &log("Run $identifiers[$i] (previously $contid[$i]) has fitness",
	   "$fitness[$i] (previously $contfit[$i])");
    }
  }
}

my @halloffame;			# [[pop data], fitness, gen, id]

my ($ntrans, $front)
  = &maintainhalloffame(\@halloffame, \@population, \@fitness,
			0, $param, \@identifiers);
&log("Generation 0: $ntrans moved to hall of fame;",
     "best fitnesses", &frontstr($front));

# Main GA loop

for(my $gen = 1; $gen <= $ngen; $gen++) {

  my @newpopulation = &breed($rule, $param, $priors, $constraints,
			     \@population, \@fitness);

  my @newidentifiers;
  my @newfitness = &calcfitness($spommFile, $occupFile, $search,
				\@biophysFiles, \@newpopulation,
				\@newidentifiers);


  my ($ntrans, $front)
    = &maintainhalloffame(\@halloffame, \@newpopulation, \@newfitness,
			  $gen, $param, \@newidentifiers);

  &log("Generation $gen: $ntrans moved to hall of fame;",
     "best fitnesses", &frontstr($front));

  &select($rule, $param, $gen,
	  \@population, \@fitness, \@newpopulation, \@newfitness,
	  \@identifiers, \@newidentifiers);
}

&savega($wd, \@population, \@fitness, \@identifiers,
	\@saveargv, $search, \@biophysFiles);

for(my $i = 0; $i <= $#halloffame; $i++) {
  &log("Hall of fame position ", $i + 1, "fitness",
       &fitnessstr($halloffame[$i][1]),
       "found in generation $halloffame[$i][2] (run $halloffame[$i][3]):",
       &popstr($halloffame[$i][0], $search, \@biophysFiles));
}

&log("Stopped");

exit 0;

###############################################################################
#
# Subroutines
#
###############################################################################

# savega(wd, population, fitness, identifiers, argv, search, biophysfiles)
#
# Save the GA's state to a file

sub savega {
  my ($wd, $population, $fitness, $identifiers,
      $argv, $search, $biophysfiles) = @_;

  my $statefile = "$wd/state-$$.txt";

  open(STATE, ">$statefile")
    or die "Cannot create state file $statefile: $!\n";

  print STATE "CWD\n";
  print STATE getcwd(), "\n";
  print STATE "$wd\n";

  print STATE "ARGV\n";
  for(my $i = 0; $i <= $#$argv; $i++) {
    print STATE "$$argv[$i]\n";
  }
  print STATE "POPULATION\n";
  for(my $i = 0; $i <= $#$population; $i++) {
    my $popstr = &popstr($$population[$i], $search, $biophysfiles);
    my $fitstr = &fitnessstr($$fitness[$i]);
    print STATE "$$identifiers[$i]: $popstr : $fitstr\n";
  }
}

# loadga(file, population, fitness, identifiers, argv) -> cwd
#
# Load the GA's state from a file

sub loadga {
  my ($file, $population, $fitness, $identifiers, $argv) = @_;

  open(STATE, "<$file") or die "Cannot read state file $file: $!\n";

  # Read the cwd

  my $line = <STATE>;
  chomp $line;
  if($line ne 'CWD') {
    die "Error in state file $file: Expecting CWD, found $line\n";
  }
  my $cwd = <STATE>;
  chomp $cwd;
  my $wd = <STATE>;
  chomp $wd;

  # Read the argv

  $line = <STATE>;
  chomp $line;
  if($line ne 'ARGV') {
    die "Error in state file $file: Expecting ARGV, found $line\n";
  }
  while($line = <STATE>) {
    chomp $line;
    last if $line eq 'POPULATION';
    push(@$argv, $line);
  }

  # Read the population

  while($line = <STATE>) {
    chomp $line;
    if($line =~ /^(.*): (.*) : (.*)$/) {
      my $id = $1;
      my $popstr = $2;
      my $fitstr = $3;
      
      push(@$identifiers, $id);
      push(@$population, &parsepopstr($popstr));
      push(@$fitness, &parsefitnessstr($popstr));
    }
    else {
      die "Error in state file $file: Expecting population data, found ",
      "$line\n";
    }
  }

  return ($cwd, $wd);
}

# popstr(sample, search, biophysFiles) -> string
#
# Return a string containing the state of the population

sub popstr {
  my ($sample, $search, $biophysFiles) = @_;

  my $str = "[";
  my $i;
  for($i = 0; $i <= $#$search; $i++) {
    $str .= ", " if $i > 0;
    $str .= "$$search[$i] = $$sample[$i]";
  }
  for(my $j = 0; $j <= $#$biophysFiles; $j++, $i++) {
    $str .= ", $$biophysFiles[$j] = $$sample[$i]";
  }
  $str .= "]";

  return $str;
}

# parsepopstr(str) -> sample
#
# Return a reference to an array containing the state of the population

sub parsepopstr {
  my ($str) = @_;

  my @sample;

  if($str =~ /^\[.*\]$/) {
    $str = $1;
  }
  else {
    die "Invalid population string: $str\n";
  }
  $str =~ s/\s//g;
  my @samplestr = split(/,/, $str);
  for(my $i = 0; $i <= $#samplestr; $i++) {
    my ($search, $value) = split(/=/, $samplestr[$i]);

    if(!defined($value) || $value eq "") {
      die "Invalid sample string for population: $samplestr[$i]\n";
    }
    push(@sample, $value);
  }

  return \@sample;
}

# frontstr(front) -> string
#
# Return a formatted string of a pareto front of fitnesses

sub frontstr {
  my ($front) = @_;

  my $str = "{";
  for(my $i = 0; $i <= $#$front; $i++) {
    $str .= "; " if $i > 0;
    $str .= &fitnessstr($$front[$i]);
  }
  $str .= "}";

  return $str;
}

# fitnessstr(fitness) -> string
#
# Return a string containing a printable fitness

sub fitnessstr {
  my ($fitness) = @_;

  my $str = "(";
  for(my $i = 0; $i <= $#$fitness; $i++) {
    $str .= ", " if $i > 0;
    $str .= "$$fitness[$i]";
  }
  $str .= ")";

  return $str;
}

# parsefitnessstr(str) -> fitness
#
# Return a reference to a fitness array parsed from a string

sub parsefitnessstr {
  my ($str) = @_;

  if($str =~ /^\(.*\)$/) {
    $str = $1;
  }
  else {
    die "Invalid fitness string: $str\n";
  }
  $str =~ s/\s//g;
  my @fitness = split(/,/, $str);

  return \@fitness;
}

# readga(gafile, biophysFiles)
#
# The GA file contains parameters for the GA, arranged as follows:
#
# population <population>
# generations <generations>
# rule <rule>
# parameters {
#   <param name> = <param value>
# }
# searchSPOMM {
#   <SPOMM file> / <SPOMM parameter> : <prior> <constraint>
# }
# searchBETA {
#   <name[,name...]|*> : <prior> <constraint>
# }

sub readga {
  my ($gafile, $biophysFiles) = @_;

  open(GA, "<$gafile") or die "Cannot open GA file $gafile: $!\n";
  
  my $npop = &readnpop(*GA, $gafile);
  my $ngen = &readngen(*GA, $gafile);
  my $rule = &readrule(*GA, $gafile);
  my %param = &readparam(*GA, $gafile);
  my @search;
  my @priors;
  my @constraints;
  &readsearchSPOMM(*GA, $gafile, \@search, \@priors, \@constraints);
  &readsearchBETA(*GA, $gafile, $biophysFiles, \@priors, \@constraints);

  close(GA);

  my @paramtxt;
  while(my ($key, $value) = each(%param))  {
    push(@paramtxt,"$key = $value");
  }

  &log("population size $npop, number of generations $ngen, rule $rule",
       join(", ", @paramtxt), "-- searching", join(", ", @search), "and betas",
       join(", ", @$biophysFiles), "priors",
       join(", ", @priors), "constraints", join(", ", @constraints));

  return ($npop, $ngen, $rule, \%param, \@search, \@priors, \@constraints);
}

# buildpop(npop, priors, constraints) -> population
#
# Build a random new population with npop members, using the priors to sample
# and the constraints to reject samples.

sub buildpop {
  my ($npop, $priors, $constraints) = @_;

  my @population;

  for(my $pop = 0; $pop < $npop; $pop++) {
    $population[$pop] = &samplepop($priors, $constraints);
  }

  return @population;
}

# calcfitness(spommFile, occupFile, search, biophysFiles, population) -> f'ness
#
# Compute the fitness of each member of the population. This is done in
# separate steps to enable parallel execution. First all the runs are started,
# then we wait for the runs to complete, then we compute the fitness.

sub calcfitness {
  my ($spommFile, $occupFile, $search, $biophysFiles, $population,
      $identifiers) = @_;

  my @fitness;

  for(my $pop = 0; $pop <= $#$population; $pop++) {
    $$identifiers[$pop] = &startrun($spommFile, $search, $biophysFiles,
				   $$population[$pop]);
  }

  &waitforruns($identifiers);

  for(my $pop = 0; $pop <= $#$population; $pop++) {
    my $arr = &runfitness($$identifiers[$pop], $occupFile);
    if($conflatefitness ne 'no') {
      my $conf;
      if($conflatefitness eq 'sum') {
	$conf = 0;
      }
      elsif($conflatefitness =~ /^wsum\(.*\)$/) {
	$conf = 0;
      }
      elsif($conflatefitness eq 'product') {
	$conf = 1;
      }
      elsif($conflatefitness =~ /^wproduct\(.*\)$/) {
	$conf = 1;
      }
      else {
	die "Unrecognised value for -conflate argument: $conflatefitness\n";
      }

      for(my $i = 0; $i <= $#$arr; $i++) {
	if($conflatefitness eq 'sum') {
	  $conf += $$arr[$i];
	}
	elsif($conflatefitness =~ /^wsum\((.*)\)$/) {
	  my @weights = split(/,/, $1);
	  $conf += $$arr[$i] * $weights[$i];
	}
	elsif($conflatefitness =~ /^wproduct\((.*)\)$/) {
	  my @weights = split(/,/, $1);
	  $conf *= $$arr[$i] ** $weights[$i];
	}
	else {			# product
	  $conf *= $$arr[$i];
	}
      }
      $arr = [$conf];
    }
    $fitness[$pop] = $arr;

    &log("Fitness of run $$identifiers[$pop] is", &fitnessstr($fitness[$pop]));
  }

  return @fitness;
}

# breed(rule, param, priors, constraints, population, fitness) -> new pop
#
# Generate a new population of samples to explore

sub breed {
  my ($rule, $param, $priors, $constraints, $population, $fitness) = @_;

  my @newpopulation;

  for(my $i = 0; $i <= $#$population; $i++) {
    my $sample;

    do {
      if($rule =~ /^DEMCZ/) {
	$sample = &breedDEMCZ($rule, $param, $priors,
			      $population, $fitness, $i);
      }
      elsif($rule =~ /^DEMC/) {
	$sample = &breedDEMC($rule, $param, $priors,
			     $population, $fitness, $i);
      }
      elsif($rule =~ /^DE/) {
	$sample = &breedDE($rule, $param, $priors, $population, $fitness, $i);
      }
      elsif($rule =~ /^GA/) {
	my @sorted = &sortbyfitness($fitness);
	$sample = &breedGA($rule, $param, $priors, $population, $fitness, $i,
			   \@sorted);
      }
      else {
	die "Unrecognised breeder rule: $rule\n";
      }
    } while(!&constrainedall($sample, $constraints));

    $newpopulation[$i] = $sample;
  }

  return @newpopulation;
}

# select(rule, param, population, fitness, newpopulation, newfitness)
#
# Use a proposed new population to replace the current one, applying various
# different rules to determine how this is done.

sub select {
  my($rule, $param, $gen, $population, $fitness, $newpopulation, $newfitness,
     $identifiers, $newidentifiers) = @_;

  if($rule eq 'DEStrict' || $$param{'select'} eq 'posreplacebetter') {
    for(my $i = 0; $i <= $#$population; $i++) {
      my $cmp = &cmpvector($$newfitness[$i], $$fitness[$i]);
      if($cmp > 0) {
	$$population[$i] = $$newpopulation[$i];
	$$identifiers[$i] = $$newidentifiers[$i];
	$$fitness[$i] = $$newfitness[$i];
      }
    }
  }
  elsif($rule =~ /^DEMCStrict/ || $$param{'select'} eq 'posprobreplacebetter') {
    for(my $i = 0; $i <= $#$population; $i++) {
      my $r = &meanratiovector($$newfitness[$i], $$fitness[$i]);
      if(rand() < $r) {
	$$population[$i] = $$newpopulation[$i];
	$$identifiers[$i] = $$newidentifiers[$i];
	$$fitness[$i] = $$newfitness[$i];
      }
    }
  }
  elsif($rule =~ /^DEMCRlen/ || $$param{'select'} eq 'lenprobreplacebetter') {
    for(my $i = 0; $i <= $#$population; $i++) {
      my $r = &lengthratiovector($$newfitness[$i], $$fitness[$i]);
      if(rand() < $r) {
	$$population[$i] = $$newpopulation[$i];
	$$identifiers[$i] = $$newidentifiers[$i];
	$$fitness[$i] = $$newfitness[$i];
      }
    }
  }
  elsif($rule =~ /^DEMCZ/) {
    if($gen > 0 && $gen % $$param{'K'} == 0) {
      for(my $i = 0; $i <= $$param{'chains'}; $i++) {
	push(@$population, $$newpopulation[$i]);
	push(@$identifiers, $$newidentifiers[$i]);
	push(@$fitness, $$newfitness[$i]);
      }
    }
    else {
      for(my $i = 0; $i <= $$param{'chains'}; $i++) {
	my $r;
	if($rule eq 'DEMCZStrict') {
	  $r = &meanratiovector($$newfitness[$i], $$fitness[$i]);
	}
	else {
	  $r = &lengthratiovector($$newfitness[$i], $$fitness[$i]);
	}
	if(rand() < $r) {
	  $$population[$i] = $$newpopulation[$i];
	  $$identifiers[$i] = $$newidentifiers[$i];
	  $$fitness[$i] = $$newfitness[$i];
	}
      }
    }
  }
  elsif($$param{'select'} eq 'lottery') {
    &selectlottery($population, $fitness, $newpopulation, $newfitness,
		   $identifiers, $newidentifiers);
  }
  elsif($$param{'select'} eq 'lotteryall') {
    &selectlotteryall($population, $fitness, $newpopulation, $newfitness,
		      $identifiers, $newidentifiers);
  }
  elsif($$param{'select'} eq 'top') {
    &selecttop($population, $fitness, $newpopulation, $newfitness,
	       $identifiers, $newidentifiers);
  }
  elsif($$param{'select'} eq 'new') {
    for(my $i = 0; $i <= $#$population; $i++) {
      $$population[$i] = $$newpopulation[$i];
      $$identifiers[$i] = $$newidentifiers[$i];
      $$fitness[$i] = $$newfitness[$i];
    }
  }
  else {
    die "Unable to determine selection process from rule $rule and select ",
    "parameter $$param{'select'}\n";
  }
}

# samplepop(priors, constraints) -> genome
#
# Sample a single member of the population from the prior distribution, given
# the constraints. See definitions for samplegene() and constrained() to find
# out more about how these are defined.

sub samplepop {
  my ($priors, $constraints) = @_;

  my @genome;

  for(my $i = 0; $i <= $#$priors; $i++) {
    $genome[$i] = &samplegene($$priors[$i], $$constraints[$i]);
  }

  return \@genome;
}

# selectlottery(population, fitness, newpopulation, newfitness)
#
# Build the population from the new population only, allocating tickets in
# a lottery according to fitness rank.

sub selectlottery {
  my ($population, $fitness, $newpopulation, $newfitness,
      $identifiers, $newidentifiers) = @_;

  my @sorted = &sortbyfitness($newfitness);

  for(my $i = 0; $i <= $#$population; $i++) {
    my $j = &sampleranklottery(\@sorted);
    $$population[$i] = $$newpopulation[$j];
    $$identifiers[$i] = $$newidentifiers[$j];
    $$fitness[$i] = $$newfitness[$j];
  }
}

# sortbyfitness(fitness) -> array of indexes
#
# Build an array of indexes on fitness that would be used to create an array
# of fitness sorted in ascending order. Since we now have a partial ordering,
# this is done by successively finding the pareto front and 'unshifting' that
# on to the front of the sorted array of indexes. The pareto front array must
# be shuffled so there is no bias in which unordered solutions are prefered.

sub sortbyfitness {
  my ($fitness) = @_;

  my @ary;
  my @sorted;

  for(my $i = 0; $i <= $#$fitness; $i++) {
    $ary[$i] = $i;
  }

  my @front;
  do {
    @front = &shuffle(&paretofront($fitness, @sorted));
    unshift(@sorted, @front);
  } while(scalar(@front) > 0);

  return @sorted;
}

# shuffle(arr) -> shuffled arr
#
# Shuffle the elements in the array. This relies on the array having no two
# items that are equal.

sub shuffle {
  my (@arr) = @_;

  my @shuffler;
  my %indexes;
  for(my $i = 0; $i <= $#arr; $i++) {
    $indexes{$arr[$i]} = $i;
    
    $shuffler[$i] = rand();
  }

  my @shuffled = sort { $shuffler[$indexes{$a}]
			  <=> $shuffler[$indexes{$b}] } @arr;
  return @shuffled;
}

# paretofront(fitness, notarr) -> pareto front
#
# Find the set of incomparable best fitnesses in the fitness array, with the
# exception of those in notarr. These may be assumed to be fitnesses that have
# already been found. The algorithm is possibly not the most efficient.

sub paretofront {
  my ($fitness, @notarr) = @_;

  my @front;
  my %not;

  foreach my $notthis (@notarr) {
    $not{$notthis} = $notthis;
  }

  my $i;
  for($i = 0; $i <= $#$fitness; $i++) {
    if(!defined($not{$i})) {
      push(@front, $i);
      last;
    }
  }
  for(; $i <= $#$fitness; $i++) {
    next if defined($not{$i});
    my $allincomparable = 1;
    my $added = 0;
    for(my $j = 0; $j <= $#front; $j++) {
      my $cmp = &pocmpvector($$fitness[$i], $$fitness[$front[$j]]);
      if($cmp ne 'incomparable') {
	$allincomparable = 0;
	if($cmp > 0) {
	  if(!$added) {
	    $front[$j] = $i;
	    $added = 1;
	  }
	  else {
	    splice(@front, $j, 1);
	  }
	}
      }
    }
    push(@front, $i) if $allincomparable;
  }

  return @front;
}
    

# sampleranklottery(sorted) -> sample index
#
# Select from an array of sorted indexes using a lottery where i tickets are
# given to the ith member (thus member N, the best member, gets N tickets, 
# whilst the worst member gets 1 ticket).

sub sampleranklottery {
  my ($sorted) = @_;

  my $popsize = scalar(@$sorted);
  my $ntickets = $popsize * ($popsize + 1) / 2.0;

  my $ticket = rand() * $ntickets;
  my $tticket = 0;
  for(my $j = 0; $j <= $#$sorted; $j++) {
    $tticket += $j + 1;
    if($tticket >= $ticket) {
      return $$sorted[$j];
    }
  }
  return $$sorted[$#$sorted];
}

# selecttop(population, fitness, newpopulation, newfitness
#
# The new population consists of the top N members of population and
# newpopoulation, where N is the population size.

sub selecttop {
  my ($population, $fitness, $newpopulation, $newfitness,
      $identifiers, $newidentifiers) = @_;

  my @wholepopulation = (@$population, @$newpopulation);
  my @wholefitness = (@$fitness, @$newfitness);
  my @wholeidentifiers = (@$identifiers, @$newidentifiers);

  my @sorted = &sortbyfitness(\@wholefitness);

  for(my $i = 0; $i <= $#$population; $i++) {
    $$population[$i] = $wholepopulation[$sorted[$i]];
    $$identifiers[$i] = $wholeidentifiers[$sorted[$i]];
    $$fitness[$i] = $wholefitness[$sorted[$i]];
  }
}

# breedDEMCZ(rule, param, priors, population, fitness, i) -> baby
#
# Use the breeder rule (2) in ter Braak (2008) to determine the new sample
# to explore

sub breedDEMCZ {
  my ($rule, $param, $priors, $population, $fitness, $i) = @_;

  my $popsize = scalar(@$population) - $$param{'chains'};
  my $R1 = int(rand($popsize - 1));
  $R1++ if $R1 >= $i;
  my $R2 = int(rand($popsize - 2));
  $R2++ if $R2 >= ($i < $R1 ? $i : $R1);
  $R2++ if $R2 >= ($i > $R1 ? $i : $R1);

  $R1 += $$param{'chains'};
  $R2 += $$param{'chains'};

  my $gene = &subvector($$population[$R1], $$population[$R2]);
  $gene = &scalevector($gene, $$param{'gamma'});
  $gene = &addvector($$population[$i], $gene);
  if($rule =~ /Normal/) {
    $gene = &perturbnormalvector($gene, $$param{'b'});
  }
  else {
    $gene = &perturbvector($gene, $$param{'b'});
  }

  if(defined($$param{'CR'})) {
    $gene = &crossoverDEMC($gene, $$population[$i], $$param{'CR'});
  }

  $gene = &normalisevector($gene, $priors);

  return $gene;
}

# breedDEMC(rule, param, priors, population, fitness, i) -> baby
#
# Use the breeder rule (2) in ter Braak (2004) to determine the new sample
# to explore

sub breedDEMC {
  my ($rule, $param, $priors, $population, $fitness, $i) = @_;

  my $popsize = scalar(@$population);
  my $R1 = int(rand($popsize - 1));
  $R1++ if $R1 >= $i;
  my $R2 = int(rand($popsize - 2));
  $R2++ if $R2 >= ($i < $R1 ? $i : $R1);
  $R2++ if $R2 >= ($i > $R1 ? $i : $R1);

  my $gene = &subvector($$population[$R1], $$population[$R2]);
  $gene = &scalevector($gene, $$param{'gamma'});
  $gene = &addvector($$population[$i], $gene);
  if($rule =~ /Normal/) {
    $gene = &perturbnormalvector($gene, $$param{'b'});
  }
  else {
    $gene = &perturbvector($gene, $$param{'b'});
  }

  if(defined($$param{'CR'})) {
    $gene = &crossoverDEMC($gene, $$population[$i], $$param{'CR'});
  }

  $gene = &normalisevector($gene, $priors);

  return $gene;
}

# breedDE(rule, param, priors, population, fitness, i) -> baby
#
# Use the more traditional differential evolution rule (1) in ter Braak (2004)
# to build the proposal

sub breedDE {
  my ($rule, $param, $priors, $population, $fitness, $i) = @_;

  my $popsize = scalar(@$population);
  my $R0 = int(rand($popsize));
  my $R1 = int(rand($popsize));
  my $R2 = int(rand($popsize));

  my $gene = &subvector($$population[$R1], $$population[$R2]);
  $gene = &scalevector($gene, $$param{'gamma'});
  $gene = &addvector($$population[$R0], $gene);

  $gene = &normalisevector($gene, $priors);

  return $gene;
}

# breedGA(rule, param, priors, population, fitness, i) -> baby
#
# Use a more traditional GA breeding rule, which generates a new sample by
# applying genetic operators to selected parents.

sub breedGA {
  my ($rule, $param, $priors, $population, $fitness, $i, $sorted) = @_;

  my $S1 = &sampleranklottery($sorted);
  my $S2 = &sampleranklottery($sorted);

  my $gene = ((rand() < $$param{'pcrossover'})
	      ? &crossover($$population[$S1], $$population[$S2])
	      : (rand() < 0.5 ? $$population[$S1] : $$population[$S2]));
  $gene = &mutate($gene, $$param{'pmutate'}, $priors);

  if(defined($$param{'b'})) {
    my $pperturb = defined($$param{'pperturb'}) ? $$param{'pperturb'} : 1;
    $gene = &perturb($gene, $pperturb, $$param{'b'});
  }

  return $gene;
}

# crossoverDEMC(parent1, parent2, cr) -> child
#
# Simple crossover from ter Braak (2004): take element j from parent1 with
# probability cr, otherwise from parent2.

sub crossoverDEMC {
  my ($parent1, $parent2, $cr) = @_;

  my @gene;

  for(my $j = 0; $j <= $#$parent1; $j++) {
    $gene[$j] = (rand() < $cr) ? $$parent1[$j] : $$parent2[$j];
  }

  return \@gene;
}

# crossover(parent1, parent2) -> child
#
# Traditional crossover: choose a crossover point randomly in the genome and
# take all the elements from parent1 before that point, and all the elements
# from parent2 thereafter.

sub crossover {
  my ($parent1, $parent2) = @_;

  my $point = int(rand(scalar(@$parent1) + 1));

  my @child;

  for(my $i = 0; $i <= $#$parent1; $i++) {
    $child[$i] = ($i < $point) ? $$parent1[$i] : $$parent2[$i];
  }

  return \@child;
}

# mutate(gene, mp, priors) -> mutant
#
# Mutate the gene by resampling some of its elements from their corresponding
# prior distributions.

sub mutate {
  my ($gene, $mp, $priors) = @_;

  my @mutant;

  for(my $i = 0; $i <= $#$gene; $i++) {
    $mutant[$i] = (rand() < $mp) ? &samplegene($$priors[$i], 'X') : $$gene[$i];
  }

  return \@mutant;
}

# perturb(gene, pp, amount) -> perturbed
#
# Perturb the gene a little by adding a random uniform quantity in the range
# -amount to +amount to some of its elements.

sub perturb {
  my ($gene, $pp, $amount) = @_;

  my $fullyperturbed = &perturbvector($gene, $amount);
  my @perturbed;

  for(my $i = 0; $i <= $#$gene; $i++) {
    $perturbed[$i] = (rand() < $pp) ? $$fullyperturbed[$i] : $$gene[$i];
  }

  return \@perturbed;
}

# subvector(vec1, vec2) -> difference
#
# Create a new vector equal to vec1 - vec2, which are assumed to be of the
# same dimensionality.

sub subvector {
  my ($vec1, $vec2) = @_;

  my @ans;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    $ans[$i] = $$vec1[$i] - $$vec2[$i];
  }
  
  return \@ans;
}

# addvector(vec1, vec2) -> added
#
# Create a new vector equal to vec1 + vec2, which are assumed to be of the
# same dimensionality.

sub addvector {
  my ($vec1, $vec2) = @_;

  my @ans;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    $ans[$i] = $$vec1[$i] + $$vec2[$i];
  }

  return \@ans;
}

# mulvector(vec1, vec2) -> dot product
#
# Return the dot product of the two vectors

sub mulvector {
  my ($vec1, $vec2) = @_;

  my $ans = 0;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    $ans += $$vec1[$i] * $$vec2[$i];
  }

  return $ans;
}

# meanratiovector(vec1, vec2) -> mean ratio
#
# Return the mean ratio of the elements of the two vectors

sub meanratiovector {
  my ($vec1, $vec2) = @_;

  my $tratio = 0;
  my $n = 0;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    $tratio += $$vec1[$i] / $$vec2[$i];
    $n++;
  }

  return $tratio / $n;
}

# lengthratiovector(vec1, vec2) -> length of the vector of ratios
#
# Return the length of the vector of ratios

sub lengthratiovector {
  my ($vec1, $vec2) = @_;

  my $ssratio = 0;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    my $ratio = $$vec1[$i] / $$vec2[$i];
    $ssratio += $ratio * $ratio;
  }

  return sqrt($ssratio);
}

# cmpvector(vec1, vec2) -> comparison
#
# Return a comparison of two vectors. This is -1 if all elements of vec1 are
# less than or equal to their corresponding element in vec2; +1 if the relation
# is greater than or equal to; and 0 if the vectors are equal or incomparable.
# This method may not be suitable for using in sorting algorithms, as the
# semantics of a zero return value would break the properties of == (i.e. that
# if A == B and A == C then B == C; A could be incomparable to B and C, and 
# B > C, say)

sub cmpvector {
  my ($vec1, $vec2) = @_;

  my $ans = pocmpvector($vec1, $vec2);
  return $ans eq 'incomparable' ? 0 : $ans;
}

# pocmpvector(vec1, vec2) -> comparison
#
# Returns a partial order comparison of two vectors. If all elements of vec1
# are equal to their corresponding element in vec2, 0 is returned; if all
# elements of vec1 are greater than or equal to their corresponding element
# in vec2 then +1; if less than, -1. Otherwise, the string 'incomparable'
# is returned.

sub pocmpvector {
  my ($vec1, $vec2) = @_;

  my $ans = 0;
  for(my $i = 0; $i <= $#$vec1; $i++) {
    my $cmp = $$vec1[$i] <=> $$vec2[$i];
    if(($ans < 0 && $cmp > 0) || ($ans > 0 && $cmp < 0)) {
      return 'incomparable';
    }
    elsif($ans == 0) {
      $ans = $cmp;
    }
  }

  return $ans;
}

# scalevector(vec, scale) -> scaled
#
# Create a new vector equal to scale * vec

sub scalevector {
  my ($vec, $scale) = @_;

  my @ans;

  for(my $i = 0; $i <= $#$vec; $i++) {
    $ans[$i] = $scale * $$vec[$i];
  }

  return \@ans;
}

# perturbvector(vec, amount) -> perturbed
#
# Create a new vector equal to vec + [U(-amount, +amount)]^d, where d is the
# number of dimensions of vec

sub perturbvector {
  my ($vec, $amount) = @_;

  my @ans;

  for(my $i = 0; $i <= $#$vec; $i++) {
    $ans[$i] = $$vec[$i] + &sampleuniform(-$amount, $amount);
  }

  return \@ans;
}

# perturbnormalvector(vec, var) -> perturbed
#
# Create a new vector equal to vec + [N(0, var)]^d, where d is the
# number of dimensions of vec

sub perturbnormalvector {
  my ($vec, $var) = @_;

  my @ans;

  for(my $i = 0; $i <= $#$vec; $i++) {
    $ans[$i] = $$vec[$i] + &samplenormal(0, $var);
  }

  return \@ans;
}

# eqvector(vec1, vec2) -> boolean
#
# Return whether or not the two vectors (assumed to have the same number of
# dimensions) are numerically equal.

sub eqvector {
  my ($vec1, $vec2) = @_;

  for(my $i = 0; $i <= $#$vec1; $i++) {
    return 0 if($$vec1[$i] != $$vec2[$i]);
  }

  return 1;
}

# normalisevector(vec, priors) -> normalised
#
# Apply normalisation parameters (if any) from the priors to a sample

sub normalisevector {
  my ($vec, $priors) = @_;

  my @normalised;

  for(my $i = 0; $i <= $#$vec; $i++) {
    $normalised[$i] = &normalise($$vec[$i], $$priors[$i]);
  }

  return \@normalised;
}

# normalise(value, prior) -> normalised
#
# Normalise a single genome according to its prior

sub normalise {
  my ($value, $prior) = @_;

  my ($mode, $priord, $param) = &getnormparams($prior);

  if($mode eq 'none') {
    return $value;
  }
  else {
    return &normalisesample($value, $mode, $param);
  }
}

# getnormparams(prior) -> (mode, prior, param)
#
# Extract the normalisation mode and parameters for it from the prior,
# and return the contained prior too

sub getnormparams {
  my ($prior) = @_;

  my $normalise = 'none';
  my @normalparam;

  my $priord = $prior;

  if($prior =~ /^([STM])\[(.*);(-?\d*\.?\d+);(-?\d*\.?\d+)\]$/) {
    $normalise = $1;
    $priord = $2;
    push(@normalparam, $3, $4);
  }
  elsif($prior =~ /([UL][VEe])\[(.*);(-?\d*\.?\d+)\]$/) {
    $normalise = $1;
    $priord = $2;
    push(@normalparam, $3);
  }

  return($normalise, $priord, \@normalparam);
}

# samplegene(prior, constraint) -> sample
#
# Provide a sample for a single gene from the prior. The prior is a string
# formatted U(min,max) or N(mean,var) from which samples are taken, and
# the constraint is a satisfaction criterion for population membership,
# specifying a range the value of the gene may take (see comments to
# constrained())

sub samplegene {
  my ($sprior, $constraint) = @_;

  my ($normalise, $prior, $param)  = &getnormparams($sprior);

  my $gene;
  do {
    if($prior =~ /^U\((-?\d*\.?\d+),(-?\d*\.?\d+)\)$/) {
      $gene = &sampleuniform($1, $2);
    }
    elsif($prior =~ /^N\((-?\d*\.?\d+),(-?\d*\.?\d+)\)$/) {
      $gene = &samplenormal($1, $2);
    }
    elsif($prior =~ /^E\((-?\d*\.?\d+)\)$/) {
      $gene = &sampleexponential($1);
    }
    else {
      die "Invalid format for prior: $prior (from $sprior)\n";
    }
    if($normalise ne 'none') {
      $gene = &normalisesample($gene, $normalise, $param);
    }
  } while(!&constrained($gene, $constraint));

  return $gene;
}

# sampleuniform(min, max) -> sample
#
# Return a sample from a uniform distribution.

sub sampleuniform {
  my ($min, $max) = @_;

  die "Invalid uniform sample range U($min,$max)\n" if($max < $min);

  return (rand() * ($max - $min)) + $min;
}

# samplenormal(mean, var) -> sample
#
# Return a sample from a normal distribution, using the Box-Muller 
# transform.

BEGIN {
  # This is the Perl idiom for the equivalent in C of static variables in
  # function definition:

  my $normalsampleswitch;
  my $savednormalsample;


  sub samplenormal {
    my ($mean, $var) = @_;

    if($normalsampleswitch) {
      $normalsampleswitch = 0;
      return (sqrt($var) * $savednormalsample) + $mean;
    }

    my $sample1 = rand();
    $sample1 = 1 if $sample1 == 0;

    my $sample2 = rand();
    $sample2 = 1 if $sample2 == 0;

    $normalsampleswitch = 1;
    $savednormalsample = sqrt(-2.0 * log($sample1)) * sin(pi * 2 * $sample2);

    return (sqrt($var) * sqrt(-2.0 * log($sample1))
	    * cos(pi * 2 * $sample2)) + $mean;
  }
}

# sampleexponential(lambda) -> sample
#
# Return a sample from an exponential distribution

sub sampleexponential {
  my ($lambda) = @_;

  my $usample = rand();
  $usample = 1 if $usample == 0;
  return (-log($usample)) / $lambda;
}

# normalisesample(value, method, nparam) -> normalised sample
#
# Use various methods to impose bounds on samples mathematically. This is
# useful e.g. when sampling from a normal distribution.

sub normalisesample {
  my ($value, $method, $nparam) = @_;
  
  if($method eq 'none') {
    return $value;
  }
  if($method eq 'T' || $method eq 'S' || $method eq 'M') {
    my($min, $max) = @$nparam;

    if($min >= $max) {
      die "Invalid normalisation parameters for sigmoid: ($min, $max)\n";
    }

    if($method eq 'T') {
      return $min + (0.5 + ((($max - $min) / 2) * tanh($value)));
    }
    elsif($method eq 'M') {
      return ($value < $min) ? $min : (($value > $max) ? $max : $value);
    }
    else {
      return $min + (($max - $min) / (1.0 + exp(-$value)));
    }
  }
  if($method eq 'LV' || $method eq 'LE' || $method eq 'Le') {
    my($min) = @$nparam;

    if($method eq 'LV') {
      return $min + abs($value);
    }
    elsif($method eq 'Le') {
      return $min + exp(-$value);
    }
    else {
      return $min + exp($value);
    }
  }
  elsif($method eq 'UV' || $method eq 'UE' || $method eq 'Ue') {
    my($max) = @$param;

    if($method eq 'UV') {
      return $max - abs($value);
    }
    elsif($method eq 'Ue') {
      return $max - exp(-$value);
    }
    else {
      return $max - exp($value);
    }
  }

  die "Invalid normalisation method: $method\n";
}

# constrained(value, constraint) -> boolean
#
# Return 1 if the value is constrained by the constraint, and 0 otherwise.
# The constraint has the format 'none' or 'X' if there are no constraints,
# otherwise <fp><cmpl>X<cmpl><fp>, X<cmpl><fp>, or X<cmpg><fp>, where <cmpl>
# is one of < or <=, <cmpg> one of > or >=, and <fp> is a floating point
# number (without scientific notation).

sub constrained {
  my ($value, $constraint) = @_;

  return 1 if($constraint eq 'none' || $constraint eq 'X');

  if($constraint !~ /^(-?\d*\.?\d+<=?)?X(<=?-?\d*\.?\d+)?$/
     && $constraint !~ /^X>=?-?\d*\.?\d+$/) {
    die "Invalid constraint format: $constraint\n";
  }

  return 0 if($constraint =~ /X<(-?\d*\.?\d+)$/ && $value >= $1);
  return 0 if($constraint =~ /X<=(-?\d*\.?\d+)$/ && $value > $1);

  return 0 if($constraint =~ /^(-?\d*\.?\d+)<X/ && $value <= $1);
  return 0 if($constraint =~ /^(-?\d*\.?\d+)<=X/ && $value < $1);

  return 0 if($constraint =~ /^X>(-?\d*\.?\d+)$/ && $value <= $1);
  return 0 if($constraint =~ /^X>=(-?\d*\.?\d+)$/ && $value < $1);

  return 1;
}

# contrainedall(sample, constraints) -> boolean
#
# Check whether a whole sample is constrained as specified

sub constrainedall {
  my ($sample, $constraints) = @_;

  for(my $i = 0; $i <= $#$sample; $i++) {
    return 0 if(!&constrained($$sample[$i], $$constraints[$i]));
  }

  return 1;
}

# maintainhalloffame(halloffame, population, fitness, gen, param) -> ntrans
#
# Keep a record of all the best samples found, their fitness, and the
# generation they were first discovered. Return the number of members of the
# population transferred to the hall of fame.
#
# This method works by finding the pareto front in the current hall of fame
# and the population taken together. The hall of fame is then made to consist
# of at least that pareto front. If the size of the pareto front is less than
# the size of the hall of fame, then old members of the hall of fame are
# retained.

sub maintainhalloffame {
  my ($halloffame, $population, $fitness, $gen, $param, $identifiers) = @_;

  # Build arrays containing all information about the population and the
  # hall of fame

  my @wholepopulation;
  my @wholefitness;
  my @wholegen;
  my @wholeid;

  my $i;
  for(my $j = 0, $i = 0; $j <= $#$population; $j++, $i++) {
    push(@wholepopulation, $$population[$j]);
    push(@wholefitness, $$fitness[$j]);
    push(@wholegen, $gen);
    push(@wholeid, $$identifiers[$j]);
  }

  my $hofstart = $i;		# Remember where we started adding the hall
				# of fame

  for(my $j = 0; $j <= $#$halloffame; $j++) {
    push(@wholepopulation, $$halloffame[$j][0]);
    push(@wholefitness, $$halloffame[$j][1]);
    push(@wholegen, $$halloffame[$j][2]);
    push(@wholeid, $$halloffame[$j][3]);
  }

  # Get the pareto front of the hall of fame union the population

  my @front = &shuffle(&paretofront(\@wholefitness, ()));
				# Shuffle the front so when members are
				# removed from the hall of fame there is no
				# bias
  my @fitnessfront;
  for(my $j = 0; $j <= $#front; $j++) {
    $fitnessfront[$j] = $wholefitness[$front[$j]];
  }

  # Make sure the hall of fame directory is there if required

  if($keephof && !-e "$wd/hof") {
    mkdir("$wd/hof") or die "Cannot create directory $wd/hof: $!\n";
  }

  # Find out the number of members of the front that are in the hall of fame
  # and remember which they are

  my $ninhof = 0;
  my %inhof;
  for(my $k = 0; $k <= $#front; $k++) {
    if($front[$k] >= $hofstart) {
      $ninhof++;
      my $hofref = $front[$k] - $hofstart;
      $inhof{$hofref} = $k;
    }
  }

  # Calculate the hall of fame size. This will be bigger than the 'keep'
  # parameter if the latter is smaller than the size of the pareto front

  my $nhalloffame = (defined($$param{'keep'})
		     ? $$param{'keep'}
		     : scalar(@$population));

  $nhalloffame = scalar(@front) if scalar(@front) > $nhalloffame;

  # Remove members of the hall of fame from the front of the array unless
  # they are in the pareto front, until the hall of fame is small enough
  # to add members of the population in the pareto front to the hall of fame

  $i = 0;
  my @halloffamefront;
  while(scalar(@$halloffame) + scalar(@front)
	+ scalar(@halloffamefront) - $ninhof > $nhalloffame) {
    my $hofmember = shift(@$halloffame);
    if(!defined($inhof{$i})) {
      &rm("$wd/hof/$$hofmember[3]") if $keephof;
    }
    else {
      push(@halloffamefront, $hofmember);
    }
    $i++;
  }
  unshift(@$halloffame, @halloffamefront);

  # Add members of the population in the pareto front to the hall of fame,
  # keeping track of which and how many there were

  my %transhof;
  my $ntrans = 0;
  for(my $k = 0; $k <= $#front; $k++) {
    if($front[$k] < $hofstart) {
      &mvhof($$identifiers[$front[$k]]) if $keephof;
      $transhof{$front[$k]} = $front[$k];
      $ntrans++;

      push(@$halloffame, [$$population[$front[$k]],
			  $$fitness[$front[$k]],
			  $gen, $$identifiers[$front[$k]]]);
    }
  }

  # Remove data stored for members of the population not in the hall of fame,
  # if required

  if($keephof) {
    for(my $j = 0; $j <= $#$identifiers; $j++) {
      &rm("$wd/$$identifiers[$j]") unless defined($transhof{$j});
    }
  }

  # Return the number of members of the population transferred to the hall of
  # fame

  return ($ntrans, \@fitnessfront);
}

# mvhof(identifer)
#
# Save a run to the hall of fame

sub mvhof {
  my ($identifier) = @_;

  my $hofname = "$wd/hof/$identifier";
  my $c = -1;
  while(-e "$hofname") {
    $c++;
    $hofname = "$wd/hof/$identifier-$c";
  }
  &mv("$wd/$identifier", $hofname);
}

# readnpop(fp, gafile)
#
# Read the population size from the GA file pointer

sub readnpop {
  my ($fp, $gafile) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading population size\n";
  }
  chomp $line;
  if($line =~ /^population\s+(\d+)$/) {
    return $1;
  }
  else {
    die "Error in GA file $gafile: expecting \"population N\" found $line\n";
  }
}

# readngen(fp, gafile)
#
# Read the number of generations from the GA file pointer

sub readngen {
  my ($fp, $gafile) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading number of ",
    "generations\n";
  }
  chomp $line;
  if($line =~ /^generations\s+(\d+)$/) {
    return $1;
  }
  else {
    die "Error in GA file $gafile: expecting \"generations N\" found $line\n";
  }
}

# readrule(fp, gafile)
#
# Read the rule from the GA file pointer

sub readrule {
  my ($fp, $gafile) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading breeding rule\n";
  }
  chomp $line;
  if($line =~ /^rule\s+(\S+)$/) {
    return $1;
  }
  else {
    die "Error in GA file $gafile: expecting \"rule S\" found $line\n";
  }
}

# readparam(fp, gafile)
#
# Read the GA rule parameters from the GA file pointer

sub readparam {
  my ($fp, $gafile) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading breeding rule\n";
  }
  chomp $line;
  if($line !~ /^parameters\s+\{$/) {
    die "Error in GA file $gafile: expecting \"parameters {\" found $line\n";
  }
  my %param;
  while($line = <$fp>) {
    chomp $line;
    last if($line eq '}');
    my @words = split(" ", $line);
    if($words[1] ne '=' || scalar(@words) != 3) {
      die "Error in GA file $gafile: expecting \"<param> = <value>\" found ",
      "$line\n";
    }
    $param{$words[0]} = $words[2];
  }
  if($line ne '}') {
    die "Error in GA file $gafile: expecting \"}\", found EOF\n";
  }
  return %param;
}

# readsearchSPOMM(fp, gafile, search, priors, constraints)
#
# Read SPOMM search parameters from the GA file pointer

sub readsearchSPOMM {
  my ($fp, $gafile, $search, $priors, $constraints) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading SPOMM search ",
    "parameters\n";
  }
  chomp $line;
  if($line !~ /^searchSPOMM\s+\{$/) {
    die "Error in GA file $gafile: expecting \"searchSPOMM {\" found $line\n";
  }

  while($line = <$fp>) {
    chomp $line;
    last if $line eq '}';
    my @words = split(" ", $line);
    if(scalar(@words != 6) || $words[1] ne '/' || $words[3] ne ':') {
      die "Erorr in GA file $gafile: expecting \"<SPOMM file> / <SPOMM ",
      "parameter> : <prior> <constraint>\" found $line\n";
    }
    push(@$search, "$words[0]/$words[2]");
    push(@$priors, $words[4]);
    push(@$constraints, $words[5]);
  }

  if($line ne '}') {
    die "Error in GA file $gafile: expecting \"}\" found EOF\n";
  }
}

# combinations(n, r) -> c
#
# Return n! / r!(n - r)!
#
# According to this webpage: http://blog.plover.com/math/choose.html
# This is best done observing that
#
# C(n + 1, r + 1) = ((n + 1) / (r + 1)) * C(n, r)
#
# With C(n, 0) = 1, and c(0, r) = 0

sub combinations {
  my ($n, $r) = 0;

  return 0 if($r > $n);

  my $c = 1;

  $r = $n - $r if($r > $n / 2);

  for(my $i = $r; $i >= 1; $i--, $n--) {
    $c *= $n / $i;	     
				# The website recommends multiplying by n
				# then dividing by i to avoid rounding errors
				# I don't really care here...
  }

  return $c;
}

# readsearchBETA(fp, gafile, biophysFiles, priors, constraints)
#
# Read beta vector search parameters from the GA file pointer

sub readsearchBETA {
  my ($fp, $gafile, $biophysFiles, $priors, $constraints) = @_;

  my $line = <$fp>;
  if(!$line) {
    die "Unexpected EOF in GA file $gafile while reading beta vector search ",
    "parameters\n";
  }
  chomp $line;

  my $mode;
  my $intercept;
  if($line =~ /^searchBETA\s+\{$/) {
    $mode = 1;
    $intercept = 1;
  }
  elsif($line =~ /^searchGaussianBETA\s+\{$/) {
    $mode = 2;
    $betaMode = 'Gaussian';
    $intercept = 0;
  }
  elsif($line =~ /^searchPoly(\d+)BETA\s+\{$/) {
    $mode = $1;
    $betaMode = "polynomial,$mode";
    $intercept = 1;
  }
  elsif($line =~ /^searchSpline(\d+)BETA\s+\{$/) {
    $mode = $1;
    $betaMode = "spline,$mode";
    $intercept = 0;
  }
  else {
    die "Error in GA file $gafile: expecting \"searchBETA {\" or ",
      "\"searchGaussianBETA\", found $line\n";
  }
  # The number of BETAs there should be is $mode * scalar(@$biophysfiles) +
  # $intercept

  my $nBETA = 0;
  my %remainingBETAs;		# Betas for which priors and constraints have
				# not yet been defined
  for(my $j = 1; $j <= $mode; $j++) {
    for(my $i = 0; $i <= $#$biophysFiles; $i++) {
      # IMPORTANT AS IT AFFECTS BETA ORDERING
      #
      # This line (and that below) sets up how the betas are ordered:
      #
      # mode index = 1, biophys files 1 to N
      # mode index = 2, biophys files 1 to N
      # ...
      # mode index = M, biophys files 1 to N
      # intercept
      #
      # -- where N is the number of biophysical variables, and M is 
      # the mode number $mode. See comments later for Gaussian
      # (mode index 1 == mean; mode index 2 == std. dev)
      #
      # When reconstructed in buildpatches(), this ordering needs to
      # be respected so that the priors and constraints are respected
      $remainingBETAs{$$biophysFiles[$i], $j}
      = (($j - 1) * scalar(@$biophysFiles)) + $i;
    }
  }
  if($intercept) {
    # IMPORTANT: SEE ABOVE
    $remainingBETAs{'intercept', 0} = $mode * scalar(@$biophysFiles);
  }
  my @betapriors;
  my @betaconstraints;
  while($line = <$fp>) {
    chomp $line;
    last if $line eq '}';

    if($nBETA == $mode * scalar(@$biophysFiles) + $intercept) {
      die "Error in GA file $gafile: expecting \"}\" as priors and ",
      "constraints defined for all $nBETA beta parameters, found $line\n";
    }

    my @words = split(" ", $line);
    if(scalar(@words) != 4 || $words[1] ne ':') {
      die "Error in GA file $gafile: expecting <name[,name...]|*> : <prior> ",
      "<constraint> found $line\n";
    }
    
    if($words[0] eq '*') {
      for(my $j = 1; $j <= $mode; $j++) {
	for(my $i = 0; $i <= $#$biophysFiles; $i++) {
	  if(defined($remainingBETAs{$$biophysFiles[$i], $j})) {
	    my $k = $remainingBETAs{$$biophysFiles[$i], $j};

	    $betapriors[$k] = $words[2];
	    $betaconstraints[$k] = $words[3];

	    delete($remainingBETAs{$$biophysFiles[$i], $j});
	    $nBETA++;
	  }
	}
      }
      if($intercept && defined($remainingBETAs{'intercept', 0})) {
	my $k = $remainingBETAs{'intercept', 0};

	$betapriors[$k] = $words[2];
	$betaconstraints[$k] = $words[3];

	delete($remainingBETAs{'intercept', 0});
	$nBETA++;
      }
    }
    else {
      my @betas = split(/,/, $words[0]);
      foreach my $beta (@betas) {
	my $ix;
	if($mode == 1) {
	  $ix = 1;
	}
	elsif($beta =~ /^(.*)^M$/) { # gauss -- mean (mode ix == 1)
	  $beta = $1;
	  $ix = 1;
	}
	elsif($beta =~ /^(.*)^S$/) { # gauss -- st.dev (mode ix == 2)
	  $beta = $1;
	  $ix = 2;
	}
	elsif($beta =~ /^(.*)^(\d+)$/) { # poly -- power
	  $beta = $1;
	  $ix = $2;
	}
	elsif($beta =~ /^(.*)^S(\d+)$/) { # spline -- number
	  $beta = $1;
	  $ix = $2;
	}
	elsif($beta eq 'intercept' && $intercept) {
	  $ix = 0;
	}
	else {
	  die "Error in GA file $gafile: expecting <BETA>^",
	    (($mode == 2 && $intercept == 0) ? "M|S" : "<n>"),
	      ", found $beta\n";
	}
	if(defined($remainingBETAs{$beta, $ix})) {
	  my $i = $remainingBETAs{$beta, $ix};
	  $betapriors[$i] = $words[2];
	  $betaconstraints[$i] = $words[3];
	  delete($remainingBETAs{$beta, $ix});
	  $nBETA++;
	}
	else {
	  die "Error in GA file $gafile: unrecognised biophysical property ",
	  "file name $beta, or priors and constraints already defined for ",
	  "this beta\n";
	}
      }
    }
  }

  if($line ne '}') {
    die "Error in GA file $gafile: expecting \"}\", found EOF\n";
  }

  if($nBETA != $mode * scalar(@$biophysFiles) + $intercept) {
    my $files = join(", ", keys(%remainingBETAs));
    die "Error in GA file $gafile: priors and constraints not defined for ",
    "the following betas: $files\n";
  }

  push(@$priors, @betapriors);
  push(@$constraints, @betaconstraints);
}

# waitforruns(identifers)
#
# Wait for all runs to terminate. In case it helps, a list of identifiers
# is passed as argument

sub waitforruns {
  my ($identifiers) = @_;

  &wait(1);
}

# getdatetime()
#
# Return the date and time as a string to use in an identifier

sub getdatetime() {
  my ($s, $mi, $h, $d, $mo, $y) = localtime(time());

  $mo++;
  $y %= 100;
  
  return sprintf("%02d%02d%02d%02d%02d%02d", $y, $mo, $d, $h, $mi, $s);
}

{
  my %children;			# Local variable for startrun and wait
  my %exitstatus;		# Local variable for wait and runfitness

  # startrun(spommFile, search, biophysFiles, sample) -> run ID
  # 
  # Fork to start a run, waiting until there's a spare slot

  sub startrun {
    my ($spommFile, $search, $biophysFiles, $sample) = @_;

    &wait($maxproc);
  FORK:
    {
      my $dt = &getdatetime();
      my $pid;
      if($pid = fork()) {	# parent
	&log("Sample", &popstr($sample, $search, $biophysFiles),
	     "process $pid stored in directory $wd/$pid-$dt");
	$children{$pid} = "$wd/$dt-$pid";
	return "$dt-$pid";
      }
      elsif(defined($pid)) {	# child
	&run("$wd/$dt-$$", $spommFile, $search, $biophysFiles, $sample);
	exit 1;			# shouldn't get here
      }
      elsif($! == EAGAIN) {
	sleep 5;
	redo FORK;
      }
      else {
	die "Can't fork: $!\n";
      }
    }
  }
  
  # wait(max_nchildren)
  #
  # Wait until the number of child processes is less than the argument

  sub wait {
    my ($max_nchildren) = @_;
    
    while(scalar(keys(%children)) >= $max_nchildren) {
      my $pid = wait();
      if($pid != -1 && defined($children{$pid})) {
	$exitstatus{$children{$pid}} = $?;
	&log("Run $pid stopped with exit status $?");
	delete $children{$pid};
      }
      elsif($pid == -1) {
	die "Expecting ", scalar(keys(%children)), " child processes, but ",
	"there don't seem to be any\n";
      }
      elsif(!defined($children{$pid})) {
	warn "Child process $pid is not one I knew about!\n";
      }
    }
  }
  
  # runfitness(pid, occupFile)
  # 
  # Compute the run fitness. Contrary to what was promised earlier, the
  # occupancy file is assumed to take the format generated by
  # setup-speciesdist.pl.

  sub runfitness {
    my ($id, $occupFile) = @_;

    my @fitness = (0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
    shift(@fitness) if $fitnessnooccup;
    shift(@fitness) unless $nohalf;

    if($exitstatus{$id} != 0) {
				# Run failed: small fitness
      delete $exitstatus{$id};
      return \@fitness;
    }
    delete $exitstatus{$id};

    my $rundir = "$wd/$id";

    if(-e "$rundir/InvalidBetaVector") {
				# Beta vector produced habitat out of range
      open(ZERO, "<$rundir/InvalidBetaVector");
      my @lines = <ZERO>;
      close(ZERO);
      my $nvalid = $lines[$#lines - 1];
      my $npatches = $lines[$#lines];
      chomp $nvalid;
      chomp $npatches;
      
      &rm($rundir) if $delete;

      $fitness[0] = $nvalid + 1;
				# Penultimate line of this file should have the
				# number of valid patches on it
      return \@fitness;
    }

    my %spommparam = &readspomm("$rundir/spomm-$id.spom");
    $fitness[0] = $spommparam{'nPatches'} + 2;

    # Read the observed file

    my @correct;
    my @x;
    my @y;
    open(OCCUP, "<$occupFile")
      or die "Cannot open base occupancy file $occupFile: $!\n";
    <OCCUP>;
    while(my $line = <OCCUP>) {
      chomp $line;
      my @words = split(" ", $line);
      push(@correct, $words[$#words]);
      if(scalar(@words) == 3) {
	push(@x, $words[0]);
	push(@y, $words[1]);
      }
    }
    close(OCCUP);

    # Read the results file

    my $startstep = int(($spommparam{'nStep'} - 1) / 2);
    my $midstep = $startstep + int($startstep / 2);

    my $psppfile = "$rundir/occupancy-$id.csv";

    my @results;
    my @occupancy;
    my @q3results;
    my @q4results;

    open(PSPP, "<$psppfile")
      or die "Cannot open result occupancy file $psppfile: $!\n";
    <PSPP>;

    my $n = 0;
    my $q3n = 0;
    my $q4n = 0;
    my $last_step = -1;
    while(my $line = <PSPP>) {
      chomp $line;
      my ($step, $patch, $x, $y, $occup) = split(/,/, $line);
      if($step == $startstep) {
	$results[$patch - 1] = $occup;
	$n++ if($step != $last_step);
	$q3results[$patch - 1] = $occup;
	$q3n++ if($step != $last_step);
      }
      elsif($step > $startstep) {
	$results[$patch - 1] += $occup;
	$n++ if($step != $last_step);

	if($step == $midstep) {
	  $q4results[$patch - 1] = $occup;
	  $q4n++ if($step != $last_step);
	}
	elsif($step < $midstep) {
	  $q3results[$patch - 1] += $occup;
	  $q3n++ if($step != $last_step);
	}
	elsif($step > $midstep) {
	  $q4results[$patch - 1] += $occup;
	  $q4n++ if($step != $last_step);
	}
      }

      if($step != $last_step) {
	$occupancy[$step] = $occup;
	$last_step = $step;
      }
      else {
	$occupancy[$step] += $occup;
      }
    }
    close(PSPP);

    if(scalar(@correct) != scalar(@results)) {
      die "Base occupancy file $occupFile contains ", scalar(@correct),
      " patches, whilst results file $psppfile contains ", scalar(@results),
      " patches\n";
    }

    # Check for extinction

    my $extinctionFile = $spommparam{'exctinctionFile'};
    open(EXTINCT, "<$rundir/$extinctionFile")
      or die "Cannot open $rundir/$extinctionFile: $!\n";
    my @lines = grep(/\S/, <EXTINCT>);
    close(EXTINCT);

    if(scalar(@lines) > 1) {
      my $line = chomp $lines[$#lines];
      my ($spp, $year) = split(/,/, $line);

      $fitness[1] = $year + 1;
      return \@fitness unless $useextinct;
    }
    else {
      $fitness[1] = $spommparam{'nStep'} + 1;
    }

    # Check for equilibrium. This is done by dividing the latter half of
    # the simulation into two parts, and testing the difference between
    # the distributions of occupancy at the 5% level.

    my $mean1 = 0;
    my $var1 = 0;
    my $mean2 = 0;
    my $var2 = 0;
    my $m = 0;
    for(my $i = $startstep, my $j = $midstep; $i < $midstep; $i++, $j++) {
      $mean1 += $occupancy[$i];
      $mean2 += $occupancy[$j];
      $m++;
    }
    $mean1 /= $m;
    $mean2 /= $m;

    for(my $i = $startstep, my $j = $midstep; $i < $midstep; $i++, $j++) {
      my $diff1 = $occupancy[$i] - $mean1;
      my $diff2 = $occupancy[$j] - $mean2;
      $var1 += $diff1 * $diff1;
      $var2 += $diff2 * $diff2;
    }
    $var1 /= $m - 1;
    $var2 /= $m - 1;

    if($var1 > 0 && $var2 > 0) {
				# If both variances are zero (e.g. when there
				# is total extinction), then equilibirum has
				# been reached.
      my $z = ($mean1 - $mean2) / sqrt(($var1 / $m) + ($var2 / $m));

      if(abs($z) > 1.96) {

	$fitness[2] = 1 + ($m / abs($z));
				# Not reached equilibrium. There's an
				# argument for saying we should increase
				# nStep and rerun here, as in the general
				# case, we won't know how long it takes to
				# reach equilibrium with the right parameters.
	return \@fitness unless $usenoequi;
      }
      else {
	$fitness[2] = 2 + ($m / 1.96);
      }
    }
    else {
      $fitness[2] = 2 + ($m / 1.96);
    }

    # Check for patch-level equilibrium

    $fitness[3] = 1;
    for(my $i = 0; $i <= $#q3results; $i++) {
      my $q3np = $q3results[$i];
      my $q3nq = $q3n - $q3np;
      my $q4np = $q4results[$i];
      my $q4nq = $q4n - $q4np;

      if(($q3np == 0 || $q3nq == 0) && ($q4np == 0 || $q4nq == 0)) {
	$fitness[3]++ if($q3np == $q4np || $q3nq == $q4nq);
				# One of these two will be true if the same
				# equilibrium is reached even if q3n != q4n
      }
      elsif($q3np >= 5 && $q3nq >= 5 && $q4np >= 5 && $q4nq >= 5) {
				# According to my stats book (p. 207), the
				# normal approximation to the binomial is
				# good in this case
	my $q3npq = ($q3np * $q3nq) / $q3n;
	my $q4npq = ($q4np * $q4nq) / $q4n;

	my $z = ($q3np - $q4nq) / sqrt(($q3npq / $q3n) + ($q4npq / $q4n));
				# Use De Moivre's approximation to test the
				# difference between the two distributions
	$fitness[3]++ if(abs($z) <= 1.96);
      }
      else {			# Do it the hard way...
	my $P = 0;		# Compute the probability they are the same
	my $q3p = $q3np / $q3n;
	my $q3q = 1 - $q3p;
	my $q4p = $q4np / $q4n;
	my $q4q = 1 - $q4p;
	for(my $i = 0; $i <= $q3n && $i <= $q4n; $i++) {
	  my $q3P = &combinations($q3n, $i);
	  my $q4P = ($q3n == $q4n) ? $q3P : &combinations($q4n, $i);
	  
	  $q3P *= ($q3p ** $i) * ($q3q ** ($q3n - $i));
	  $q4P *= ($q4p ** $i) * ($q4q ** ($q4n - $i));

	  $P += $q3P * $q4P;
	}
	$fitness[3]++ if($P >= 0.05);
      }
    }

    # Compute the fitness

    open(SPPDIST, ">$rundir/sppdist-$id.txt")
      or die "Cannot create $rundir/sppdist-$id.txt: $!\n";

    if(scalar(@x) > 0) {
      print SPPDIST "x y p.occupied.$n\n";
    }
    else {
      print SPPDIST "patch.id p.occupied.$n\n";
    }

    my $matchfitness = 0;
    my $tpoccup = 0;
    my $toccup = 0;
    my $halfdiffoccup = 0;
    for(my $i = 0; $i <= $#correct; $i++) {
      my $poccup = $results[$i] / $n;
      my $error = $correct[$i] - ($poccup);
      $matchfitness += 1 - ($error * $error);

      $tpoccup += $poccup;
      $toccup += $correct[$i];
      $halfdiffoccup += 2 * abs($poccup - 0.5);

      if(scalar(@x) > 0) {
	print SPPDIST "$x[$i] $y[$i] $poccup\n";
      }
      else {
	print SPPDIST ($i + 1), " $poccup\n";
      }
    }
    close(SPPDIST);

    $fitness[4] = $matchfitness + 1;

    my $occupfitness = $spommparam{'nPatches'} - abs($toccup - $tpoccup);

    $fitness[5] = $occupfitness + 1 unless $fitnessnooccup;

    if($nohalf) {
      $fitness[$fitnessnooccup ? 5 : 6] = $halfdiffoccup + 1;
    }

    # Remove the directory if required, otherwise prepare a 'blob' plot

    if($delete) {
      &rm($rundir);
    }
    elsif($doblob) {
      if(!-x "$blob") {
	$doblob = 0;
	&log("Can't find blob.R at $blob -- stopping blob creation");
      }
      else {
	my $sppdist = "$rundir/sppdist-$id.txt";
	my $blobfile = "$rundir/sppdist-$id.pdf";
	if(system("$blob $sppdist $blobkey $blobfile") != 0) {
	  $doblob = 0;
	  &log("Running $blob failed ($!) or had nonzero exit status ($?)",
	       "-- stopping blob creation");
	}
      }
    }

    return \@fitness;
  }
}

# run(rundir, spommFile, search, biophysFiles, sample)
#
# Prepare the files in the working directory and run the SPOMM

sub run {
  my ($rundir, $spommFile, $search, $biophysFiles, $sample) = @_;

  my $spomm = &preparefiles($rundir, $spommFile, $search,
			    $biophysFiles, $sample);
  chdir $rundir or die "Cannot cd to $rundir: $!";
  exec "$cmd -s -b -p $spomm > stdout-$$.txt 2> stderr-$$.txt";
  die "exec $cmd -s -b -p $spomm from $wd failed: $!\n";
}

# preparefiles(rundir, spommFile, search, biophysFiles, sample
#
# rundir (scalar): directory to run the experiment from, not created.
#
# spommFile (scalar): SPOMM parameter file on which to base the runs, all
#        files referred to therein should be in the current working directory
#
# search (array): array of SPOMM parameters to search
#
# biophysFiles (array): array of files containing biophysical data
#
# sample (array): array of parameters in the gene, the length of which
#        should be the same as the sum of the lengths of the search and
#        biophysFiles arrays.

sub preparefiles {
  my ($rundir, $spommFile, $search, $biophysFiles, $sample) = @_;

  if(!-e "$rundir") {
    mkdir($rundir) || die "Cannot create directory $rundir: $!\n";
  }
  else {
    # This actually isn't very helpful unless the new directory
    # is returned as a different identifier for the population
    # member somehow. Hopefully this will not happen now that
    # the date and time are included in the directory name.
    my $kk = 0;
    while(-e "$rundir-$kk") {
      $kk++;
    }
    &mv($rundir, "$rundir-$kk");
    mkdir($rundir) || die "Cannot create directory $rundir: $!\n";
  }

  my @spommsample;
  my @betasample;
  for(my $i = 0; $i <= $#$sample; $i++) {
    if($i < scalar(@$search)) {
      push(@spommsample, $$sample[$i]);
    }
    else {
      push(@betasample, $$sample[$i]);
    }
  }

  open(GENE, ">", "$rundir/gene.txt")
    or die "Cannot create gene file gene.txt in directory $rundir: $!\n"

  my %spommparam = &readspomm($spommFile);
  my %sppparam;
  for(my $i = 0; $i <= $#$search; $i++) {
    my ($file, $parameter) = split("/", $$search[$i]);

    print GENE "$$search[$i] = $spommsample[$i]\n";
    
    if($file eq 'SPOMM') {
      $spommparam{$parameter} = $spommsample[$i];
    }
    elsif($file eq 'spp' || $file eq 'species') {
      $sppparam{$parameter} = $spommsample[$i];
    }
    else {
      die "Unrecognised SPOMM parameter file in GA file: $file. Expecting ",
      "SPOMM or species\n";
    }
  }
  close(GENE);

  $spommparam{'speciesPerPatchFile'} = "occupancy-$$.csv";
  $spommparam{'nStepListSpecies'} = 1;

  my $patchFile = &buildpatches($rundir, $biophysFiles,
				$spommparam{'nPatches'}, \@betasample);
  $spommparam{'patchFile'} = $patchFile;

  my $sppFile = &buildspecies($rundir, $spommparam{'speciesFile'}, \%sppparam);
  $spommparam{'speciesFile'} = $sppFile;

  return &buildspomm($rundir, \%spommparam);
}

# readspomm(spommFile)
#
# Read the spomm file and return the parameter values as a hash

sub readspomm {
  my ($spommFile) = @_;

  open(SPOMM, "<$spommFile") or die "Cannot open spomm file $spommFile: $!\n";

  my %spommparam;

  my $beginFound = 0;
  while(my $line = <SPOMM>) {
    chomp $line;
    if($line eq '@begin') {
      $beginFound = 1;
      next;
    }
    elsif(!$beginFound) {
      if($line !~ /^\#/ && $line !~ /^\s*$/) {
	die "Invalid SPOMM file $spommFile: expecting \"\@begin\" found ",
	"$line\n";
      }
      next;
    }

    if($line eq '@end') {
      last;
    }

    my ($param, $value) = split(" ", $line);

    $spommparam{$param} = $value;
  }
  close(SPOMM);

  return %spommparam;
}

# buildpatches(rundir, biophysFiles, npatches, betasample)
#
# Build the patches file. This doesn't, as promised earlier, do anything
# sophisticated to guess the biophysical file format, but assumes a two
# column csv file with patch id first and data value second. For splines,
# remaining columns are assumed to contain the spline variables.

sub buildpatches {
  my ($rundir, $biophysFiles, $npatches, $betasample) = @_;

  my $n_expected_files = 0;


  open(GENE, ">>", "$rundir/gene.txt")
    or die "Cannot append to gene file gene.txt in directory $rundir: $!\n"

  if($betaMode eq "linear") {
    $n_expected_files = scalar(@$betasample) - 1;
    for(my $i = 0; $i <= $#$biophysFiles; $i++) {
      print GENE "BETA/$$biophysFiles[$i] = $$betasample[$i]\n";
    }
    print GENE "BETA/intercept = $$betasample[$#$betasample]\n";
  }
  elsif($betaMode eq "Gaussian") {
    $n_expected_files = scalar(@$betasample) / 2;
    for(my $i = 0; $i <= $#$biophysFiles; $i++) {
      print GENE "MU/$$biophysFiles[$i] = $$betasample[$i]\n";
      print GENE "SIGMA/$$biophysFiles[$i] = ",
	"$$betasample[$n_expected_files + $i]\n";
    }
  }
  elsif($betaMode =~ /^polynomial,(\d+)$/) {
    my $max = $1;
    $n_expected_files = (scalar(@$betasample) - 1) / $max;
    for(my $power = 1; $power <= $max; $power++) {
      for(my $var = 0; $var < $n_expected_files; $var++) {
	my $k = (($power - 1) * $n_expected_files) + $var;

	print GENE "BETA^$power/$$biophysFiles[$var] = $$betasample[$k]\n";
      }
    }
    print GENE "BETA/intercept = $$betasample[$#$betasample]\n";
  }
  elsif($betaMode =~ /^spline,(\d+)$/) {
    my $nsplines = $1;
    $n_expected_files = scalar(@$betasample) / $nsplines;
    for(my $i = 0; $i < $nvars; $i++) {
      for(my $j = 1; $j <= $nsplines; $j++) {
	print GENE "SPLINE_$j/$$biophysFiles[$i] = ",
	  $$betasample[(($j - 1) * $nvars) + $i], "\n";
      }
    }
  }
  else {
    die "Unrecognised beta mode $betaMode\n";
  }

  close(GENE);

  if(scalar(@$biophysFiles) != $n_expected_files) {
    die "Building patches from different number of biophysical files (",
    scalar(@$biophysFiles), ") than expected ($n_expected_files) ",
    "given beta mode ($betaMode) and beta vector length (",
    scalar(@$betasample), ")\n";
  }

  my @biophys;
  for(my $i = 0; $i <= $#$biophysFiles; $i++) {
    open(BIOPHYS, "<$$biophysFiles[$i]")
      or die "Cannot open biophysical file $$biophysFiles[$i]: $!\n";
    <BIOPHYS>;			# Chuck away the first line
    my $patchcount = 0;
    my @data;
    while(my $line = <BIOPHYS>) {
      chomp $line;
      last if $line =~ /^\s*$/;
      $patchcount++;
      my ($patchid, @value) = split(/,/, $line);
      push(@data, \@value);
    }
    close(BIOPHYS);
    if($patchcount != $npatches) {
      die "Error in biophysical file $$biophysFiles[$i]: expecting $npatches ",
      "of patch data, but got $patchcount\n";
    }
    push(@biophys, \@data);
  }

  my $patchFile = "$rundir/patches-$$.csv";

  open(PATCHES, ">$patchFile")
    or die "Cannot create patch file $patchFile: $!\n";

  print PATCHES "Patch number,Habitat,Species\n";

  my $ninvalid = 0;
  for(my $i = 0; $i < $npatches; $i++) {
    my $habitat;

    if($betaMode eq 'linear') {
      $habitat = $$betasample[$#$betasample]; # intercept
      for(my $j = 0; $j < $#$betasample; $j++) {
	$habitat += $biophys[$j][$i][0] * $$betasample[$j];
      }
    }
    elsif($betaMode eq 'Gaussian') {
      $habitat = 1;
      for(my $j = 0; $j < scalar(@$betasample) / 2; $j++) {
	$habitat *= exp(-1.0 * ((($biophys[$j][$i][0] - $$betasample[$j]) ** 2)
			/ $$betasample[$j + (scalar(@$betasample) / 2)] ** 2));
      }
    }
    elsif($betaMode =~ /^polynomial,(\d+)$/) {
      my $max = $1;
      $habitat = $$betasample[$#$betasample]; # intercept
      my $nvars = $#$betasample / $max;

      for(my $power = 1; $power <= $max; $power++) {
	for(my $var = 0; $var < $nvars; $var++) {
	  my $k = (($power - 1) * $nvars) + $var;

	  $habitat += ($biophys[$var][$i][0] ** $power) * $$betasample[$k];
	}
      }
    }
    elsif($betaMode =~ /^spline,(\d+)$/) {
      my $nspline = $1;
      $habitat = 0;
      my $nvars = scalar(@$betasample) / $nspline;
      for(my $j = 0; $j < $nvars; $j++) {
	for(my $k = 1; $k <= $nspline; $k++) {
	  # 20140516: I don't understand why $k goes from 1-$nspline
	  # instead of 0 to $nspline - 1. biophys[j][i][0] is never
	  # going to be used... this must depend somehow on the format
	  # of the biophys file
	  $habitat += $biophys[$j][$i][$k]
	    * $$betasample[(($k - 1) * $nvars) + $j];
	  # 20140516: (k - 1) changed from original k. Again, I can't
	  # think why I would have had it as k, which would mean the
	  # first nvars of betasample are never used.
	}
      }
    }
    else {
      die "Unrecognised beta mode $betaMode (should not get here!)\n";
    }

    if($habmodifier eq 'sigmoid' || $betaMode =~ /^spline,/) {
      $habitat = 1 / (1 + exp(-$habitat));

      if($betaMode =~ /^spline,/
	 && ($habmodifier ne 'none' && $habmodifier ne 'sigmoid')) {
	die "Habitat modifier (argument to -habmod) $habmodifier is invalid ",
	"for beta mode spline (sigmoid selected by default)\n";
      }
    }
    elsif($habmodifier eq 'pct') {
      $habitat /= 100;
    }
    elsif($habmodifier ne 'none') {
      die "Unrecognised habitat modifier (argument to -habmod): ",
      "$habmodifier\n";
    }

    if($habitat < 0 || $habitat > 1) {
      open(ZERO, ">>$rundir/InvalidBetaVector")
	or die "Cannot write to $rundir/InvalidBetaVector\n";
      print ZERO "Patch $i, habitat = $habitat\n";
      close(ZERO);
  
      $ninvalid += ($habitat < 0 ? exp($habitat) : exp(1 - $habitat));
				# Tell the GA how close we were to having a
				# valid habitat (So, this isn't really
				# counting the number of invalid patches.)
    }
    print PATCHES $i + 1, ",$habitat,", ($habitat > 0 ? 1 : 0), "\n";
  }

  close(PATCHES);

  if($ninvalid > 0) {
    open(ZERO, ">>$rundir/InvalidBetaVector")
      or die "Cannot write to $rundir/InvalidBetaVector\n";
    print ZERO "$ninvalid\n$npatches\n";
    close(ZERO);
    exit 0;
  }

  return $patchFile;
}

# buildspecies(rundir, speciesFile, param)
#
# Build the species file for the run from the original species file using
# the species parameters

sub buildspecies {
  my ($rundir, $speciesFile, $param) = @_;

  open(SPPIN, "<$speciesFile")
    or die "Cannot open species file $speciesFile: $!\n";

  my $headers = <SPPIN>;
  my $values = <SPPIN>;

  close(SPPIN);

  if(!$headers || !$values) {
    die "Unexpected EOF in file $speciesFile\n";
  }

  chomp $headers;
  chomp $values;

  my @keys = split(/,/, $headers);
  my @values = split(/,/, $values);
  my %sppparam;

  for(my $i = 0; $i <= $#keys; $i++) {
    if(defined($$param{$keys[$i]})) {
      $sppparam{$keys[$i]} = $$param{$keys[$i]};
    }
    else {
      $sppparam{$keys[$i]} = $values[$i];
    }
  }

  my $runsppFile = "$rundir/species-$$.csv";

  open(SPPOUT, ">$runsppFile")
    or die "Cannot create run species file $runsppFile: $!\n";

  print SPPOUT "$headers\n";
  for(my $i = 0; $i <= $#keys; $i++) {
    print SPPOUT "," if $i > 0;
    print SPPOUT "$sppparam{$keys[$i]}";
  }
  print SPPOUT "\n";
  close(SPPOUT);

  return $runsppFile;
}

# buildspomm(rundir, param)
#
# Build the spomm file

sub buildspomm {
  my ($rundir, $param) = @_;

  my $spommFile = "$rundir/spomm-$$.spom";
  open(SPOMM, ">$spommFile") 
    or die "Cannot create SPOMM file $spommFile: $!\n";

  print SPOMM '@begin', "\n";
  while(my ($key, $value) = each(%$param)) {
    if($key eq 'changeHabitatFile'
       || $key eq 'autoCorrelatedFieldFile'
       || $key eq 'goodBadYearFile'
       || $key eq 'csvLocalizedField_file'
       || $key eq 'landUseHabitatFile'
       || $key eq 'predationFile'
       || $key eq 'habitatSpecificMuFile'
       || $key eq 'sinkHabitatPropertieFile'
       || $key eq 'occupiedPatchesPerSpeciesOutputFile') {
      if(!-e "$value") {
	open(TOUCH, ">$rundir/$value")
	  or die "Cannot touch $rundir/$value: $!\n";
	close(TOUCH);
      }
      else {
	open(READ, "<$value") or die "Cannot open $value: $!\n";
	open(WRITE, ">$rundir/$value")
	  or die "Cannot write $rundir/$value: $!\n";
	print WRITE <READ>;
	close(WRITE);
	close(READ);
      }
    }
    print SPOMM "$key $value\n";
  }
  print SPOMM '@end', "\n";

  close(SPOMM);

  return $spommFile;
}

# mv(from, to)
#
# Moves one directory to another

sub mv {
  my ($from, $to) = @_;

  mkdir($to) or die "Cannot create directory $to: $!\n";

  opendir(DIR, "$from") or die "Cannot read directory $from: $!\n";
  foreach my $file (readdir(DIR)) {
    next if($file eq '.' || $file eq '..');
    link("$from/$file", "$to/$file")
      or die "Cannot link $from/$file to $to/$file: $!\n";
    unlink("$from/$file") or die "Cannot unlink $from/$file: $!\n";
  }
  closedir(DIR);
  
  rmdir($from) or die "Cannot remove directory $from: $!\n";

  &log("Moved $from to $to");
}

# cp(from, to)
#
# Copies one directory to another

sub cp {
  my ($from, $to) = @_;

  mkdir($to) or die "Cannot create directory $to: $!\n";

  opendir(DIR, "$from") or die "Cannot read directory $from: $!\n";
  foreach my $file (readdir(DIR)) {
    next if($file eq '.' || $file eq '..');
    open(FROM, "<$from/$file") or die "Cannot read file $from/$file: $!\n";
    open(TO, ">$to/$file") or die "Cannot create file $to/$file: $!\n";
    print TO <FROM>;
    close(TO);
    close(FROM);
  }
  closedir(DIR);

  &log("Copied $from to $to");
}

# rm(dir)
#
# Removes a directory

sub rm {
  my ($dir) = @_;

  opendir(DIR, "$dir") or die "Cannot read directory $dir: $!\n";
  foreach my $file (readdir(DIR)) {
    next if($file eq '.' || $file eq '..');
    unlink("$dir/$file") or die "Cannot unlink $dir/$file: $!\n";
  }
  closedir(DIR);

  rmdir($dir) or die "Cannot remove directory $dir: $!\n";

  &log("Removed $dir");
}

# log(data)
#
# Logs information

sub log {
  my @data = @_;

  my ($sec, $min, $hour, $mday, $mon, $year) = localtime;

  open(LOG, ">>$log") or die "Cannot append to log file $log: $!\n";
  print LOG "SPOMM-MGAv3 ($$) ", $year + 1900,
  sprintf("%02d%02d", $mon + 1, $mday), "T",
  sprintf("%02d%02d%02d", $hour, $min, $sec), ": ", join(" ", @data), "\n";
  close(LOG);
}
