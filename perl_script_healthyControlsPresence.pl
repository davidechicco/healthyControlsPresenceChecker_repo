#!/usr/bin/env perl

#import standard Perl modules
use strict;
use warnings;
use Term::ANSIColor;
use Getopt::Std;
use LWP::Simple;
use File::Basename;
use File::HomeDir;


my $RscriptFile = './healthyControlsPresentInputParams.r';
unless (-e $RscriptFile) {

    # The script downloads the R script if it does not exist in the current folder
    print "\n The \"$RscriptFile\" file doesn't exist and will be downloaded now\n";
    system `curl https://raw.githubusercontent.com/davidechicco/healthyControlsPresenceChecker_repo/main/healthyControlsPresentInputParams.r -o healthyControlsPresentInputParams.r`;

}

# We need a random number for the temporary file name
srand(42);
my $randomNumber = int rand(1000);

my $GEOcode = "GSE30174";
my $rScriptOutputFile="temp_R_output_rand$randomNumber.txt";
print "Created file $rScriptOutputFile\n";

# Call to the R script that checks the presence of healthy controls and prints its output in the $rScriptOutputFile
system("Rscript healthyControlsPresentInputParams.r $GEOcode TRUE >> $rScriptOutputFile");

# Call to awk in bash to read the last word of the last line of the temporary file, that can be TRUE or FALSE
my $outcome = `awk 'END {print \$NF}' $rScriptOutputFile`;

# Print the final outcome: TRUE if healhty controls were found, FALSE otherwise.
print "Outcome of healthyControlsPresent.r for $GEOcode: $outcome";

# Remove the temporary file created earlier
system("rm $rScriptOutputFile;");
print "Removed file $rScriptOutputFile\n";


