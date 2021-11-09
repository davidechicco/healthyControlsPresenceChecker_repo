#!/usr/bin/env perl

#import standard Perl modules
use strict;
use warnings;
use Term::ANSIColor;
use Getopt::Std;
use LWP::Simple;
use File::Basename;
use File::HomeDir;

# We need a random number for the temporary file name
srand(42);
my $randomNumber = int rand(1000);

# my $GEOcode = "GSE30174";
my $GEOcode = "GSE38043";

print "> > > lynx phase started\n";

# lynx part
my $lynxOutputFile="temp_lynx_output_${GEOcode}_rand${randomNumber}.txt";
 print "Created file $lynxOutputFile\n";

my $GEOurl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GEOcode}";
print "Executing lynx --dump $GEOurl > $lynxOutputFile\n";
system("lynx --dump $GEOurl > $lynxOutputFile");

my $output_lynx = `grep "healthy control" $lynxOutputFile`;
my $result_lynx_search = "";
if ($output_lynx eq "")
{
        $result_lynx_search = "FALSE";
}
else
{
        $result_lynx_search = "TRUE";
}

print "Outcome of healthy controls search in the file downloaded by lynx for $GEOcode: $result_lynx_search \n";

# Remove the temporary file created earlier
system("rm $lynxOutputFile;");
print "Removed file $lynxOutputFile\n";
print "< < <  lynx phase finished\n\n\n";

print "> > > healthyControlsPresence.r phase started\n";


if($result_lynx_search eq "FALSE"){

    my $RscriptFile = './healthyControlsPresentInputParams.r';
    unless (-e $RscriptFile) {

        # The script downloads the R script if it does not exist in the current folder
        print "\n The \"$RscriptFile\" file doesn't exist and will be downloaded now\n";
        system `curl https://raw.githubusercontent.com/davidechicco/healthyControlsPresenceChecker_repo/main/healthyControlsPresentInputParams.r -o healthyControlsPresentInputParams.r`;

    }

    my $rScriptOutputFile="temp_R_output_${GEOcode}_rand${randomNumber}.txt";
    print "Created file $rScriptOutputFile\n";
    
    # Call to the R script that checks the presence of healthy controls and prints its output in the $rScriptOutputFile
    system("Rscript healthyControlsPresentInputParams.r $GEOcode TRUE > $rScriptOutputFile");

    # Call to awk in bash to read the last word of the last line of the temporary file, that can be TRUE or FALSE
    my $outcome_rscript = `awk 'END {print \$NF}' $rScriptOutputFile`;
    $outcome_rscript =~ s/\n//;

    # Print the final outcome: TRUE if healhty controls were found, FALSE otherwise.
    print "Outcome of healthyControlsPresent.r for $GEOcode: $outcome_rscript\n";

    # print "outcome_rscript=$outcome_rscript\n";

    # Remove the temporary file created earlier
    system("rm $rScriptOutputFile;");
    print "Removed file $rScriptOutputFile\n";
    print "< < < healthyControlsPresence.r phase finished\n\n\n";



}



