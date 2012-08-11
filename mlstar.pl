#!/usr/bin/perl
#
# Multi Loci Sequence Typing Analysis helpeR
#
# This script matches and returns an ST type number for a given base reference sequence, MLST fragment files (all fasta 
# formatted), and ST type table file.
#
# Usage: mlstar.pl -c <config file> [-quiet] [-debug] [-l <program log>] [-fast]
# 
# -c <config file = REQUIRED parameter with the name and path of the config file to use for this run. See below for details.
# -quiet = Optional flag to reduce amount of output to standard out.  See the log file for full output.
# -debug = Optional flag to generate additional debugging message (verbose).  Also stored in log file.
# -l <program log> = Optional parameter to specify the log file to write to during the program run.
# -fast = Tells mlstar to stop looking after finding the first matching allele of each gene fragment, and st type table.
# 
# By default, mlstar will look for matches in all alleles of all gene fragments and all entries in the st type table,
# throwing an error if duplicates are found.  Use -fast if you are sure you have no such duplicates in your data.
#
# The <config file> is a plain text file containing information in the following format (similar to a sectioned INI file):
#
# ------------------------------------------
# # Comments and blank lines are ignored
#
# [ST ALLELE TABLE FILES]
# st_allele_table.txt
#
# [ALLELE FILES]
# arcc_alleles.fasta
# aroe_alleles.fasta
#
# [REF FILES]
# NC007793_corrected.fsa
#
# -------------------------------------------
#
# In short, the config file has three sections, and the start of each section is indicated by a [SECTION] tag:
# [ST ALLELE TABLE FILES]
# [ALLELE FILES]
# [REF FILES]
#
# The non-blank, non-comment lines (beginning with #) after each section tag should be file names (and paths) 
# related to the given section.
#
# Files in the [ST ALLEL TABLE FILES] section should contain the path and filenames to the ST allele table, with the
# expected content format of:
# <ST_type> <arcc_allele_#> <aroe_allele_#> <glpf_allele_#> <gmk_allele_#> <pta_allele_#> <tpi_allele_#> <yqil_allele_#>
# (basically a table of numbers)
#
# Files in the [ALLELE FILES] section should contain FASTA formatted files with the allele sequences for the housekeeping genes.
#
# IMPORTANT NOTE: The top down order of the files in this section should match up against the left-to-right columnar order in the
# ST allele table file.  So if you have files listed as:
#
# arcc_alleles.fasta
# aroe_alleles.fasta
# gmk_allel.fasta
#
# in the allele files section, then the table file columns should be in the order: 
#
# ST_Type   arc_allele_#   aroe_allele_#   gmk_allel_# ...
# 
# Files in the [REF FILES] section contains the base reference sequence files of the isolate that you want to determine 
# the ST type of.
# 
# Milstar will process each base reference file in the [REF FILES] section, matching against the files listed in the [ALLELE FILES]
# section, and looking up the final ST type in the table files listed in the [ST ALLELE TABLE FILES] section.  If there is more
# than one file in the [ST ALLELE TABLE FILES] section, mlstar will use them in top-down order, stopping when it finds the first
# successful ST type lookup.
#
# ST analysis output is logged to REF_FILE.mlstar.out (i.e. one for each reference file).  Program output is logged to
# mlstar.log (or can be specified with the -l <logfile> parameter).
#
# 02/09/2012 AJP Initial development.

use strict;
use v5.10;

use constant VERSION => "1.1";

print "Mlstar (Multi Loci Sequence Typing Analysis helpeR)\nv".VERSION."\n\n";

if ($#ARGV < 1){
  print "Usage: $0 -c <config file> [-debug] [-quiet] [-l <log filename>] [-fast]\n";
	print "See script header comments on config file information.\n\n";
	exit 1;
}

my $debug="N";
my $quiet="N";
my $fast="N";
my $arg_num=0;
my $configfile;
my $configFilename;
my @referenceFiles;
my @alleleFiles;
my @sttableFiles;
my $logFilename="mlstar.log";
my $logfile;

my @alleleGlobalData;
my %sttableGlobalData;

my $noteCount = 0;
my $warningCount = 0;

#
# MAIN
#

process_parameters();
load_allele_data();
load_st_table_data();
process_reference_files();


#############
# FUNCTIONS
#############

sub get_reverse_complement{
	# Taken from Jeremiah Faith code listing in his article"
	# http://code.izzid.com/2011/08/25/How-to-reverse-complement-a-DNA-sequence-in-perl.html
	# Supports IUPAC nucleotide codes.

	my $dna = $_[0];

	# reverse the DNA sequence
  my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
  $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
  return $revcomp;
}

sub process_reference_files{
	# Processes each reference file in @referenceFiles looking for alleles defined in @alleleGlobalData, then
  # looks for the matching ST Type ID from %sttableGlobalData.
  #
  # The results for each reference file is written to <ref filename>.mlstar.out

  my $refFilename;

  foreach $refFilename (@referenceFiles){
		my $reffile;
		my $refoutfile;
		my $refoutFilename;
		my $refID;
		my @reffileLines;
		my $line=0;
		my $refsequence = "";
		my $refsttableRE = "^([0-9]+)"; # regular expression to match against the ST Type table data.
		my $alleleID = "";
		my $alleleidSeq = "";
		my $alleleName = "";
		my $partialSTLookup="N"; # If we don't find a matching allele, could end up matching multiple sttypes.  Report on that.

    print_all("Processing $refFilename...\n");

		open($reffile,"<",$refFilename) or die "Unable to read reference file $refFilename! $!";
		open($refoutfile,">",$refFilename.".mlstar.out") or die "Unable to write output file $refFilename.mlstar.out! $!";

    # Read in the base reference sequence, and join into one big long string.
		while (<$reffile>){
			$line++;
			chomp;
			# Save off the ID from the first line, which should be prefixed with >.
			if (/^>(.+)/){
				$refID=$1;
				print_debug("Found reference genome id $refID\n");
				next;
			}

      # The rest of the file should be base pair sequences
			if (/^[AGCT]+$/){
				push(@reffileLines,$_);
      }
			else{
				print_all("*** ERROR: Invalid format in reference file on line $line (invalid base pair character?): $_\n");
				exit 1;
			}
		}

    # Reconstruct the sequence as one big long line.
		$refsequence = join("",@reffileLines);
		#print_debug("Full ref sequence is $refsequence\n");

		print $refoutfile "Reference genome id: $refID\n\n";
		print_all("Found Reference genome id: $refID\n\n");

    # Look for the alleles in the ref sequence, for all the allele gene fragments defined in @alleleGlobalData.
    # as we find the matching allele for a gene fragment, build the search regex (in $refsttableRE) we will use to 
    # do the lookup in the sttableGlobalData hash table.

    # NOTE: The order in which we look for the alleles is the same order as listed in the allele files section of
    # the config file, which should match the left-to-right column order in the ST type table.

    my $i=0;
		while ($i <= $#alleleGlobalData){
      my %alleleHash = %{$alleleGlobalData[$i]};
			my $alleleNumber = "";
			my $alleleKey = "";
			my $alleleName = "";
			my $is_reversed = "N";
			my $foundAllele = "N";
			foreach $alleleKey (keys %alleleHash){
			  print_debug("Looking for $alleleKey in sequence.\n");

				# Get the allele id name and number -- assumed to be like foo123 (so name is foo and id is 123 in this made-up example).
				if ($alleleKey =~ m/(.+?)([0-9]+)$/){
					$alleleName=$1;
					$alleleNumber=$2;
					print_debug("Allele name appears to be $alleleName. Allele id number appears to be $alleleNumber.\n");
				}
				else{
					print_all("*** ERROR: Cannot determine allele id number from $alleleKey!\n");
					exit 1
				}

			  $is_reversed = "N";

        # Get the reverse complement sequence, since need to check for that too (gene fragments could be encoded in 5'-3' or 3'-5').
        my $revcomprefsequence = get_reverse_complement($alleleHash{$alleleKey});
				if (($refsequence =~ /$alleleHash{$alleleKey}/) || (($refsequence =~ /$revcomprefsequence/) && ($is_reversed="Y"))){
				  print_all("Found $alleleKey sequence at position $-[0] of reference.\n");
				  print $refoutfile "Found $alleleKey sequence at position $-[0] of reference.\n";

          if ($is_reversed eq "Y"){
            print_all("^^^NOTE: Match was found on the reverse complement of this allele sequence.\n\n");
            print $refoutfile "^^^NOTE: Match was found on the reverse complement of this allele sequence.\n\n";
						$noteCount++;
					}
					
					# See if we already found an allele for this gene fragment.
					if ($foundAllele ne "N"){
					  $alleleidSeq =~ /^(?: [0-9]+)* ([0-9])+$/;
					  my $lastMatchedAlleleNumber = $1;
					  print_all("*** ERROR: Reference base already has an allele of this gene.  Previously matched allele number is $alleleNumber\n");
					  print $refoutfile "*** ERROR: Reference base already has an allele of this gene.  Previously matched allele number is $alleleNumber\n";
						exit 1;

					}
					else{
						$alleleidSeq .= " $alleleNumber";
					}

          # Indicate we have now found an allele in this given file.
					$foundAllele = "Y";

          # Stop looking if -fast flag was set.
					if ($fast eq "Y"){
					  last;
				  }
        }
		  }

			# If no matches were found, record that fact, too.
			if ($foundAllele eq "N"){
			  print_all("\n*** WARNING: No matches were found for any of the alleles belonging to $alleleName!\n\n");
			  print $refoutfile "\n*** WARNING: No matches were found for any of the alleles belonging to $alleleName!\n\n";
				$warningCount++;
				$alleleidSeq .= " *";
				$partialSTLookup="Y";
			}
      $i++;
		}

		# Report the final allele id sequence
		print_all("\nFinal allele id sequence for this reference is $alleleidSeq\n");
		print $refoutfile "\nFinal allele id sequence for this reference sequence: $alleleidSeq\n";

    if ($partialSTLookup eq "Y"){
			print_all("\n*** WARNING: One or more gene fragment alleles were not found.  Mlstar will attempt to show all possible ST Types.\n");
			print $refoutfile "\n*** WARNING: One or more gene fragment alleles were not found.  Mlstar will attempt to show all possible ST Types.\n";
			$warningCount++;
		}

		# Do ST Type lookup here.
    print_all("\nLooking up ST Type identifier...\n");

		# Change the allele id sequence into a regular expression that can match against the st table file.
		$alleleidSeq =~ s/ /\\s+/g;
		$alleleidSeq =~ s/\*/\.\*\?/g;
		$refsttableRE .= $alleleidSeq;
		$refsttableRE .= "\\s\*\$"; #trailing spaces ok.
		print_debug("\nLooking up ST type with regex $refsttableRE.\n");

    my $stKey;
		my $foundKeys = 0;
		foreach $stKey (%sttableGlobalData){
      if ($sttableGlobalData{$stKey} =~ /$refsttableRE/){
			  # Error on finding multiple st types?
				if (($foundKeys != 0) && ($partialSTLookup eq "N")){
          print_all("*** ERROR: Multiple matches against st table data!  Matched $stKey when another had already been matched!\n");
          print $refoutfile "*** ERROR: Multiple matches against st table data!  Matched $stKey when another had already been matched!\n";
					exit 1;
				}
			  $foundKeys++;
        print_all("Found matching ST type, for ST type id $stKey.\n"); 
		    print $refoutfile "Matched ST type for this reference sequence: $stKey\n";

				# Stop looking if -fast flag was set. Ignore if looking for partial matches due to incomplete allele lookups.
				if (($fast eq "Y") && ($partialSTLookup eq "N")){
					last;
				}
			}
		}

    close($reffile);
	  close($refoutfile);
		print_all("\n");
	}

  if ($warningCount > 0){
		print_all("\nThere were $warningCount warnings for this run.\n");
	}

  if ($noteCount > 0){
		print_all("\nThere were $noteCount notes for this run.\n");
  }

	print_all("\nLog file for this run is in $logFilename.\n");
	print_all("\nDone.\n");
}

sub load_st_table_data{
  # Loads the ST type table file data into memory for fast lookups.
  # The table file columns representing the allele id should be in the same order (top to down) as listed in 
	# the alelle files section: 
  #
  # ST_Type   arc_allele_#   aroe_allele_#   gmk_allel_# ...
  #
  # We store this in memory in a hash table, keyed by the ST_Type ID with values being the entire row.
  # Duplicate ST_Type IDs are considered an error (i.e. like if the same ST_Type ID exists in the same or different
  # files but has different allele id values in the row).

	my $sttableFilename;
	my $sttablefile;
	my $stID;

  %sttableGlobalData = ();
  foreach $sttableFilename (@sttableFiles){
		print_all("Loading ST Type data from $sttableFilename.\n");
	  open($sttablefile,"<","$sttableFilename") or die "Unable to open allele file $sttableFilename for writing! $!\n";

    while (<$sttablefile>){
			chomp;

      # Ignore blank lines.
      if (/^ *$/){
				next;
      }

      # Look for expected format (should be a row of space delimited numbers, may or may not have trailing spaces).
			if (/^([0-9]+)(?:\s+[0-9]+)+\s*$/){
				$stID=$1;

				print_debug("Found ST Type ID $stID\n");
				
				# Make sure the ST Type ID doesn't already exist (and more importantly, does not have a different allele ID 
        # sequence defined.
        
				if (exists($sttableGlobalData{$stID})){
					print_all("*** WARNING: Duplicate ST Type ID found: $stID\n");
					print_all("The duplicated line is        : $_\n");
					print_all("The previously defined line is: $sttableGlobalData{$stID}\n");
					$warningCount++;

          # If the lines are the same, can ignore. First normalize away any extraneous whitespace between values.
					my $compare1=$_;
					my $compare2=$sttableGlobalData{$stID};
					$compare1 =~ s/\s+/ /g;
					$compare2 =~ s/\s+/ /g;
					if ($compare1 eq $compare2){
            print_all("ST Type definitions appear to be the same.  No conflict.  Continuing...\n");
					}
					else{
            print_all("*** ERROR: ST Type defintions do not match!  Conflict!  Quitting...\n");
						exit 1;
					}
			  }	
				else{
					print_debug("Adding new st type def $_\n");
          $sttableGlobalData{$stID}=$_;
				}
      }
		  else{
        print_all("*** ERROR: In file $sttableFilename, badly formatted line: $_ (expect just space delimited numbers)\n");
				exit 1;
			}
    }
		close($sttablefile);
	}
	print_all("\n");
	print "\n";
}

sub load_allele_data{
  # Loads all the allele information into the @alleleData array variable from the allele gene files defined in the cfg file.
  #
  # NOTE: The order in which the alleles are defined is important since it matches the same left-to-right
  # column order in the st_type table.  So, need to track this.
  # 
  # So the data structure we will use is:
  # @alleleGlobalData[order_num] = %alleleLocalData{"allele_name"} = $allele_sequence
  #
  # order_num is implicit (we are using push to add to the global data array).
  #
  # %alleleLocalData is a hash table keyed by allele name  and with the sequence as the value.
  #
  # $allele_sequence is a string value, but since the base sequence is stored line by line within the file, the file lines
  # are stored in a temporary array while being read, and then joined when being assigned to its hash key.  I.e.
  # there is no actual scalar $allele_sequence -- it is the join statement.
  #
  # In other words, it is an array of hash tables.  The array maintains the order of the allele files being loaded,
  # and each hash table at the array index holds the allele sequence variants, with the allele name as the key and
  # its base sequence as the value.

	my $alleleFilename;
	my $allelefile;
	my %alleleLocalData = ();
	my $currentAllele;
	my @currentAlleleData;


  foreach $alleleFilename (@alleleFiles){
		print_all("Loading allele data from $alleleFilename.\n");
	  open($allelefile,"<","$alleleFilename") or die "Unable to open allele file $alleleFilename for writing! $!\n";

	  @currentAlleleData = ();
		$currentAllele = "";
	  %alleleLocalData = ();

    while (<$allelefile>){
			chomp;

      # Ignore blank lines.
      if (/^ *$/){
				next;
      }

			# Look for line starting with >.  Should be the allele name, formatted alpha# with the # matching the value in
			# allele's corresponding column in the st_type table file.
			if (/^>(.+)/){

        # Join all the previous allele's data into one long sequence string and store it in this allele's hash table.
				if ($currentAllele ne ""){
					print_debug("Saving allele $currentAllele data.\n");
          $alleleLocalData{$currentAllele}=join('',@currentAlleleData);
					#print_debug("$alleleLocalData{$currentAllele}\n");
					@currentAlleleData=(); # Clear the working array to store the next allele data.
				}

				$currentAllele=$1;
				print_debug("Found allele $currentAllele\n");
      }
			else{
        # Store the line of base sequences into the array to be joined later.
				push(@currentAlleleData,$_);
			}
    }
		# Join the final allele's data into one long sequence string and store it.
		if ($currentAllele ne ""){
			print_debug("Saving allele $currentAllele data.\n");
			$alleleLocalData{$currentAllele}=join('',@currentAlleleData);
			#print_debug("$alleleLocalData{$currentAllele}\n");
			@currentAlleleData=(); # Clear the working array to store the next allele data.
		}

		close($allelefile);

    # At this point, all the allele data in a given allele file has been processed, so store in the global array.
    # Note that this is stored as a reference to a COPY of the local hash table (because the local hash table var
		# will be reused when processing the next allele file).
    push(@alleleGlobalData,{%alleleLocalData});

		print_debug("Keys of alleleLocalData are: \n");
		if ($debug eq "Y"){
      my $debugKey;
			foreach $debugKey (keys(%alleleLocalData)){
        print_debug("$debugKey\n");
			}
		}

	}
}

sub process_parameters{
  # Process main program parameters and verify command line arguments.

	$arg_num=0;
	while ($arg_num <= $#ARGV){
		given ($ARGV[$arg_num]){
			when ("-debug"){
				$debug="Y";
				print "Producing debug output.\n";
			}
		  when ("-quiet") {
				print "Producing quiet (no stdout) output.  Log file is still generated.\n";
				$quiet="Y";
			}
		  when ("-c"){
				$arg_num++;
				$configFilename="$ARGV[$arg_num]";
				if ( ! -r $configFilename ){
					print "Unable to read config file $configFilename!";
					exit 1;
				}
			}
		  when ("-l"){
				$arg_num++;
				$logFilename="$ARGV[$arg_num]";
				print "Writing program log to $logFilename\n";
			}
			when ("-fast"){
        $fast="Y";
				print "Will use fast evaluation (won't check for duplicates during processing).\n";
			}
		  default{
				print "*** ERROR: Unknown command line argument or file $ARGV[$arg_num].\n";
				exit 1;
		  }
    }
		$arg_num++;
	}

	open($logfile,">","$logFilename") or die "Unable to open log file $logFilename for writing! $!\n";
	open($configfile,"<","$configFilename") or die "Unable to open config file $configFilename for reading! $!\n";

  # Process the config file, line by line.

	my $currentConfigSection="";
	my $configLine=0;
	my $configLineText="";
	while (<$configfile>){

		$configLine++;

		chomp;
		$configLineText=$_;

		# Ignore blank lines and comments.
		if ((/^#/) || (/^ *$/)){
			print_debug("Ignoring config comment: $_\n");
			next;
		}
		# Look for section tags.
		elsif (/^\[ST ALLELE TABLE FILES\]$/){
			$currentConfigSection="TABLE";
			print_debug("Now in section $currentConfigSection\n");
			next;
		}
		elsif (/^\[ALLELE FILES\]$/){
			$currentConfigSection="ALLELE";
			print_debug("Now in section $currentConfigSection\n");
			next;
		}
		elsif (/^\[REF FILES\]$/){
			$currentConfigSection="REF";
			print_debug("Now in section $currentConfigSection\n");
			next;
		}
		# Else the line should be a file name (with or without path).  Check if it can be read and push it into the appropriate section variable array.
		elsif ( -r ){
			print_debug("Now checking file $_\n");

			if ($currentConfigSection eq ""){
				print_all("*** ERROR: Unknown line in config (line $configLine) -- file entry found outside of a known section?\n");
				print_all("Line is: $_\n");
				exit 2;
			}

			given ($currentConfigSection){
				when ("TABLE"){ push(@sttableFiles,$configLineText); }
				when ("ALLELE"){ push(@alleleFiles,$configLineText); }
				when ("REF"){ push(@referenceFiles,$configLineText); }
				default { print_all("*** ERROR: Unknown section $currentConfigSection!\n"); }
			}

		}
		else{
			print_all("*** ERROR: Unknown line in config (line $configLine) -- not a known section or readable file or ignorable comment!\n");
			print_all("Line is: $_\n");
			exit 2;
		}
	}
}

sub print_debug{
  # Prints output to STDOUT and log file if debug flag is set.  Otherwise nothing.
  # If the quiet function is also set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($debug eq "Y"){
    if ($quiet eq "N"){
      print "$_[0]";
    }
    print $logfile "$_[0]";
  }
}

sub print_all{
  # Silly MUX'ed function to write output to both standard out and logfile in one call.
  # If the quiet function is set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($quiet eq "N"){
    print "$_[0]";
  }
  print $logfile "$_[0]";
}


