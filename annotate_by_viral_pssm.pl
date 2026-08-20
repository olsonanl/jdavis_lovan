#! /usr/bin/env perl

use strict;
use JSON::XS;
use File::Slurp;
use Getopt::Long;
use Cwd;
use gjoseqlib;

my $default_data_dir = $ENV{LOWVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

#
# Exit codes.
#
# Before these existed the script returned 1 for both "this genome cannot be placed" and
# "you passed a bad flag", and 0 for "classified, but nothing scored" -- so a caller could
# not tell a legitimate empty result from a misuse or a crash, and the GTO wrapper had to
# grep stderr for two English sentences to find out.  These give the caller the answer
# directly.
#
# The two "nothing to annotate" codes are the point: both are legitimate outcomes for a
# 200 bp fragment, not errors, and a caller that runs the downstream stages should skip
# them rather than fail the job.
#
# 10 and up rather than 3 and up because Perl's `die` exits with $! when $! happens to be
# set, which reaches the low single digits routinely (ENOENT is 2).  Anything outside this
# table -- an uncaught die, 255, a signal at 128+n -- is an unexpected failure.
#
use constant EXIT_OK            => 0;    # features called
use constant EXIT_USAGE         => 1;    # bad flags, unknown taxon, unreadable config
use constant EXIT_INTERNAL      => 2;    # the script itself broke
use constant EXIT_NO_REFERENCE  => 10;   # no reference contig scored above -mcb
use constant EXIT_NO_FEATURES   => 11;   # taxon chosen, but no PSSM cleared its bit_cutoff
use constant EXIT_OUT_OF_BOUNDS => 12;   # input length outside -min/-max

#
# Fail with a fixed code instead of `die`, whose status is whatever $! held.
#
sub bail
{
	my($code, @msg) = @_;
	print STDERR @msg;
	exit($code);
}

my $usage = 'annotate_by_viral_pssm.pl [options] -i subject_contig(s).fasta 

		-h   help
		-i   Input subject contigs in fasta format
		-t   Declare a temp file (d = random)
		-tax Declare a taxonomy id (D = 10239)
		-vtax Declare the annotation taxon, skipping BLASTn-based detection.
		     Must be one of the names printed by -list-vtax (a PSSM directory name).
		     The classification BLASTn still runs and warns if it disagrees, and
		     the -mcb rejection is downgraded to a warning.
		     NOTE: this value also becomes viral_family in the GTO, which selects
		     the Transcript-Editing/ and Splice-Variants/ dirs used downstream.
		-list-vtax  Print the valid -vtax values to STDOUT and exit.
		-vfam Declare the annotation *family* when the genus is not known, and let
		     BLASTn pick the genus within it.  The reference sweep is restricted to
		     the reps of the taxa in that family, so a record whose lineage stops at
		     family ("Bat coronavirus", "Paramyxoviridae sp.") can still be placed.
		     Use -list-vfam for the valid names.  Unlike -vtax this still has to make
		     a call, so -mcb applies normally and a genome that clears nothing is
		     still rejected -- pair it with a lower -mcb to work the 50-150 bit band,
		     where a family-restricted call is right about 93% of the time against
		     roughly 60% for an unrestricted one.  Cannot be combined with -vtax or
		     -skip-classification.  -margin works and is more meaningful here: the
		     runner-up is a sibling genus rather than an unrelated family.
		-list-vfam  Print the valid -vfam values to STDOUT and exit.
		-skip-classification  Only BLASTn the reference contigs belonging to the
		     declared -vtax taxon (1-14 files) instead of every rep in -c.  Requires -vtax.  Cuts the
		     dominant cost of a run; the tradeoff is that no cross-check against
		     the other taxa is possible, so a mistaken -vtax will not be caught.
		     Intended for bulk reruns where the taxon is already known to be right.
		-g   Genome name (D = Viruses);
		-p   Prefix for the output files (.ffn, .faa, and .tbl files)
		-s   Append sequences to the feature table (default is that they are left off)
		-ks  If a pssm match extends beyond a stop codon it will keep the entire region of the match including the stop. 
		-threads number of blast threads in the blastn and tblastn

		-min Minimum contig	length (d = 300)  # otherwise the genome is rejected
		-max Maximum contig length (d = 35000) # for reference Measles is 15,894 and Beilong is 19,212
	
		-mcb Minimum contig bitscore to enable annotation (d = 150) #otherwise the genome is rejected.

		-margin Minimum ratio between the bit score of the winning taxon and that of
		     the runner-up taxon (e.g. -margin 1.3).  OFF by default: without it
		     nothing is computed, warned about, or written.  When given, the
		     classification is scored for confidence as well as absolute strength --
		     a genome that clears -mcb but beats the runner-up by only a few bits is
		     a near-tie, and a near-tie in the wrong direction is a mis-annotation
		     rather than a failed job.  The genome is still annotated: this warns on
		     STDERR and writes the numbers to <prefix>.classification, it does not
		     reject.  Not available with -skip-classification, which searches the
		     reps of one taxon only and so has no runner-up to compare against.

        -j   Full path to the options file in JSON format which carries data for a match (D = $default_data_dir/Viral_PSSM.json)
		-c   Representative contigs directory (D = $default_data_dir/Viral-Rep-Contigs)
		-pssm   Base directory of PSSMs   (D = $default_data_dir/Viral-PSSMs)
	           Note that this is set up as a directory of pssms
	           right now this is hardcoded as: "virus".pssms within this directory.

	      
	   Debug parms (turns off output file generation)
	   	-tmp keep temp dir
		-no  no output files generated, to be used in conjunction with one of the following:
		    -dna print only genes to STDOUT 
		    -aa print proteins to STDOUT
		    -tbl print only feature table to STDOUT
		    -ctbl [file name] concatenate table results to a file (for use with many genomes)


 	## Output Table columns are the following:
 			0	genome_id	
 			1	genome_name	
 			2	accession (original contig id)
 			3	annotation_source (always LV)	
 			4	feature_type	
 			5	a made up feature id	
 			6	gene symbol	
 			7	start	
 			8	end	
 			9	strand	
 			10	na_length	
 			11	annotation taxon (BLASTn-detected, declared with -vtax, or picked
 			    within a declared family with -vfam)
 			12	closest genome (best rep within the annotation taxon)
 			13	closest genome bit score (BLASTn of contigs)
 			14  closest genome id(s)
 			15  closest genome name
 			16	PSSM
 			17	NT sequence (if -s append sequences)	
 			18	AA sequence (if -s append sequences)


 	## Exit status.  10 and up are the "nothing was annotated" outcomes: all three are
 	## legitimate results for a short or divergent fragment rather than errors, and a
 	## caller that chains the downstream stages should skip them rather than fail the job.

 			0	features were called
 			1	usage or configuration error (bad flag, unknown taxon, unreadable JSON)
 			2	internal error
 			10	no reference contig scored above -mcb, so no taxon could be chosen
 			11	a taxon was chosen, but no PSSM cleared its bit_cutoff
 			12	input length outside -min/-max
';

my ($help, $opt_file, $contig_file, $tmp, $tax, $vtax, $list_vtax, $vfam, $list_vfam, $skip_class, $keep_stop, $genome_name, $cdir, $pdir, $keep_temp, $min_len, $max_len, $aa_only, $dna_only, $tbl_only, $no_out, $ctbl, $prefix, $append_seqs, $threads, $min_contig_bit, $min_margin);

my $opts = GetOptions( 'h'         => \$help,
                       'tmp'       => \$keep_temp,
                       'no'        => \$no_out,
                       'aa'        => \$aa_only,
                       'dna'       => \$dna_only,
                       'tbl'       => \$tbl_only,
                       'ctbl=s'    => \$ctbl,
                       'i=s'       => \$contig_file,
                       'tax=s'     => \$tax,
                       'vtax=s'    => \$vtax,
                       'list-vtax' => \$list_vtax,
                       'vfam=s'    => \$vfam,
                       'list-vfam' => \$list_vfam,
                       'skip-classification' => \$skip_class,
                       'threads=s' => \$threads,
                       't=s'       => \$tmp,
                       'g=s'       => \$genome_name,
                       'c=s'       => \$cdir,
                       'pssm=s'    => \$pdir,
                       'min=s'     => \$min_len,
                       'max=s'     => \$max_len,
                       'mcb=f'     => \$min_contig_bit,
                       'margin=f'  => \$min_margin,
                       'j=s'       => \$opt_file,
                       'ks'        => \$keep_stop,
                       'p=s'       => \$prefix,
                       's'         => \$append_seqs); 

if ($help){print "$usage\n"; exit(EXIT_OK);}   # asking for help is not an error
unless ($contig_file || $list_vtax || $list_vfam){bail(EXIT_USAGE, "must declare an input subject file with -i \n\n$usage\n");}

# Name the Temp file:
# generates a random 20 digit string of 0-9a-f
unless ($tmp){$tmp .= sprintf("%x", rand 16) for 1..20;}
unless ($min_len){$min_len = 300; }
unless ($max_len){$max_len = 35000; }
unless ($min_contig_bit){$min_contig_bit = 150;}
unless ($tax){$tax = "10239"; }
unless ($genome_name){$genome_name = "Viruses"; }
unless ($prefix){$prefix = "Viral_Annotation";}
unless ($cdir){$cdir = "$default_data_dir/Viral-Rep-Contigs"; }
unless ($pdir){$pdir = "$default_data_dir/Viral-PSSMs"; }
unless ($opt_file){$opt_file = "$default_data_dir/Viral_PSSM.json"; }
unless ($threads){$threads = 8}


## This is the json file with all of the protein-specific information for 
## Customizing the annotation.
open (IN, "<$opt_file") or bail(EXIT_USAGE, "Cant find JSON options file, use -opt\n");
my $options = eval { decode_json(scalar read_file(\*IN)) };
bail(EXIT_INTERNAL, "Could not parse JSON options file $opt_file: $@") if $@;
close IN;

## Resolve the annotation taxon.  By default it is detected by BLASTn against
## Viral-Rep-Contigs (below); -vtax lets the user declare it instead.
## Do this before makeblastdb / the temp dir / the output filehandles so that a
## bad name fails without leaving anything behind on disk.

my @valid_vtax = valid_vtax_list($pdir, $options);
my %taxon_family = taxon_family_map($options, \@valid_vtax);
my @valid_vfam   = do { my %s; $s{$_}++ for values %taxon_family; sort keys %s };

if ($list_vtax)
{
	print "$_\n" foreach @valid_vtax;
	exit(EXIT_OK);
}

if ($list_vfam)
{
	print "$_\n" foreach @valid_vfam;
	exit(EXIT_OK);
}

# -vfam and -vtax are different claims about how much the caller knows, and they cannot
# both be true: -vtax says "annotate as this, do not decide", -vfam says "decide, but only
# among these".  Combining them would silently ignore one.
# Checked before the -skip-classification rule below so the more specific message wins.
if ((defined $vfam) && (defined $vtax))
{
	bail(EXIT_USAGE,
	    "-vfam and -vtax cannot be combined: -vtax declares the annotation taxon outright,\n"
	  . "-vfam asks BLASTn to choose one within a family. Use whichever matches what is known.\n");
}
if ((defined $vfam) && $skip_class)
{
	bail(EXIT_USAGE,
	    "-vfam and -skip-classification cannot be combined: -vfam already restricts the sweep,\n"
	  . "and the genus still has to be chosen from the reps it leaves.\n");
}

if ($skip_class && (! defined $vtax))
{
	bail(EXIT_USAGE, "-skip-classification only makes sense with -vtax: there is no declared taxon to restrict the sweep to.\n");
}

# The margin is winner/runner-up across taxa, so it needs the full sweep.  Under
# -skip-classification only one taxon's reps are searched and there is no runner-up
# to divide by; fail loudly rather than silently reporting an infinite margin.
if (defined $min_margin)
{
	bail(EXIT_USAGE, "-margin needs a ratio greater than 1 (got $min_margin); 1 would accept any winner.\n")
		if $min_margin <= 1;
	bail(EXIT_USAGE,
	    "-margin cannot be combined with -skip-classification: only the declared taxon's\n"
	  . "reference contigs are searched, so there is no runner-up taxon to compare against.\n")
		if $skip_class;
}

if (defined $vtax)
{
	$vtax =~ s/^\s+|\s+$//g;
	$vtax =~ s/\.pssms$//;   # tolerate a pasted PSSM directory name

	my %exact = map {($_ => 1)} @valid_vtax;
	my %fold  = map {(lc($_) => $_)} @valid_vtax;

	# Exact match first.  This is what keeps -vtax Orthopneumovirus from ever
	# resolving to Orthopneumovirus_muris -- never use prefix/substring matching here.
	# It also settles the eight Bunyavirales names, which are both annotation taxa and
	# families: as taxa they win here, and -vfam on them would mean the same thing
	# anyway because each is the only member of its own family.
	unless ($exact{$vtax})
	{
		my %fam_fold = map {(lc($_) => $_)} @valid_vfam;

		if (exists $fold{lc $vtax})
		{
			print STDERR "Interpreting -vtax '$vtax' as '$fold{lc $vtax}'\n";
			$vtax = $fold{lc $vtax};
		}
		elsif (exists $fam_fold{lc $vtax})
		{
			# A family name given to -vtax is the natural thing to type when the
			# lineage stopped at family.  Honour it as -vfam rather than dying on a
			# technicality, but say so: the semantics differ (this one still decides).
			$vfam = $fam_fold{lc $vtax};
			print STDERR "-vtax '$vtax' is a family, not an annotation taxon; "
			           . "treating it as -vfam $vfam (BLASTn will choose the genus within it).\n";
			undef $vtax;
		}
		else
		{
			bail(EXIT_USAGE,
			    "Unknown -vtax value '$vtax'.\n"
			  . "Valid annotation taxa (PSSM dirs in $pdir that also have a $opt_file entry):\n"
			  . join("", map {"\t$_\n"} @valid_vtax)
			  . "Use -list-vtax to print this list, or -vfam / -list-vfam to declare only a family.\n");
		}
	}
}

# -vfam: restrict rather than declare.  Resolved the same way as -vtax -- exact first,
# case-fold as a courtesy -- against the family names carried in the JSON options file.
if (defined $vfam)
{
	$vfam =~ s/^\s+|\s+$//g;

	my %exact_f = map {($_ => 1)} @valid_vfam;
	my %fold_f  = map {(lc($_) => $_)} @valid_vfam;

	unless ($exact_f{$vfam})
	{
		if (exists $fold_f{lc $vfam})
		{
			print STDERR "Interpreting -vfam '$vfam' as '$fold_f{lc $vfam}'\n";
			$vfam = $fold_f{lc $vfam};
		}
		else
		{
			bail(EXIT_USAGE,
			    "Unknown -vfam value '$vfam'.\n"
			  . "Valid families (the \"family\" field of each taxon in $opt_file):\n"
			  . join("", map {"\t$_\n"} @valid_vfam)
			  . "Use -list-vfam to print this list.\n");
		}
	}
}


# Version the Genome (another random alpha numeric)
my $version;
{$version .= sprintf("%x", rand 16) for 1..6;}  ## 11M combos

# Open file handles for the outputs

unless($no_out)
{
	if ($prefix)
	{
		open (AA,  ">$prefix.faa");
		open (DNA, ">$prefix.ffn");
		open (TBL, ">$prefix.feature.tbl");
	}
	else
	{
		open (AA, ">$tax.$version.faa");
		open (DNA, ">$tax.$version.ffn");
		open (TBL, ">$tax.$version.feature.tbl");
	}
}

if ($ctbl)
{
	open(CTBL, ">>", "$ctbl")  
}

# Make a hash out of the subject contigs
open (IN, "<$contig_file");
my @seqs = &gjoseqlib::read_fasta(\*IN);
close IN;

# Make sure it falls within the min and max contig length.
# Count up the size so that it can handle multiple contigs.
my $len = 0;
my %contigH;
my @contig_order;
for my $i (0..$#seqs)
{
	$contigH{$seqs[$i][0]} = uc ($seqs[$i][2]);
	$len += length($seqs[$i][2]);
	push @contig_order, $seqs[$i][0];
}
if (($len < $min_len) || ($len > $max_len))
{
	bail(EXIT_OUT_OF_BOUNDS, "Unexpected input genome size equal to: $len\nMin = $min_len\tMax = $max_len\n");	
}	

# Make the temp dir.
my $base = getcwd;
mkdir ($tmp); 
system "cp $contig_file $tmp";
my $s_file = $contig_file;
$s_file =~ s/.+\///g; 

chdir ($tmp);
system "makeblastdb -dbtype nucl -in $s_file >/dev/null";


#   Select the correct species based on the contig blastN
#   If this begins to break down, new reference contigs can be added to the directory
#   Formatting them as genus.version.dna

opendir (DIR, "$cdir");
my @reps = sort grep{$_ !~ /^\./}readdir(DIR);   # sort: reproducible reference pick on ties (gist #19)
closedir(DIR);

# -skip-classification: drop every rep that does not belong to the declared taxon.
# The remaining blastn is no longer a classification -- it exists only to pick the
# closest reference within $vtax for the close-genome columns.  Same taxon-name
# derivation as the loop below, so Orthopneumovirus cannot pull in Orthopneumovirus_muris.
if ($skip_class)
{
	my $before = scalar @reps;
	@reps = grep { my $t = $_; $t =~ s/\..+//g; $t eq $vtax } @reps;
	unless (@reps)
	{
		bail(EXIT_USAGE, "-skip-classification: no reference contigs for -vtax $vtax in $cdir\n");
	}
	printf STDERR "Restricting BLASTn sweep to %d of %d reference contigs (-skip-classification, -vtax %s)\n",
	       scalar @reps, $before, $vtax;
}

# -vfam: same filter, widened from one taxon to one family.  What is left is still a
# classification -- the winner among these reps becomes the annotation taxon -- so unlike
# -skip-classification this does not weaken the -mcb decision, it only narrows the field
# it is made over.  Restricting can never raise a bit score, so this cannot rescue a
# genome that scores below -mcb against everything; what it buys is that the call made in
# the band above -mcb is constrained to the family the submitter already told us.
if (defined $vfam)
{
	my $before = scalar @reps;
	@reps = grep { my $t = $_; $t =~ s/\..+//g;
	               defined $taxon_family{$t} && $taxon_family{$t} eq $vfam } @reps;
	unless (@reps)
	{
		bail(EXIT_USAGE, "-vfam: no reference contigs for any taxon of family $vfam in $cdir\n");
	}
	my %in_scope = map {my $t = $_; $t =~ s/\..+//g; ($t => 1)} @reps;
	printf STDERR "Restricting BLASTn sweep to %d of %d reference contigs (-vfam %s: %s)\n",
	       scalar @reps, $before, $vfam, join(", ", sort keys %in_scope);
}

my $best_contig_bit = 0;
my $best_virus_match;
my $best_rep;
my %taxon_best;   # taxon => [best bit, best rep] -- lets -vtax report its own taxon's score
foreach (@reps)
{
	my $rep = $_;
	my $virus = $_;
	$virus =~ s/\..+//g; 
	
	my $rep_file = "$cdir/$rep";
	open (IN, "blastn -query $rep_file -subject $s_file -evalue 0.5 -reward 2 -penalty -3 -word_size 11 -outfmt 15 -soft_masking false -num_threads $threads 2>/dev/null |") or die "Could not run blastn: $!"; 
	my $blastn = decode_json(scalar read_file(\*IN));	
	close IN;

	my $match_bit = get_blastn_bit($blastn);  
	unless ($match_bit < $best_contig_bit)
	{
		$best_contig_bit = $match_bit;
		$best_virus_match = $virus;
		$best_rep = $rep;
	}

	# Track the best rep per taxon as well, so a declared -vtax can report the
	# closest reference *within* the taxon it is actually annotating with.
	# !exists guarantees an entry even when every rep for a taxon scores 0.
	if ((! exists $taxon_best{$virus}) || ($match_bit >= $taxon_best{$virus}->[0]))
	{
		$taxon_best{$virus} = [$match_bit, $rep];
	}
}	

# The taxon we annotate with: a user-declared -vtax wins over the BLASTn call.
# This is the only place the annotation taxon is chosen; everything below keys off $virus.

my $virus = defined $vtax ? $vtax : $best_virus_match;

# Columns 12-13 (closest genome / bit score) track the taxon we are actually annotating
# with, not the overall BLASTn winner.  Reporting an out-of-taxon reference here would
# put a permanently mislabeled record in the GTO close_genomes list, which carries no
# taxon field of its own.  The cross-taxon result is preserved in the warning below.

my ($anno_bit, $anno_rep) = defined $vtax
                          ? @{ $taxon_best{$vtax} || [0, undef] }
                          : ($best_contig_bit, $best_rep);

# here we make it fail gracefully if we don't enounter a reference genome with a high 
# enough blastn bit score.  If the user declared the taxon they have taken that call,
# so warn and continue instead of rejecting the genome.

if ($best_contig_bit < $min_contig_bit)
{
	if (defined $vtax)
	{
		my $scope = $skip_class ? "best within $vtax" : "best overall";
		print STDERR "WARNING: no reference contig scored above -mcb $min_contig_bit ($scope = "
		           . (defined $best_rep ? $best_rep : "none") . ", bit = $best_contig_bit); "
		           . "continuing because -vtax $vtax was declared.\n";
	}
	elsif (defined $vfam)
	{
		# The family was declared but nothing in it scored.  Say which reps were
		# actually searched, so this is not confused with the unrestricted failure:
		# a -vfam run that fails here would have failed unrestricted too.
		print STDERR "No reference contig of family $vfam scored above -mcb $min_contig_bit "
		           . "(best = " . (defined $best_rep ? $best_rep : "none") . ", bit = $best_contig_bit).\n";
		exit(EXIT_NO_REFERENCE);
	}
	else
	{
		print STDERR "No matching reference contigs with bit score greater than $min_contig_bit\n";
		exit(EXIT_NO_REFERENCE);
	}
}

# Not reachable under -skip-classification: only $vtax's own reps were searched, so
# $best_virus_match is either $vtax or undef.  Guarded explicitly rather than relying on that.
if ((! $skip_class) && (defined $vtax) && (defined $best_virus_match) && ($best_virus_match ne $vtax))
{
	print STDERR "WARNING: -vtax declared $vtax, but BLASTn classification chose $best_virus_match "
	           . "(rep $best_rep, bit = $best_contig_bit). Best hit within declared $vtax = "
	           . (defined $anno_rep ? "$anno_rep, bit = $anno_bit" : "none")
	           . ". Annotating as $vtax.\n";
}

# Margin scoring (-margin).  Off unless a ratio was given.
#
# -mcb asks "is the best hit strong enough"; the margin asks "is it better than the
# alternative".  Those come apart: a genome can clear 150 comfortably and still be a
# near-tie between two taxa, and that is the shape of a wrong call.  The reference
# additions of 2026-08-19 produced exactly one wrong accepted call in a 500-genome
# validation set, and its margin (1.26) was the smallest of all 363 accepted calls
# against a median of 15.4 -- unremarkable in bit score, an outlier in margin.
#
# %taxon_best already holds the per-taxon best from the sweep, so the runner-up costs
# no extra BLAST.  This warns and records; it never rejects, because for the failure
# set a wrong-but-flagged annotation is more useful than another dead job.

#
# The sidecar is also written under -vfam, with no -margin.  A family-restricted call is
# a decision the caller asked us to make on partial information, usually in the band where
# -mcb was relaxed, so the losing siblings and their scores are the evidence for it and
# have to survive the run.  It is also what lets the GTO wrapper recover the chosen genus
# when the annotation comes back empty -- there is no feature table to read it off.
#
my ($margin, $runner_up, $runner_up_bit);
if ((defined $min_margin) || (defined $vfam))
{
	# rank taxa by their best bit score; the winner is $best_virus_match by construction
	my @ranked = sort {$taxon_best{$b}->[0] <=> $taxon_best{$a}->[0]} keys %taxon_best;
	if (@ranked >= 2)
	{
		$runner_up     = $ranked[1];
		$runner_up_bit = $taxon_best{$runner_up}->[0];
	}

	# No runner-up, or a runner-up that scored zero, means nothing else was in
	# contention.  That is maximum confidence, not a division by zero.
	$margin = (defined $runner_up_bit && $runner_up_bit > 0)
	        ? $best_contig_bit / $runner_up_bit
	        : undef;

	my $shown = defined $margin ? sprintf("%.2f", $margin) : "inf";
	if (! defined $min_margin)
	{
		# -vfam without -margin: report, do not judge.  Within a family the runner-up
		# is a sibling genus, so this ratio is a more meaningful number than the
		# cross-family one -- but no threshold has been calibrated for it either.
		print STDERR "Within-family margin $shown"
		           . (defined $runner_up ? " (runner-up $runner_up, bit = $runner_up_bit)" : " (no runner-up taxon in $vfam)")
		           . "\n";
	}
	elsif (defined $margin && $margin < $min_margin)
	{
		printf STDERR "WARNING: low classification margin %.2f (< -margin %s): %s scored %s but "
		            . "%s scored %s. The winning taxon is only %.0f%% ahead of the runner-up; "
		            . "treat this call as unconfirmed.\n",
		            $margin, $min_margin, $best_virus_match, $best_contig_bit,
		            $runner_up, $runner_up_bit, 100 * ($margin - 1);
	}
	else
	{
		print STDERR "Classification margin $shown"
		           . (defined $runner_up ? " (runner-up $runner_up, bit = $runner_up_bit)" : " (no runner-up taxon)")
		           . "\n";
	}

	# Sidecar, so the numbers survive the run.  Written whenever -margin is given,
	# not only on a warning -- a margin is evidence either way, and the GTO wrapper
	# reads this file to attach it to close_genomes.
	#
	# $base, not the cwd: we are inside the temp dir at this point (chdir above), and
	# the wrapper deletes it.  This is where the .tbl/.ffn/.faa outputs land too --
	# they are written after the chdir back, so an absolute -p must not be re-rooted.
	my $cls_file = ($prefix =~ m,^/,) ? "$prefix.classification" : "$base/$prefix.classification";
	if (open(my $cls, ">", $cls_file))
	{
		print $cls "taxon\t$best_virus_match\n";
		print $cls "reference\t" . (defined $best_rep ? $best_rep : "-") . "\n";
		print $cls "bit\t$best_contig_bit\n";
		print $cls "runner_up\t" . (defined $runner_up ? $runner_up : "-") . "\n";
		print $cls "runner_up_bit\t" . (defined $runner_up_bit ? $runner_up_bit : 0) . "\n";
		print $cls "margin\t$shown\n";
		if (defined $min_margin)
		{
			print $cls "margin_threshold\t$min_margin\n";
			print $cls "below_threshold\t" . ((defined $margin && $margin < $min_margin) ? 1 : 0) . "\n";
		}
		# Scope of the sweep the numbers above came out of.  Without this a -vfam
		# ranking reads like a full one that happened to find only sibling genera.
		if (defined $vfam)
		{
			print $cls "scope\tfamily\n";
			print $cls "family\t$vfam\n";
			print $cls "min_contig_bit\t$min_contig_bit\n";
		}
		print $cls "annotated_as\t$virus\n";
		# full ranking, so a near-tie can be inspected without rerunning the sweep
		for my $t (@ranked)
		{
			last if $taxon_best{$t}->[0] <= 0;
			print $cls "score\t$t\t$taxon_best{$t}->[0]\t$taxon_best{$t}->[1]\n";
		}
		close $cls;
	}
	else
	{
		print STDERR "WARNING: could not write $cls_file: $!\n";
	}
}

#get closest genome data, for the taxon we are annotating with
my ($best_rep_ids, $best_rep_name) = ("", "");
if (defined $anno_rep)
{
	$best_rep_ids  = $options->{$virus}->{close_genomes}->{$anno_rep}->{genome_ids};
	$best_rep_name = $options->{$virus}->{close_genomes}->{$anno_rep}->{genome_name};
}

my $how = defined $vtax ? ($skip_class ? "\t(taxon declared with -vtax, sweep restricted)"
                                       : "\t(taxon declared with -vtax)")
        : defined $vfam ? "\t(genus chosen within declared family $vfam)"
        :                 "";
print STDERR "-----------------------\nAnnotating as $virus\tBit = $anno_bit$how"
           . "\n-----------------------\n";
opendir (DIR, "$pdir/$virus.pssms");
my @pssm_dirs = sort grep{$_ !~ /^\./}readdir(DIR);  # reads the directory of PSSM dirs (sorted: gist #19)
closedir(DIR);

my @all_seqs;
my %positions;
#my ($upstream_ext, $downstream_ext, $keep_stop, $bit_cutoff, $cov_cutoff);

my $non_pssm_feat = {};
foreach (@pssm_dirs)  #Each PSSM dir contains one or more PSSMs for a given homolog
{
	my $pssmdir = $_; 
	opendir (DIR, "$pdir/$virus.pssms/$pssmdir");
	my @pssm_files = sort grep{$_ !~ /^\./}readdir(DIR); # gets the pssms (sorted: gist #19)
	closedir(DIR);

	# retrieve the protein specific options for the corresponding PSSM.
	unless ($options->{$virus}->{features}->{$pssmdir})
	{
		print STDERR "No data in JSON file for VIRUS: $virus\tPSSM: $pssmdir\n"; 
		next;
	}
	my $upstream_ext   = $options->{$virus}->{features}->{$pssmdir}->{upstream_ext};
	my $downstream_ext = $options->{$virus}->{features}->{$pssmdir}->{downstream_ext};
	my $bit_cutoff     = $options->{$virus}->{features}->{$pssmdir}->{bit_cutoff};
	my $cov_cutoff     = $options->{$virus}->{features}->{$pssmdir}->{coverage_cutoff};
	my $start_to_met   = $options->{$virus}->{features}->{$pssmdir}->{start_to_met};
	my $feature_type   = $options->{$virus}->{features}->{$pssmdir}->{feature_type};
	my $anno           = $options->{$virus}->{features}->{$pssmdir}->{anno};
	my $symbol         = $options->{$virus}->{features}->{$pssmdir}->{gene_symbol};
		
	print STDERR "\t$virus\t$pssmdir\t$anno\tbit\t$bit_cutoff\tcov\t$cov_cutoff\tkeep_stop\t$keep_stop\tupstream_ext\t$upstream_ext\tdownstream_ext\t$downstream_ext\n"; 		

	#   Select the best pssm per protein
	#   If this begins to break down, new reference pssms can be added to the 
	#   "Taxon.pssms/protein/"   directory
	#   Adding more pssms should theoretically improve the accuracy

	my ($best_pssm, $best_results, $pssm_name);
	my $best_bit = 0;

	foreach (@pssm_files)  ##  <---- This is where I would need to lookup the rules for each protein
	{		
		my $pssm_file = "$pdir/$virus.pssms/$pssmdir/$_";
		my $name = $_;
		open (IN, "tblastn -outfmt 15 -db $s_file -in_pssm $pssm_file -seg no -num_threads $threads |") or die "Could not run tblastn: $!";
		my $pssm_blast = decode_json(scalar read_file(\*IN));		
		close IN; 	
		
		my  ($results, $hsp_best_bit) = matching_tblastn_hsps_json($pssm_blast, $bit_cutoff, $cov_cutoff);
		print STDERR "\t$name\t$hsp_best_bit\n"; 

		unless ($hsp_best_bit < $best_bit)   #evaulate each pssm blast based on the bit score, and pick the best one
		{
			$best_bit = $hsp_best_bit;
			$best_pssm = $name;
			$best_results = $results;
		}
	}
	print STDERR "\n\tChoosing $best_pssm\t$best_bit\n\n"; 
	my $nhsps = scalar @$best_results;
	for my $i (0..($nhsps -1))
	{		
		my $contig = $best_results->[$i]->{contig};
		$contig =~ s/\s.+//g;								
		my $hseq = $best_results->[$i]->{hseq};
		$hseq =~ s/\-//g; # eliminate alignment gaps in hit seq
			
		print STDERR "\t$contig\t$anno\t$best_results->[$i]->{hit_from}\t$best_results->[$i]->{hit_to}\t$best_results->[$i]->{frame}\t$best_results->[$i]->{bit}\n";

		## Okay, Now we need to find our features.
	
	
		# lookup DNA Coordinates
		my ($from, $to);
		my $strand = "+";
		if ($best_results->[$i]->{frame} < 0)
		{
			$from = $best_results->[$i]->{hit_to};
			$to   = $best_results->[$i]->{hit_from};
			$strand = "-";
		}
		else
		{
			$from = $best_results->[$i]->{hit_from};
			$to   = $best_results->[$i]->{hit_to};
		}

		###  $from and $to are the coordinates of the protein from the blast without the stop codon.
		###  The end extension loop below, scans in the 3' direction for the next stop codon.
		###  if $downstream_ext is declared,  it will get new gene coordinates with the stop
		###  codon being included in the gene coordinates by convention. 		

		my ($gene_begin, $gene_end);
		if ((! $downstream_ext) || ($hseq =~ /\*/))
		{
			$gene_begin = $from;
			$gene_end = $to;
		}
		else
		{
			($gene_begin, $gene_end)= scan_to_stop_codon ($from, $to, $contigH{$contig});
			print STDERR "\tScanning to Stop codon, Original Coords: $from to $to\tNew Coords: $gene_begin to $gene_end\n"; 	

		}
		#print STDERR "COV = $best_results->[$i]->{cov}\n";

		
		## If a match contains an internal stop, everything after the stop is cropped by default.		
		my $new_end;
		if (($hseq =~ /\*/) && (! $keep_stop))
		{
			print STDERR "\tMatch contains a stop codon, cropping\n";
			$new_end = crop_to_stop_codon($gene_begin, $gene_end, $hseq);
		
			#check the coverage.
			# Coverage of the RETAINED (pre-stop) portion of the hit vs the query.
			# Prior versions divided \$len (whole-genome length) here, so this check
			# never fired.  Corrected below -- NOTE this can now drop low-coverage
			# features that previously slipped through (gist #4).
			(my $hseq_pre_stop = $hseq) =~ s/\*.+//s;
			my $cov2 = length($hseq_pre_stop) / (abs(($best_results->[$i]->{q_to} - $best_results->[$i]->{q_from})+1));

			if ($cov2 < $cov_cutoff)
			{
				print STDERR "\tCoverage prior to stop codon \($cov2\) is lower than cutoff: $cov_cutoff\n";
				next;
			}
			print STDERR "\tCropping to coordinates: $gene_begin\tto\t$new_end\n"; 			
		}	
		if ($new_end){$gene_end = $new_end};	
		

		# When the pssm match doesn't start with an M, it will search to the left until it finds the first one.
		# This must be turned off if there is a non-AUG start.

		my ($new_gene_begin, $new_gene_end); 
		if (($hseq !~ /^M/i) && ($upstream_ext))
		{
			($new_gene_begin, $new_gene_end) = scan_to_met_start( $gene_begin, $gene_end, $contigH{$contig});
			print STDERR "\tScanning for Met start, Original Coords: $gene_begin to $gene_end\tNew Coords: $new_gene_begin to $new_gene_end\n";

		}
		if ($new_gene_begin){$gene_begin = $new_gene_begin};	


		my $gene = &gjoseqlib::DNA_subseq($contigH{$contig}, $gene_begin, $gene_end );	
		my $protein = &gjoseqlib::translate_seq( $gene );

		if ($start_to_met){$protein =~ s/^[A-Z]/M/i;}   # M not m: /i affects the pattern, not the replacement (gist #5)
	
		# Set up calling non-pssm features that are anchored to pssm coordinates
		if (exists $options->{$virus}->{features}->{$pssmdir}->{non_pssm_partner})
		{
			my @non_pssms = @{$options->{$virus}->{features}->{$pssmdir}->{non_pssm_partner}};
			
			for my $i (0..$#non_pssms)
			{
				my $start_coord;
				my $stop_coord;
				my $feat = $non_pssms[$i];
				#if the current pssm match feature corresponds with the non-pssm start site
				if ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{START_OFFSET} = $options->{$virus}->{features}->{$feat}->{begin}->{begin_offset};
					if ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm_loc} =~ /START/)
					{	
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_begin;
					}
					elsif ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_end;
					}
				}
				elsif ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{STOP_OFFSET} = $options->{$virus}->{features}->{$feat}->{end}->{end_offset};
					if ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm_loc} =~ /START/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_begin;
					}
					elsif ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_end;
					}
				}			
				$non_pssm_feat->{$feat}->{ANNO}     = $options->{$virus}->{features}->{$feat}->{anno};
				$non_pssm_feat->{$feat}->{MIN}      = $options->{$virus}->{features}->{$feat}->{min_len};
				$non_pssm_feat->{$feat}->{MAX}      = $options->{$virus}->{features}->{$feat}->{max_len};			
				$non_pssm_feat->{$feat}->{AA}       = $options->{$virus}->{features}->{$feat}->{translate};
				$non_pssm_feat->{$feat}->{TYPE}     = $options->{$virus}->{features}->{$feat}->{feature_type};
				$non_pssm_feat->{$feat}->{SYMBOL}   = $options->{$virus}->{features}->{$feat}->{gene_symbol};
			}
		}
			push @all_seqs, ([$best_results->[$i]->{contig}, $gene_begin, $gene_end, $anno, $strand, $best_pssm, $gene, $protein, $feature_type, $symbol]); 
	}
	print STDERR "-----------------------\n"; 
}


#add any non-pssm features that are anchored to PSSM coordinates
if ($non_pssm_feat)
{
	my @tuples = call_non_pssm_features($non_pssm_feat, \%contigH);   # pass a ref, not a flattened hash (gist #13)
	push @all_seqs, @tuples;
}


#Tuple is:	
# 0 contig
# 1 start
# 2 end
# 3 anno
# 4 strand
# 5 pssm  (this is the full pssmfile name)
# 6 gene
# 7 protein
# 8 feature_type	
# 9 symbol

# $virus is the best virus match
# Sort the output in order of contig and then start position. 

my (@prot_seqs, @gene_seqs);
my $count = 0;
foreach (@contig_order)
{
	my $contig = $_;
	my @features = ();
	for my $i (0..$#all_seqs)
	{
		if ($all_seqs[$i][0] eq $contig)
		{
			push @features, $all_seqs[$i];		
		}		
	}
	my @sorted =  sort { $a->[1] <=> $b->[1] } @features;
	foreach(@sorted)
	{
		$count ++; 
		# this is just a silly place holder unitl bob can help me get set up with official
		# ID generation
		my $prot_id = "lv\|$tax\.$version\.$_->[8]\.$count"; 
		
		#Trimming the star from the end of the protein sequence. 
		#This could be removed if its unwanted.
		my $protein = $_->[7];
		$protein =~ s/\*$//g;

		push @gene_seqs, ([$prot_id, $_->[3], $_->[6]]);
		push @prot_seqs, ([$prot_id, $_->[3], $protein]);
		
		# Original Feature Table looks like this:
		# genome_id	genome_name	accession	annotation_source	feature_type	patric_id	refseq_locus_tag	start	end	strand	na_length	gene	product	plfam_id	pgfam_id
		# For now, I will keep:
		# genome_id	genome_name	accession(contig)	annotation_source	feature_type	madeupid	Gene_symbol	start	end	strand	na_length	Fam	Closest_Genome	best_contig_bit pssm	product	NT_Seq	AA_seq


		my $na_len = length $_->[6];
		unless($no_out)
		{
			if ($append_seqs)
			{
				print TBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";
			}
			else
			{
				print TBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";
			}
		}
		if($tbl_only)
		{
			if ($append_seqs)
			{
				print "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";
			}
			else
			{
				print "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";
			}		
		}
		if($ctbl)
		{
			if ($append_seqs)
			{
				print CTBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";

			}
			else
			{
				print CTBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$anno_rep\t$anno_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";

			}
		}
	}
}




unless($no_out)
{
	&gjoseqlib::print_alignment_as_fasta(\*DNA, @gene_seqs);
	&gjoseqlib::print_alignment_as_fasta(\*AA, @prot_seqs);
}
if ($dna_only)
{
	&gjoseqlib::print_alignment_as_fasta(\@gene_seqs);
}
if ($aa_only)
{
	&gjoseqlib::print_alignment_as_fasta(\@prot_seqs);
}


chdir ($base);
unless ($keep_temp){system "rm -rf $tmp";}

# The taxon was decided and every PSSM was searched; $count is how many features survived
# their bit_cutoff.  Zero is a legitimate outcome, not an error -- for a short fragment
# "nothing to annotate" is the right answer -- but it is a different outcome from a genome
# with features, and a caller running the downstream stages needs to tell them apart:
# transcript-editing, splice variants and quality all die on a GTO with no features.
exit($count ? EXIT_OK : EXIT_NO_FEATURES);



##############sub call_non_pssm_features####################
#   Calles non pssm-based features that can be called based on 
#   protein coordinates.
#
#   call_non_pssm_features ($non_pssm_feat, $contigH)
#
#-----------------------------------------------------------
sub call_non_pssm_features
{
	my ($featH, $contigH) = @_;
	my @seq_data;
		
	foreach (keys %$featH)
	{
		my $feat          = $_;
		my $anno          = $featH->{$feat}->{ANNO};
		my $min           = $featH->{$feat}->{MIN};
		my $max           = $featH->{$feat}->{MAX};
		my $aa            = $featH->{$feat}->{AA};
		my $start_offset  = $featH->{$feat}->{START_OFFSET};
		my $stop_offset   = $featH->{$feat}->{STOP_OFFSET};
		my $feature_type  = $featH->{$feat}->{TYPE};
		my $symbol        = $featH->{$feat}->{SYMBOL};


		foreach (keys %{$featH->{$feat}->{COORD}})
		{
			my $contig = $_;
			
			if (($featH->{$feat}->{COORD}->{$contig}->{START}) && ($featH->{$feat}->{COORD}->{$contig}->{STOP}))
			{	
				my @starts = @{$featH->{$feat}->{COORD}->{$contig}->{START}};
				my @stops  = @{$featH->{$feat}->{COORD}->{$contig}->{STOP}};			
			
				for my $i (0..$#starts)
				{			
					my $start = $starts[$i]; 
					for my $j (0..$#stops)
					{
						my $stop = $stops[$j];
						if ($stop)
						{
							my ($strand, $begin, $end);
							if ($start < $stop)
							{
								$strand = "+";
								$begin = $start + $start_offset;   # was += : mutated $start across STOP anchors (gist #7)
								$end   = $stop  - $stop_offset;													
							}
							elsif ($start > $stop)
							{
								$strand = "-";
								$begin = $start - $start_offset;
								$end   = $stop  + $stop_offset;						
							}
					
					
							if ((abs($begin - $end) <= $max) && (abs($begin - $end) >= $min))
							{
								print STDERR "\tCalling non-PSSM feature\t$anno\tbegin: $begin\tend: $end\n"; 
								my $nt = &gjoseqlib::DNA_subseq($contigH->{$contig}, $begin, $end); 
								my $prot;
								if ($aa){ $prot = &gjoseqlib::translate_seq( $nt );}
								if ($prot =~ /\*(?!$)/)# It will not record a position-called protein with stops
								{
									print STDERR "\tInternal stop(s) found in non-pssm feature. No assignment made for $anno\n";
								}
								else
								{
									push @seq_data, ([$contig, $begin, $end, $anno, $strand, $feat, $nt, $prot, $feature_type, $symbol]);
								}
							}
						}
					}
				}
			}
		}
	}
	return @seq_data;
}













##########################sub crop_to_stop_codon###########
#
# Returns new gene boundaries, ignoring everything after the stop codon.
# Returns gene boundaries with stop codon by convention
#
#  $gene_end = crop_to_stop_codon($gene_start, $gene_end ,$hseq)
#  hseq is the AA sequence blast match containing the stop character (*)
#
#----------------------------------------------------------
sub crop_to_stop_codon
{
	my ($from, $to, $hseq) = @_; 		
	my $gene_end;
	$hseq =~ s/\*.+//g;
	my $len = length ($hseq); 
	my $to_keep = ($len * 3); 
	if ($from < $to){$gene_end = ($from + (($to_keep + 3) - 1))}
	if ($from > $to){$gene_end = (($from - ($to_keep) - 3) + 1) }
	return $gene_end;
}
###########################################################



##########################sub scan_to_stop_codon###########
# Looks beyond the end of the called blast boundaries to find the next stop codon
# returns the new gene coordinates, which includes the stop codon by convention
# The contig is the DNA string.  It is needed to ensure that you don't fall off of either end.
#
#   ($gene_begin, $gene_end)= scan_to_stop_codon ($to, $from, $contig);
#
#----------------------------------------------------------

sub scan_to_stop_codon
{
	my ($from, $to, $contig) = @_;
	my $len = length $contig; 
	my $end = $to;
		
	if ($from < $to)
	{		
		for (my $i = $end; $i <= ($len - 3); $i += 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i +1), ($i + 3) ); #<--- these are equivalent.
			my $aa = &gjoseqlib::translate_codon( $codon ); 
			
			if ($aa =~ /\*/) 
			{
				$end = ($i + 3);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			elsif ($aa =~ /x/i) 
			{
				print STDERR "\tC-term Extension X found was: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$end = ($i + 3);
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
			}		
		}
	}

	## Reverse:

	elsif ($from > $to)
	{		
		for (my $i = $end; $i > 3; $i -= 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i -1), $i - 3); #<--- these are equivalent.
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /\*/) 
			{
				$end = ($i - 3);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$aa\t$codon\n";

				last;	
			}		
			elsif ($aa =~ /x/i) 
			{
				print STDERR "\tC-term Extension X found was: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$end = ($i - 3);
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
			}		
		}
	}	
	return ($from, $end);
}
###########################################################


##########################sub scan_to_met_start#############
#
# Looks upstream for a new Met start codon.
# Returns the start position for the new Met codon.
# Does not consider alternative start codons.
#  ($new_gene_begin, $new_gene_end) = scan_to_met_start($gene_start, $gene_end, $contig)
# 
#-----------------------------------------------------------
sub scan_to_met_start
{
	my ($from, $to, $contig) = @_;
	my $len = length $contig; 
	my $start = $from;

	if ($from < $to)
	{		
		for (my $i = ($start - 3); $i >= 0; $i -= 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, $i, ($i + 2)); 
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /m/i) 
			{
				$start = $i;  #move the end position to the end of the stop codon and quit.
				print STDERR "\tN-term Extension Met found: was: $from\tnow: $start\t$aa\t$codon\n";
				last;	
			}		
			elsif ($aa =~ /x|\*/i) 
			{
				print STDERR "\tN-term Extension stopping was: $from\tnow: $start\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$start = $i;
				print STDERR "\tN-term Extended\twas: $from\tnow: $start\t$codon\t$aa\n";
			}		
		}
	}
	
	### Reverse:
	if ($from > $to)
	{		
		
		for (my $i = ($start + 1); $i <= ($len - 2); $i += 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i + 2), $i); 
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /m/i) 
			{
				$start = ($i + 2);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tN-term Extension Met found: was: $from\tnow: $start\t$aa\t$codon\n";
				last;	
			}		
			elsif ($aa =~ /x|\*/i) 
			{
				print STDERR "\tN-term Extension stopping was: $from\tnow: $start\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$start = ($i + 2);
				print STDERR "\tN-term Extended\twas: $from\tnow: $start\t$codon\t$aa\n";
			}		
		}
	}	
	return ($start, $to);
}
###########################################################



##########################sub matching_tblastn_hsps_json#####
# Find the best hsp from a tblastn blast result in json format
  
# Start with this:
#	use JSON::XS;
#   open (IN, "tblastn -outfmt 13 -db $s_file -in_pssm $pssm_file |");
#	my $blast = decode_json(scalar read_file(\*IN));	
#
#
#  (@results, $best_bit) = matching_tblastn_hsps_json($blast, $bitscore_cutoff, $coverage_cutoff);
#
#	where @results is a hash reference of every hsp meeting the bit score and coverage cutoff
#
#----------------------------------------------------------
sub matching_tblastn_hsps_json
{
	my ($blast, $bit_cutoff, $cov_cutoff) = @_;
	my $results = {};
	my @data;
	my $best_bit = 0;
	
	my @output = @{$blast->{BlastOutput2}};   # deref: was a 1-elem list (gist #12)
	for my $i (0..$#output)
	{
		my $iterations = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations};	
		my $nitr = scalar @$iterations;

		for my $j (0..($nitr -1 ))
		{
			if ($blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{iter_num} == 1)   # == not = (gist #11)
			{
				my $hits = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits};
				my $nhits = scalar @$hits;
			
				for my $k (0..($nhits -1))
				{
					# okay for reference, if there are two contigs, the hsps will be split into two different "hits"
				
					# get the contig ID.  If i understand this correctly, the title has to be the same per HSP, so the hardcoded zero should be ok.
					my $contig = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{description}->[0]->{title};
					$contig =~ s/\s.+//g;								

					my $hsps = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps};
					my $nhsps = scalar @$hsps;
															
					for my $l (0..($nhsps -1))
					{
						$results->{contig}    = $contig;
						$results->{bit}       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{bit_score};						
						$results->{hseq}      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hseq};
						$results->{hit_from}  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_from};
						$results->{hit_to}    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_to};
						$results->{frame}     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_frame};
						$results->{q_from}    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{query_from};
						$results->{q_to}      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{query_to};
						$results->{ali_len}   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{align_len};   # was {query_to} (copy-paste); field is currently unused
						$results->{e_val}     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{evalue};
						$results->{hsp_num}   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{num};
				
						my $cov = ((length $results->{hseq})/(abs($results->{q_to} - $results->{q_from})+1));	
						$results->{cov} = $cov;

						if (($results->{bit} > $bit_cutoff) && ($cov > $cov_cutoff))
						{	
							push @data, $results;
						}
					
						if (($results->{bit} > $best_bit) && ($results->{bit} > $bit_cutoff))
						{
							$best_bit = $results->{bit};
						}
						$results = {};
					}
				}
			}
		}
	}
	return (\@data, $best_bit);
}
###########################################################



##########################sub get_blastn_bit_##############
# Dig the bit score out of the json format blastn
  
# Start with something like this:
#	use JSON::XS;
#		open (IN, "blastn -query $rep_file -subject $s_file -evalue 0.5 -reward 2 -penalty -3 -word_size 11 -outfmt 13 -soft_masking false |"); 
#		my $blastn = decode_json(scalar read_file(\*IN));	
#
#
#  $bit_score = get_blastn_bit($json);
#----------------------------------------------------------
sub get_blastn_bit 
{
    my $blast = shift @_;
    my $best_bit = 0;
    
    # Check if BlastOutput2 exists and is an array reference
    return 0 unless $blast->{BlastOutput2} && ref($blast->{BlastOutput2}) eq 'ARRAY';
    
    # Iterate through each BlastOutput2 element
    for my $output (@{$blast->{BlastOutput2}}) 
    {
        next unless $output->{report} && $output->{report}->{results};
        
        my $results = $output->{report}->{results};
        next unless $results->{bl2seq} && ref($results->{bl2seq}) eq 'ARRAY';
        
        # Iterate through bl2seq array
        for my $bl2seq_item (@{$results->{bl2seq}}) 
        {
            next unless $bl2seq_item->{hits} && ref($bl2seq_item->{hits}) eq 'ARRAY';
            next if @{$bl2seq_item->{hits}} == 0; # Skip empty hits arrays
            
            # Iterate through hits
            for my $hit (@{$bl2seq_item->{hits}}) 
            {
                next unless $hit->{hsps} && ref($hit->{hsps}) eq 'ARRAY';
                
                # Iterate through hsps
                for my $hsp (@{$hit->{hsps}}) 
                {
                    if ($hsp->{bit_score} && $hsp->{bit_score} > $best_bit) 
                    {
                        $best_bit = $hsp->{bit_score};
                    }
                }
            }
        }
    }
    
    return $best_bit;
}
###########################################################





##########################sub valid_vtax_list##############
# The canonical list of annotation taxon names accepted by -vtax.
#
# Sourced from the intersection of the PSSM directory names and the top-level
# keys of the JSON options file, because those are the only two things the
# chosen taxon is ever used to index.  Deriving it from -pssm and -j means a
# user who overrides either one automatically gets the matching list.
# Viral-Rep-Contigs is deliberately not consulted: a rep file with no JSON
# entry cannot be annotated with anyway.
#
#   my @names = valid_vtax_list($pssm_dir, $options);
#
#----------------------------------------------------------
sub valid_vtax_list
{
	my ($pdir, $options) = @_;

	opendir (my $dh, $pdir) or bail(EXIT_USAGE, "Cannot open PSSM directory $pdir: $!\n");
	my @dirs = sort grep{/\.pssms$/} readdir($dh);
	closedir($dh);

	my @names;
	foreach (@dirs)
	{
		my $name = $_;
		$name =~ s/\.pssms$//;
		push @names, $name if exists $options->{$name};
	}
	return @names;
}

##########################sub taxon_family_map##############
# taxon => family, read from the "family" field of each annotation taxon in the JSON
# options file.  This is the only place family membership is recorded: the 34 taxon names
# are a mix of ranks (8 families, 25 genera, 1 species), so it cannot be derived from the
# names, and the reference filenames carry the taxon but not its parent.
#
# Eight taxa are their own family -- the Bunyavirales entries, where LowVan resolves no
# deeper -- so a family may map to exactly one taxon.  That is not an error; -vfam on such
# a family is simply equivalent to -vtax on the taxon.
#
# Dies on a missing field rather than defaulting.  A taxon with no family would silently
# drop out of every -vfam sweep, which is the sort of thing that shows up as an
# unexplained classification failure months later.
#
#   my %family = taxon_family_map($options, \@valid_vtax);
#
#----------------------------------------------------------
sub taxon_family_map
{
	my ($options, $taxa) = @_;

	my %map;
	my @missing;
	foreach my $t (@$taxa)
	{
		my $f = $options->{$t}->{family};
		if (defined $f && $f ne "") { $map{$t} = $f }
		else                        { push @missing, $t }
	}

	if (@missing)
	{
		bail(EXIT_USAGE,
		    "No \"family\" field for these annotation taxa in the JSON options file:\n"
		  . join("", map {"\t$_\n"} @missing)
		  . "Every taxon needs one; -vfam and -list-vfam are built from them.\n");
	}
	return %map;
}
###########################################################
