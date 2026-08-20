#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use File::Temp;
use Getopt::Long::Descriptive;
use File::Copy;
use IPC::Run qw(run);
use Cwd;

# Try to load version module; fall back to "dev" if not available
my $tool_version;
eval {
    require LowVanVersion;
    $tool_version = LowVanVersion::get_version();
};
if ($@ || !$tool_version) {
    $tool_version = "dev";
}

my $default_data_dir = $ENV{LOWVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

#
# The base script's exit codes, mirrored here so this wrapper can read a result instead of
# grepping English out of stderr.  Keep in step with annotate_by_viral_pssm.pl, which is
# where they are defined and documented.
#
# 10 and up are the "nothing was annotated" outcomes.  All three are legitimate results for
# a short or divergent fragment rather than errors, which is exactly the distinction the
# caller needs: the downstream stages all die on a GTO with no features, so a pipeline
# should skip them for 10/11/12 rather than fail the job.
#
use constant EXIT_OK            => 0;
use constant EXIT_USAGE         => 1;
use constant EXIT_INTERNAL      => 2;
use constant EXIT_NO_REFERENCE  => 10;
use constant EXIT_NO_FEATURES   => 11;
use constant EXIT_OUT_OF_BOUNDS => 12;

# rc -> the metadata "status" string recorded on the analysis event.  Anything not listed
# (an uncaught die, 255, a signal) is reported as "error": unexpected, not a known outcome.
my %STATUS_FOR_RC = (
    EXIT_OK()            => "annotated",
    EXIT_USAGE()         => "usage_error",
    EXIT_INTERNAL()      => "internal_error",
    EXIT_NO_REFERENCE()  => "no_reference",
    EXIT_NO_FEATURES()   => "no_features",
    EXIT_OUT_OF_BOUNDS() => "out_of_bounds",
);

my($opt, $usage) = describe_options("%c %o",
				    ["input|i=s"       => "Input file"],
				    ["output|o=s"      => "Output file"],
				    ["prefix|x=s"      => "File Prefix", { default => "Viral_Anno" }],
				    ["remove-existing" => "Remove existing CDS, mat_peptide, and RNA features if run is successful"],	
				    ["threads|t=i"     => "Limit to this many threads", { default => 8 }],
				    ["debug|d"         => "Enable debugging"],
				    ["cdir|c=s"        => "Full path to reference contigs directory", {default => "$default_data_dir/Viral-Rep-Contigs"}],
				    ["pdir|p=s"        => "Full path to the PSSM directory", {default => "$default_data_dir/Viral-PSSMs"}],
				    ["json|j=s"        => "Full path to the JSON opts file", {default => "$default_data_dir/Viral_PSSM.json"}],
				    ["max|a=i"         => "Max contig length, default is 40000", { default => 40000 }],
				    ["min|z=i"         => "Min contig length, default is 300", { default => 300 }],
				    # No default here on purpose: leave it to the base script (-mcb 150) rather
				    # than pinning a second copy of the number that could drift out of step.
				    ["min-contig-bit=f" => "Minimum classification BLASTn bit score to accept a taxon call, passed through as -mcb (base script default 150). Chiefly for --viral-family, where restricting the sweep to one family makes a call in the 50-150 band far more trustworthy than an unrestricted one"],
				    ["viral-taxon=s"   => "Declare the annotation taxon instead of detecting it by BLASTn. See --list-viral-taxa. Also sets viral_family, which selects the Transcript-Editing/ and Splice-Variants/ dirs used by the downstream steps"],
				    ["list-viral-taxa" => "List the valid --viral-taxon values and exit"],
				    ["viral-family=s"  => "Declare the annotation family when the genus is not known, and let BLASTn choose the genus within it. For records whose lineage stops at family. Unlike --viral-taxon this still makes a call, so -mcb applies and a genome that clears nothing is still rejected. See --list-viral-families. Cannot be combined with --viral-taxon or --skip-classification"],
				    ["list-viral-families" => "List the valid --viral-family values and exit"],
				    ["fallback-viral-taxon" => "If the BLASTn classification finds no reference above -mcb, retry once with --viral-taxon derived from the lineage in the input GTO's own \"taxonomy\" field. BLASTn still gets first refusal, so a genome it can place is annotated exactly as it is today; the declared taxon is used only where there would otherwise be no annotation at all. Does nothing if the GTO has no lineage, or if the lineage resolves no deeper than a rank LowVan has no taxon for"],
				    ["skip-classification" => "With --viral-taxon, BLASTn only that taxon reference contigs (1-14) instead of every reference contig. Much faster, but no cross-check against the other taxa, so a wrong --viral-taxon will not be caught. For bulk reruns where the taxon is already known to be correct"],
				    ["margin=f"        => "Minimum ratio between the winning and runner-up taxon bit scores, e.g. --margin 1.3. Off by default. Warns when the classification is a near-tie and records the margin on the GTO close_genomes record; never rejects. Incompatible with --skip-classification, which searches one taxon only"],
				    ["propagate-exit-status" => "Exit with the base script's status (10 = no reference above -mcb, 11 = no feature cleared its bit_cutoff, 12 = length out of bounds) instead of always exiting 0. The output GTO is still written either way. Off by default because a caller that runs the annotation and the downstream stages as one pipeline would fail the whole job on a legitimately empty annotation; turn it on once the caller decides per stage. The same information is always on the analysis event, as success and metadata.status"],
				    ["version|v"       => "Show version information"],
				    ["help|h"          => "Show this help message"]);

if ($opt->version) {
    print "annotate_by_viral_pssm-GTO.pl version $tool_version\n";
    exit 0;
}
print($usage->text), exit 0 if $opt->help;

# Delegate to the base script rather than reimplementing, so there is one canonical list.
# Must run before create_from_file() below, which requires -i.
if ($opt->list_viral_taxa) {
    my $listed = run(["annotate_by_viral_pssm.pl", "-list-vtax",
		      "-pssm", $opt->pdir, "-j", $opt->json]);
    exit($listed ? 0 : 1);
}
if ($opt->list_viral_families) {
    my $listed = run(["annotate_by_viral_pssm.pl", "-list-vfam",
		      "-pssm", $opt->pdir, "-j", $opt->json]);
    exit($listed ? 0 : 1);
}

die($usage->text) if @ARGV != 0;

# --viral-taxon already declares the answer, and the base script downgrades the -mcb
# rejection to a warning under it, so there is no classification failure to fall back
# from.  --viral-family is allowed: it CAN fail -mcb, and if the lineage turns out to
# carry a genus after all, retrying with that genus is strictly better than giving up.
if ($opt->fallback_viral_taxon && $opt->viral_taxon)
{
    die "--fallback-viral-taxon and --viral-taxon cannot be combined: the taxon is already declared,\n"
      . "and under --viral-taxon a low -mcb is a warning rather than a failure, so nothing would trigger.\n";
}

chomp(my $hostname = `hostname`);

my $tempdir = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1));
my $prefix = $opt->prefix // "Viral_Anno";

print STDERR "Tempdir=$tempdir\n" if $opt->debug;

my $here = getcwd;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";



#
# Invoke annotate_by_viral_pssm.pl to annotate viral genome.
#
# We parse the feature table file that is generated by annotate_by_viral_pssm.pl. It contains three feature types.
# CDS, RNA, and mat_peptide
# It returns a GTO to stdout
#

my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
my $taxon_id       = $genome_in->{ncbi_taxonomy_id};   #I should be able to get rid of this if I don't require it in the base-script
my $name           = $genome_in->{scientific_name};        #I should also be able to get rid of this 

if (! $sequences_file){die "No sequences in the input GTO\n";}
if (! $taxon_id)      {die "No NCBI taxonomy ID in the input GTO\n";}
if (! $name)          {die "No genome name in the input GTO\n";}

my @params = ("-i",    $sequences_file,
		      "-t",    "$tempdir",
		      "-p",    $prefix,
		      "-tax",  $taxon_id,
		      "-g",    $name,
		      "-threads", $opt->threads,
		      "-c",    $opt->cdir,
		      "-pssm", $opt->pdir,
		      "-j",    $opt->json,
		      "-min",  $opt->min,
		      "-max",  $opt->max,
		      "-s",
		      "-no",
		      "-tbl",
		      "-tmp");

push @params, ("-vtax", $opt->viral_taxon)         if $opt->viral_taxon;
push @params, ("-vfam", $opt->viral_family)        if $opt->viral_family;
push @params, ("-mcb",  $opt->min_contig_bit)      if defined $opt->min_contig_bit;
push @params, ("-skip-classification")             if $opt->skip_classification;
push @params, ("-margin", $opt->margin)            if $opt->margin;

print STDERR Dumper(\@params);
my $ok = run(["annotate_by_viral_pssm.pl", @params], ">", "$here/$prefix.stdout.txt", "2>", "$here/$prefix.stderr.txt");
# Capture immediately: $? is clobbered by the next system()/backtick/run().
my $rc = $? >> 8;
my $retried = 0;

# --- fallback: BLASTn could not place the genome, so declare the taxon ourselves ------
#
# Order matters and is the point of the option: BLASTn gets first refusal, so any genome
# it can place is annotated exactly as it is without the flag.  The lineage is consulted
# only for genomes that would otherwise produce no annotation at all.
#
# The trigger is a *rejection*, not merely a failure: retrying a BLAST crash under a
# declared taxon would just crash again while attaching a provenance record claiming the
# taxon was chosen on purpose.  rc = EXIT_NO_REFERENCE now says so directly; the stderr
# match is kept as a fallback for a deployed base script predating the exit codes.
my ($fallback_taxon, $fallback_rank);
if ($opt->fallback_viral_taxon && !$ok && classification_rejected($rc, "$here/$prefix.stderr.txt"))
{
    my @valid = valid_viral_taxa($opt->pdir, $opt->json);
    ($fallback_taxon, $fallback_rank) = taxon_from_lineage($genome_in, \@valid);

    if (!defined $fallback_taxon)
    {
        my $lin = $genome_in->{taxonomy};
        print STDERR "--fallback-viral-taxon: classification found no reference above the cutoff, but "
                   . (defined $lin && $lin ne ""
                        ? "no element of the GTO lineage is a LowVan annotation taxon.\n\tlineage: $lin\n"
                        : "the input GTO has no \"taxonomy\" field to fall back to.\n");
    }
    else
    {
        print STDERR "--fallback-viral-taxon: classification found no reference above the cutoff; "
                   . "retrying with -vtax $fallback_taxon (from GTO lineage element '$fallback_rank').\n";

        # Drop -vfam if it was given: we now have a genus, which is strictly more specific
        # than the family scope, and the two cannot be passed together.
        my @retry;
        for (my $i = 0; $i <= $#params; $i++)
        {
            if ($params[$i] eq "-vfam") { $i++; next }
            push @retry, $params[$i];
        }
        push @retry, ("-vtax", $fallback_taxon);
        @params = @retry;

        # Keep the BLASTn attempt's stderr; the retry would otherwise overwrite the only
        # record of why the classification was rejected.
        move("$here/$prefix.stderr.txt", "$here/$prefix.stderr-classification.txt");

        print STDERR Dumper(\@params);
        $ok = run(["annotate_by_viral_pssm.pl", @params],
                  ">", "$here/$prefix.stdout.txt", "2>", "$here/$prefix.stderr.txt");
        $rc = $? >> 8;
        $retried = 1;
    }
}

my $status = $STATUS_FOR_RC{$rc} // "error";

if (!$ok)
{
    print STDERR "Viral Annotation run failed with rc=$rc ($status). Stderr:\n";
    copy("$here/$prefix.stderr.txt", \*STDERR);   # was hardcoded Viral_Anno.stderr.txt + mislabeled "Stdout" (gist #17)
    # Execution continues below and still writes an output GTO -- rc 10/11/12 are legitimate
    # "nothing to annotate" outcomes, not errors, and the GTO the caller gets back is the
    # input plus an analysis event saying so.  Whether that becomes a non-zero exit for the
    # caller is --propagate-exit-status, decided at the bottom of the script.
}



my $event = {
    tool_name => "LowVan Annotate",
    tool_version => $tool_version,
    execution_time => scalar gettimeofday,
    parameters => \@params,
    hostname => $hostname,
};

# The event is added before the features are parsed because each feature has to reference
# its id.  add_analysis_event stores the hashref we hand it rather than a copy, so success
# and metadata are filled in below, once there is something to report, and the change is
# reflected in the GTO that gets written.
my $event_id = $genome_in->add_analysis_event($event);

my %meta = (
    status    => $status,
    exit_code => $rc,
);
$meta{retried_with_declared_taxon} = 1 if $retried;
$meta{min_contig_bit} = $opt->min_contig_bit if defined $opt->min_contig_bit;
$meta{margin_threshold} = $opt->margin if $opt->margin;

# The -min/-max gate excluded 94.1% of the reannotation failure set, so when the answer is
# out_of_bounds the length is the whole explanation and is worth having on the record.
if (ref($genome_in->{contigs}) eq 'ARRAY')
{
    my $bp = 0;
    $bp += length($_->{dna} // "") foreach @{$genome_in->{contigs}};
    $meta{input_contigs} = scalar @{$genome_in->{contigs}};
    $meta{input_bp}      = $bp;
}


## Parse the generated peptide file. We collect the CDS, mature_peptides, and RNAs then
## add features so that we can register the counts.
##

my $viral_family;
my $features_added = 0;

my ($close_bit, $close_id, $close_name, $close_file);
if (open(my $tbl, "<", "$here/$prefix.stdout.txt"))
{
	my %features;
	while (<$tbl>)
	{
		chomp;
		my ($local_genome_id, $name, $contig, $anno_source, $type, $local_peg_id, $symbol, $start, $stop, $strand, $len, $virus, $cf, $cb, $ci, $cn, $pssm, $anno, $dna, $aa) = split /\t/; 
		$close_file = $cf;
		$close_bit  = $cb;
		$close_id   = $ci;
		$close_name = $cn;

		$viral_family = $virus;
		
		my $feature;
		if ($type =~ /(mat_peptide)|(CDS)/)
		{
			my $pssm_id = $pssm;
			$pssm_id =~ s/\.pssm$//;  # Remove .pssm suffix for family_assignments
			$feature =
			{
				type        => $type,
				contig      => $contig,
				aa_sequence => $aa,
				location    => ([[$contig, $start, $strand, $len]]),
				product     => $anno,
				pssm        => ([["LOWVAN", $pssm_id, $anno, "LowVan Annotate $tool_version"]]),
				symbol      => $symbol,
			}
		}
		elsif ($type =~ /RNA/)
		{
			$feature =
			{ 
				type        => $type,
				contig      => $contig,
				location    => ([[$contig, $start, $strand, $len]]),
				product     => $anno,
				symbol      => $symbol,
			}
		}		
		push(@{$features{$type}}, $feature);
	}
	
    if (%features && $opt->remove_existing)
    {
		my @to_del = $genome_in->fids_of_type('CDS', 'mat_peptide', 'RNA');
		print STDERR "Delete @to_del\n";
		$genome_in->delete_feature($_) foreach @to_del;
		$meta{removed_existing_features} = scalar @to_del;
    }

    
    for my $type (sort keys %features)   # sort: deterministic new_feature_id assignment across runs (gist #19)
    {
		my $feats = $features{$type};
		my $n = @$feats;
		my $id_type = $type;
	
		for my $feature (@$feats)
		{
			my $p = {
					-id	                 => $genome_in->new_feature_id($id_type),
					-type 	             => $type,
					-location 	         => $feature->{location},
					-analysis_event_id 	 => $event_id,
					-annotator           => 'LowVan Annotate',
					#-alias_pairs         => [[gene => $feature->{symbol}]],
					-protein_translation => $feature->{aa_sequence},
					-function            => $feature->{product},
					-family_assignments  => $feature->{pssm},
					};
			
			if (defined $feature->{symbol} && $feature->{symbol} ne '') 
			{
				$p->{-alias_pairs} = [[gene => $feature->{symbol}]];
   			}
			$genome_in->add_feature($p);
			$features_added++;
		}
		$meta{"features_$type"} = $n;
    }

	#add close genome:
	#$close_genome, $close_bit, $close_id, $close_name,
	# An empty feature table leaves all of these undef -- do not append a record of nulls.
	# Much more likely now that --viral-taxon can name a taxon whose PSSMs match nothing.
	# Read once.  The sidecar carries the margin fields for close_genomes and, under
	# --viral-family, the genus BLASTn actually chose -- which is needed below even when
	# there is no close genome record to hang it on.
	my $cls = read_classification("$here/$prefix.classification");

	if (defined $close_id)
	{
		my $close = { genome_id => $close_id, genome_name => $close_name, file_name => $close_file, closeness_measure => "BLASTn bit score", closeness_value => $close_bit, analysis_method => "LowVan Annotate"};

		# --margin: carry the confidence of the classification, not just its strength.
		# The base script wrote these to $prefix.classification; absent that file the
		# flag was not used and the record keeps its original shape.
		if ($cls)
		{
			$close->{margin}          = $cls->{margin}        if defined $cls->{margin};
			$close->{runner_up}       = $cls->{runner_up}     if defined $cls->{runner_up};
			$close->{runner_up_value} = $cls->{runner_up_bit} if defined $cls->{runner_up_bit};
			$close->{margin_below_threshold} = $cls->{below_threshold}
				if defined $cls->{below_threshold};

			# --viral-family: the genus here was picked from a deliberately narrowed
			# field, often with -mcb relaxed.  Say so on the record, so a consumer can
			# tell it apart from a call made against every reference in the set.
			if (defined $cls->{scope} && $cls->{scope} eq "family")
			{
				$close->{classification_scope}  = "family";
				$close->{classification_family} = $cls->{family} if defined $cls->{family};
			}
		}

		# --fallback-viral-taxon: this genome was NOT placed by BLASTn -- the taxon came
		# from its own lineage after the classification was rejected.  Say so, and say
		# which lineage element it came from: the closeness_value below is then the best
		# reference *within* the declared taxon, not the score that chose it.
		if (defined $fallback_taxon)
		{
			$close->{classification_scope}   = "lineage_fallback";
			$close->{classification_taxon}   = $fallback_taxon;
			$close->{classification_lineage} = $fallback_rank;
		}

		push(@{$genome_in->{close_genomes}}, $close);

		$meta{closest_reference}      = $close_file if defined $close_file;
		$meta{closest_reference_bit}  = $close_bit  if defined $close_bit;
		$meta{runner_up}              = $close->{runner_up}       if defined $close->{runner_up};
		$meta{runner_up_bit}          = $close->{runner_up_value}  if defined $close->{runner_up_value};
		$meta{margin}                 = $close->{margin}          if defined $close->{margin};
		$meta{margin_below_threshold} = $close->{margin_below_threshold}
			if defined $close->{margin_below_threshold};
	}

	# How the taxon was arrived at, which is not recoverable from the taxon name alone and
	# is the first thing anyone reviewing a questionable annotation wants to know.
	$meta{classification_scope} = $opt->viral_taxon    ? "declared_taxon"
	                            : defined $fallback_taxon ? "lineage_fallback"
	                            : $opt->viral_family   ? "declared_family"
	                            :                        "blastn";
	$meta{classification_lineage_element} = $fallback_rank if defined $fallback_rank;
	$meta{declared_family} = $opt->viral_family if $opt->viral_family;
	$meta{skip_classification} = 1 if $opt->skip_classification;

	# Fall back to the declared taxon: without this, a --viral-taxon run that calls zero
	# features leaves viral_family undef and every downstream step dies on its "$fam or die".
	#
	# --viral-family has no declared taxon to fall back to -- the genus was chosen inside
	# the base script -- so take it from the sidecar, which -vfam always writes.  Without
	# this an empty annotation under --viral-family reintroduces exactly the guard-2 death
	# that declaring a taxon exists to prevent.
	# --fallback-viral-taxon lands in the same hole as --viral-taxon and for the same
	# reason: the retry may call zero features, and then nothing else here knows what the
	# genome was annotated as.  $fallback_taxon is what -vtax was given on the retry.
	$genome_in->{viral_family} = $viral_family
	                          // ($cls ? $cls->{annotated_as} : undef)
	                          // $opt->viral_taxon
	                          // $fallback_taxon;
}
else
{
    warn "Could not read $here/$prefix.stdout.txt\n";
    $meta{status} = $status = "error" if $status eq "annotated";
    $meta{feature_table_unreadable} = 1;
}

# "Success" is the question the caller actually has: did this run produce an annotation?
# Not "did the process exit 0" -- rc 11 is an orderly exit that annotated nothing -- and not
# "does the GTO have features", which is true of any GTO that arrived carrying GenBank
# features and is what made the downstream guard so misleading.
$event->{success}  = $features_added ? 1 : 0;
$meta{features_called} = $features_added;
$meta{viral_taxon} = $genome_in->{viral_family} if defined $genome_in->{viral_family};

# mapping<string, string>: stringify, so a consumer never has to care whether a count
# arrived as a number or a string, and JSON round-trips identically either way.
$event->{metadata} = { map { $_ => "$meta{$_}" } grep { defined $meta{$_} } keys %meta };

$genome_in->destroy_to_file($opt->output);

# Default off: a caller that runs this and the downstream stages as one shell pipeline
# fails the whole job on any non-zero stage, and rc 11 -- annotated nothing -- is common
# and legitimate.  With the stages run separately the status is what tells the caller to
# stop after this one, which is the point.
exit($opt->propagate_exit_status ? $rc : 0);

#
# Read the <prefix>.classification sidecar written by annotate_by_viral_pssm.pl -margin.
# Returns undef when the file is absent, i.e. when --margin was not used -- that is the
# normal case and must not warn.  Format is one tab-separated key/value per line, plus
# repeated "score" lines holding the full per-taxon ranking, which we do not need here.
#
sub read_classification
{
	my($file) = @_;
	return undef unless -s $file;

	open(my $fh, "<", $file) or do { warn "Could not read $file: $!\n"; return undef; };
	my %cls;
	while (<$fh>)
	{
		chomp;
		my($key, @val) = split /\t/;
		next if $key eq "score";
		$cls{$key} = $val[0];
	}
	close $fh;

	# "inf" means nothing else was in contention; leave margin unset rather than
	# putting a non-numeric value into the GTO.
	delete $cls{margin} if defined $cls{margin} && $cls{margin} !~ /^[\d.]+$/;
	$cls{margin} += 0            if defined $cls{margin};
	$cls{runner_up_bit} += 0     if defined $cls{runner_up_bit};
	$cls{below_threshold} += 0   if defined $cls{below_threshold};
	delete $cls{runner_up} if defined $cls{runner_up} && $cls{runner_up} eq "-";
	return \%cls;
}

#
# --fallback-viral-taxon support.
#

#
# Did the base script reject the classification, as opposed to falling over?
#
# The distinction is the whole gate: retrying a BLAST crash with a declared taxon would
# crash again and then attach a provenance record claiming the taxon came from the lineage
# on purpose.  EXIT_NO_REFERENCE answers this directly.
#
# The stderr match is a compatibility shim, and only for rc == 1.  Before the exit codes the
# base script exited 1 for both a rejection and a crash, so text was the only signal -- and
# the P3 environment ships a deployed annotate_by_viral_pssm.pl that can be older than this
# wrapper.  These two strings are the only two places that script prints a cutoff rejection
# (unrestricted, and under -vfam); keep them in step if the wording ever changes.  Once no
# pre-exit-code base script is reachable this can drop to the rc test alone.
#
sub classification_rejected
{
	my($rc, $file) = @_;
	return 1 if $rc == EXIT_NO_REFERENCE;
	return 0 unless $rc == EXIT_USAGE;
	return 0 unless -s $file;
	open(my $fh, "<", $file) or return 0;
	my $hit = 0;
	while (<$fh>)
	{
		if (/No matching reference contigs/ || /No reference contig of family/)
		{
			$hit = 1;
			last;
		}
	}
	close $fh;
	return $hit;
}

#
# The 34 annotation taxa, asked of the base script rather than re-derived here.
#
# -list-vtax reads the same $pdir/$json this run is using, so a data directory with a
# different taxon set answers for itself.  Returns the empty list on any failure; the
# caller then simply does not fall back, which is the same outcome as an unmatchable
# lineage and is the safe direction to fail in.
#
sub valid_viral_taxa
{
	my($pdir, $json) = @_;
	my $out = "";
	my $ok = run(["annotate_by_viral_pssm.pl", "-list-vtax", "-pssm", $pdir, "-j", $json],
		     ">", \$out, "2>", \my $err);
	if (!$ok)
	{
		warn "--fallback-viral-taxon: could not list the valid taxa (rc=" . ($? >> 8) . "); not falling back.\n$err";
		return ();
	}
	return grep { length } map { my $t = $_; $t =~ s/^\s+|\s+$//g; $t } split /\n/, $out;
}

#
# Walk the GTO's own lineage and return the most specific element that is an annotation
# taxon, as ($taxon, $lineage_element_it_came_from).
#
# Most-specific-first, because the lineage runs general -> specific and the deepest match is
# the most informative: an Orthopneumovirus muris record must land on the species taxon and
# not stop at the genus one rank above it.
#
# Matching is exact string, never prefix or substring.  "Orthopneumovirus" is a strict
# prefix of "Orthopneumovirus_muris", so a prefix match would silently mis-assign one of
# them.  The only normalization is space -> underscore, which is how that one binomial taxon
# is spelled on disk.
#
# A Bunyavirales lineage ("... ; Hantaviridae; Orthohantavirus; Orthohantavirus xyz") simply
# finds no match at genus or species and lands on the family, which is the rank LowVan
# annotates that order at.  A lineage that stops above any annotation taxon returns nothing
# and the caller reports it -- that is the genus-less population, for which --viral-family
# is the tool, not this one.
#
sub taxon_from_lineage
{
	my($gto, $valid) = @_;
	my %valid = map { $_ => 1 } @$valid;
	return (undef, undef) unless %valid;

	my @lineage;
	my $tax = $gto->{taxonomy};
	if (ref($tax) eq 'ARRAY')          { @lineage = @$tax }
	elsif (defined $tax && $tax ne "") { @lineage = split /;/, $tax }
	elsif (ref($gto->{ncbi_lineage}) eq 'ARRAY')
	{
		# The structured form some GTOs carry instead: [[name, taxid, rank], ...].
		@lineage = map { ref($_) eq 'ARRAY' ? $_->[0] : $_ } @{$gto->{ncbi_lineage}};
	}

	for my $elt (reverse @lineage)
	{
		next unless defined $elt;
		(my $name = $elt) =~ s/^\s+|\s+$//g;
		next if $name eq "";
		(my $t = $name) =~ s/\s+/_/g;
		return ($t, $name) if $valid{$t};
	}
	return (undef, undef);
}




