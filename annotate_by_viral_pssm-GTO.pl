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

# --- fallback: BLASTn could not place the genome, so declare the taxon ourselves ------
#
# Order matters and is the point of the option: BLASTn gets first refusal, so any genome
# it can place is annotated exactly as it is without the flag.  The lineage is consulted
# only for genomes that would otherwise produce no annotation at all.
#
# The trigger is the stderr text, not the exit code.  The wrapper's own note above
# explains why: a BLAST crash and the base script's graceful "no reference match" exit(1)
# are indistinguishable by rc, and retrying a crash with a declared taxon would just crash
# again with a misleading provenance record attached.
my ($fallback_taxon, $fallback_rank);
if ($opt->fallback_viral_taxon && !$ok && classification_rejected("$here/$prefix.stderr.txt"))
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
    }
}

if (!$ok)
{
    print STDERR "Viral Annotation run failed with rc=$?. Stderr:\n";
    copy("$here/$prefix.stderr.txt", \*STDERR);   # was hardcoded Viral_Anno.stderr.txt + mislabeled "Stdout" (gist #17)
    # NOTE: execution continues below and still writes an output GTO. See caveat in the porting notes --
    # a genuine BLAST crash and the base script's graceful "no reference match" exit(1) are
    # indistinguishable by return code, so we do not hard-abort here.
    # (Exception: under --viral-taxon the base script warns instead of exiting on a low -mcb,
    # so a non-zero rc there is unambiguously a real failure.  --viral-family does NOT get
    # that exception -- it still has to choose a genus, so a low -mcb exits 1 there as it
    # does with no declaration at all.)
}


    
my $event = {
    tool_name => "LowVan Annotate",
    tool_version => $tool_version,
    execution_time => scalar gettimeofday,
    parameters => \@params,
};

my $event_id = $genome_in->add_analysis_event($event);


## Parse the generated peptide file. We collect the CDS, mature_peptides, and RNAs then
## add features so that we can register the counts.
##

my $viral_family;

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
		}
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
	}

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
}

$genome_in->destroy_to_file($opt->output);

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
# We match on the text rather than the exit code deliberately.  The base script exits 1
# both when no reference clears -mcb and when BLAST itself dies, and retrying a crash with
# a declared taxon would just crash again while attaching a provenance record claiming the
# taxon came from the lineage on purpose.  These two strings are the only two places the
# base script prints a cutoff rejection (unrestricted, and under -vfam); keep them in step
# with annotate_by_viral_pssm.pl if that wording ever changes.
#
sub classification_rejected
{
	my($file) = @_;
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
		warn "--fallback-viral-taxon: could not list the valid taxa (rc=$?); not falling back.\n$err";
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




