# Viral Annotation
This repo contains code and data for improving viral annotation.  It currently covers the  *Bunyavirales*, *Coronaviridae*, *Filoviridae*, *Orthomyxoviridae*, *Paramyxoviridae*, and *Pneumoviridae*. The overall goal is to create a low-tech solution for consistently and accurately annotating viral proteins across entire viral families and to cover cases where we do not have bespoke species-specific annotations from VIGOR4.<br>

This program is not intended to be used as a *de novo* protein or ORF discovery tool.  It is designed to find proteins that we already know to exist.  


## Covered viral taxa

### Bunyavirales:
Arenaviridae<br>
Fimoviridae<br>
Hantaviridae<br>
Nairoviridae<br>
Peribunyaviridae<br>
Phasmaviridae<br>
Phenuiviridae<br>
Tospoviridae<br>

### Coronaviridae
Alphacoronavirus<br>
Betacornoavirus<br>
Gammacoronavirus<br>
Deltacoronavirus<br>

### Filoviridae:
Orthoebolavirus<br>
Orthomarburgvirus<br>

### Orthomyxoviridae:
Influenza A virus<br>
Influenza B virus<br>
Influenza C virus<br>
Influenza D virus<br>
Isavirus<br>

### Paramyxoviridae:
Aquaparamyxovirus<br>
Ferlavirus<br>
Henipavirus<br>
Jeilongvirus<br>
Metaavulavirus<br>
Morbillivirus<br>
Narmovirus<br>
Orthoavulavirus<br>
Orthorubulavirus<br>
Paraavulavirus<br>
Pararubulavirus<br>
Respirovirus<br>

### Pneumoviridae:
Orthopneumovirus<br>
Metapneumovirus<br>
Orthopneumovirus muris<br>

## Dependencies

Unless otherwise stated, the programs described in this repo are written and tested in in perl (v5.38.0).

The script(s) have the following dependencies:<br>

External CPAN perl modules:<br>

Data::Dumper<br>
File::(including Copy, Path, SearchPath, Temp, Slurp)<br>
Getopt(Long and Descriptive)<br>
IPC::Run<br>
JSON::XS<br>
Time::HiRes<br>
Cwd<br>

It also uses `gjoseqlib.pm` which is perl module that was written by Gary Olsen at the University of Illinois.  This is used for basic sequence manipulation.  You can get the latest version of this module by downloading it from Gary's repo here: https://github.com/TheSEED/seed_gjo/  <br>

There are 2 BV-BRC supported modules as well:<br>
GenomeTypeObject.pm (https://github.com/BV-BRC/p3_core/blob/master/lib/GenomeTypeObject.pm).  This perl module contains all the necessary tooling to read and write to and from a Genome Type Object, which is a JSON representation of a genome and all analysis events.<br>

The P3DataAPI "https://www.bv-brc.org/docs/cli\_tutorial/command\_list/P3DataAPI.html) is used to generate unique ids for features in the BV-BRC annotation system.<br>


The program(s) run the blast suite of tools from NCBI.  The current version requires:<br>
`blastn: 2.16.0+`<br>
`tblastn: 2.16.0+`<br>

It is not guaranteed to work on other versions of BLAST.  It uses the JSON output of BLAST and other versions have slightly different JSON structures.  <br>

For internal users, source:<br>
`/vol/patric3/cli/ubuntu-cli/user-env.sh`<br>

## Repo Contents
`annotate_by_viral_pssm.pl` the perl script that runs the BLASTs and calls the ordinary proteins (CDSs), mature peptides, and location-based features such as RNAs. <br><br>

`annotate_by_viral_pssm-GTO.pl` this perl script runs annotate_by_viral_pssm.pl and creates a GTO as output. Note that its help options are slightly different.<br>
It is run in the following way: `annotate_by_viral_pssm-GTO.pl  -x [file_prefix] -i Input.gto -o Output.gto` with other options in the help menu.<br><br>

`get_transcript_edited_features.pl` This script reads a GTO (that has already been processed by annotate_by_viral_pssm-GTO.pl) and finds sequences that have undergone transcript editing, updating the resulting CDS or mat_peptide to ensure we get the correct protein sequence.<br><br>

`get_splice_variant_features.pl` This script reads a GTO (that has already been processed by annotate_by_viral_pssm-GTO.pl) and finds sequences that are the result of splicing, updating the resulting CDS to ensure we get the correct protein sequence.<br><br>

`viral_genome_quality.pl`  This script reads the GTO and evaluates the genome quality based on CDSs and mat_peptide features present, and their copy number.  It also evaluates the contigs based on copy number and the proteins they encode.  It is intended to be run downstream of annotate_by_viral_pssm-GTO.pl and get_transcript_edited_features.pl.<br><br>

`Viral_PSSM.json`  This file contains BLAST and ORF calling parameters per feature.<br><br>

`Viral-Rep-Contigs` This is the directory of representative contigs that guides the program to the closest set of PSSMs.<br><br>

`Viral-PSSMs` This is the directory of hand curated PSSMS per taxon. <br><br>

`Transcript-Editing` This directory contains fasta files of hand-curated transcripts (post editing).<br><br>

`Splice-Variants` This directory contains fasta hand-curated fasta files with splice site locations.<br><br>

`PSSM-Alignments`  This directory is not used by any program, but it contains the alignments that correspond to each PSSM.<br><br>

`Other-Scripts` is a directory other non-essential but useful scipts and files related to the development and management of these tools.  It currently contains several readme files s and the tools for building pssms.<br><br>

## How to run annotate_by_viral_pssm.pl
`annotate_by_viral_pssm.pl [options] -i subject_contig(s).fasta`<br><br>

Options include:
```

		-h help
		-i input subject contigs in fasta format
		-t declare a temp file (d = random)
		-tax declare a taxonomy id (D = 11158 )
		-vtax declare the annotation taxon, instead of detecting it by BLASTn (see below)
		-list-vtax print the valid -vtax values and exit
		-vfam declare only the family, and let BLASTn pick the genus within it (see below)
		-list-vfam print the valid -vfam values and exit
		-skip-classification with -vtax, only BLASTn that taxon own reference contigs (bulk reruns)
		-margin minimum winner/runner-up bit score ratio, e.g. -margin 1.3 (off by default; warns
		   and records, never rejects -- see below)
		-g Genome name (D = Paramyxoviridae);
		-k Keep internal stop codons (D = off) if you think that your genome will have stops
		   within the PSSM, but still want to make a call over that region creating a pseudo gene.

		-min minimum contig	length (d = 1000) # otherwise the genome is rejected
		-max maximum contig length (d = 25000) # for reference Measles is 15894 and beilong is 19,212

        -opt Options file in JSON format which carries data for match (D = /home/jjdavis/Viral_PSSM.json)
		-l Representative contigs directory (D = /home/jjdavis/bin/Viral-Rep-Contigs)
		-p Base directory of PSSMs   (D = /home/jjdavis/bin/Viral-PSSMs)
	      Note that this is set up as a directory of pssms
	      right now this is hardcoded as: "virus".pssms within this directory.
```
Hard-coded locations currently exist as the defaults for -opt, -l, and -p.  Since that is annoying, you 
might want to run something like:<br>
`perl -i -pe 's/\/home\/jjdavis\/bin/the path to your bin/g' annotate_by_viral_pssm.pl`, or you could edit lines 70-72 by hand (but note that these are the line numbers at the time I wrote this).<br><br>




There is also a set of debugging parameters that I use frequently:
```
 -tmp keep temp dir
 -no no output files generated, to be used in conjunction with one of the following:
  -dna print only genes to STDOUT 
  -aa print proteins to STDOUT
  -tbl print only feature table to STDOUT
  -ctbl [file name] concatenate table results to a file (for use with many genomes)
```


## Declaring the taxon yourself with -vtax

By default the taxon is chosen by BLASTn against `Viral-Rep-Contigs/` (see Step 1 below).  If you
already know what the genome is -- or it is divergent enough that no reference clears `-mcb`, or
BLASTn keeps landing on the wrong neighbouring genus -- you can declare it:

```
annotate_by_viral_pssm.pl -i contigs.fasta -p out -vtax Phenuiviridae
annotate_by_viral_pssm-GTO.pl -i in.gto -o out.gto --viral-taxon Henipavirus
```

The classification BLASTn still runs, and prints a warning to stderr if it disagrees with you,
naming both its own pick and the best reference *within* the taxon you declared.  The `-mcb`
rejection becomes a warning rather than a fatal error, since you have taken that call.  The
closest-genome columns report the best reference within your declared taxon, so the output table
and the GTO `close_genomes` record stay internally consistent.

### Speeding up bulk reruns

By default `-vtax` still runs the full classification BLASTn over every reference file, so it can
tell you when BLASTn disagrees with your declared taxon.  That sweep is flat in genome size and
linear in the number of references: ~0.15 s per reference, so ~11.5 s at the 78 references in place
when the table below was measured and ~13 s at the current 87.  For a one-off that does not matter;
for a
100k-genome rerun it dominates.

`-skip-classification` restricts the sweep to the reference contigs of the declared taxon only
(1–14 files, median 1):

```
annotate_by_viral_pssm.pl -i contigs.fasta -p out -vtax Morbillivirus -skip-classification
```

Output is byte-identical to the same run without the flag.  What you give up is the cross-check
warning: nothing will tell you that BLASTn would have picked a different taxon.  Use it only when
the taxon is already known to be right.

Measured end-to-end, small fragment inputs, at 78 references:

| Declared taxon | PSSMs | Full sweep | Restricted | Speedup |
|---|---|---|---|---|
| Morbillivirus | 26 | 15.9 s | 4.5 s | 3.5x |
| Orthopneumovirus | 66 | 22.2 s | 10.4 s | 2.1x |
| Alphainfluenzavirus | 125 | 31.3 s | 19.6 s | 1.6x |
| Betacoronavirus | 410 | 74.5 s | 63.0 s | 1.2x |

The saving is the same ~11.5 s in every case; the ratio varies because the remaining cost is the
per-PSSM tBLASTn stage, which this flag does not touch and which scales with how many PSSMs the
declared taxon has (18 to 410, median 26).

### What to type


Run `annotate_by_viral_pssm.pl -list-vtax` for the authoritative list.  There are 34 values, and
they are the same strings used for the `Viral-PSSMs/<name>.pssms/` directories, the top-level keys
of `Viral_PSSM.json`, and the `Viral-Rep-Contigs/<name>[.N].dna` prefixes.

**They are not all families.**  Each clade is split at whatever rank its PSSMs discriminate:

| Rank | N | Values |
|---|---|---|
| Family | 8 | `Arenaviridae`, `Fimoviridae`, `Hantaviridae`, `Nairoviridae`, `Peribunyaviridae`, `Phasmaviridae`, `Phenuiviridae`, `Tospoviridae` -- every Bunyavirales clade |
| Genus | 25 | `Alpha`/`Beta`/`Delta`/`Gammacoronavirus`; `Alpha`/`Beta`/`Delta`/`Gammainfluenzavirus`, `Isavirus`; `Orthoebolavirus`, `Orthomarburgvirus`; `Aquaparamyxovirus`, `Ferlavirus`, `Henipavirus`, `Jeilongvirus`, `Metaavulavirus`, `Morbillivirus`, `Narmovirus`, `Orthoavulavirus`, `Orthorubulavirus`, `Paraavulavirus`, `Pararubulavirus`, `Respirovirus`; `Metapneumovirus`, `Orthopneumovirus` |
| Species | 1 | `Orthopneumovirus_muris` |

Things that trip people up:

- There is no `Coronaviridae`, `Paramyxoviridae`, `Filoviridae`, or `Pneumoviridae` entry -- name
  the genus.  Conversely there is no `Orthohantavirus`; the Bunyavirales resolve only to family, so
  it is `Hantaviridae`.
- Names are current ICTV.  Influenza A is `Alphainfluenzavirus` (not `InfluenzaA`, `FluA`, or a
  subtype like `H1N1`), influenza D is `Deltainfluenzavirus`, and Ebola is `Orthoebolavirus` (the
  genus was renamed from `Ebolavirus`).
- `Orthopneumovirus_muris` is the only name with an underscore; it stands in for the space in the
  binomial *Orthopneumovirus muris*.  Because `Orthopneumovirus` is a strict prefix of it, matching
  is exact-string only -- a case-only mismatch is accepted with a note, nothing else is.
- A segmented genome takes one taxon for all its segments; segment identity comes from the
  per-feature `segments` block in `Viral_PSSM.json`.
- The value becomes `viral_family` in the GTO, which is what the downstream steps use to pick
  `Transcript-Editing/<name>/` and `Splice-Variants/<name>/`.  Only 17 of the 34 taxa have a
  transcript-editing directory and only 5 have a splice-variant one, so declaring the taxon can
  turn those steps on or off.

## Declaring only the family with -vfam

`-vtax` needs the genus.  A great many GenBank records do not have one: "Bat coronavirus HKU10",
"Paramyxoviridae sp.", anything submitted before its genus existed.  In the 5,000-genome failure-set
rerun, 322 records (7.2%) had a lineage that stopped above genus, so `-vtax` could not be used on
them at all.

`-vfam` takes the rank they do have.  It restricts the classification BLASTn to the reference
contigs of the taxa in that family and lets BLASTn choose among them:

```
annotate_by_viral_pssm.pl -i contigs.fasta -p out -vfam Coronaviridae -mcb 50
annotate_by_viral_pssm-GTO.pl -i in.gto -o out.gto --viral-family Coronaviridae --min-contig-bit 50
```

```
Restricting BLASTn sweep to 10 of 87 reference contigs (-vfam Coronaviridae: Alphacoronavirus,
Betacoronavirus, Deltacoronavirus, Gammacoronavirus)
Within-family margin 1.11 (runner-up Gammacoronavirus, bit = 133.176)
Annotating as Betacoronavirus   Bit = 148.252   (genus chosen within declared family Coronaviridae)
```

**It is a different kind of flag from `-vtax`, and the difference is the whole point.**  `-vtax`
declares an answer, so it never rejects.  `-vfam` still has to *make* a call, so `-mcb` applies
normally and a genome that clears nothing is rejected exactly as it would be with no declaration.

That also means restriction alone rescues nothing: throwing references away can only lower the best
bit score, so a genome that failed `-mcb 150` unrestricted still fails it under `-vfam`.  What
`-vfam` buys is *trustworthiness in the band below the cutoff*, which is where these records
actually sit -- the genus-less failures have a median best hit of 114 bits.  Pair it with a lowered
`-mcb`:

| | unrestricted | `-vfam` |
|---|---|---|
| accuracy of a call in the 50-150 bit band | ~60% | ~93% |
| what a wrong call costs | an unrelated family's PSSMs | a sibling genus' PSSMs |

An unrestricted near-miss can land anywhere; a family-restricted one is wrong only between close
relatives, whose PSSM sets overlap.  So `-vfam -mcb 50` is a reasonable second pass over records
that failed classification and have no genus, in a way that a bare `-mcb 50` is not.

`-margin` is also more informative here, because the runner-up is a sibling genus rather than an
unrelated family.  With `-vfam` the within-family margin is printed whether or not you pass
`-margin`, and the `<prefix>.classification` sidecar is written either way, with three extra lines:

```
scope           family
family          Coronaviridae
min_contig_bit  50
```

The GTO wrapper puts `classification_scope` and `classification_family` on the `close_genomes`
record, and -- importantly -- falls back to the sidecar's `annotated_as` when setting
`viral_family`.  A `-vfam` run that finds no features still emits a GTO with `viral_family` set, so
the downstream transcript-editing, splice-variant and quality steps run to completion instead of
dying on `GTO has no viral_family field`.

Other behaviour:

- Valid values come from a `family` field on each taxon in `Viral_PSSM.json`; run
  `annotate_by_viral_pssm.pl -list-vfam` for the list.  There are **13**: `Coronaviridae`,
  `Filoviridae`, `Orthomyxoviridae`, `Paramyxoviridae`, `Pneumoviridae`, and the eight Bunyavirales
  families, which are annotation taxa in their own right.
- For those eight, `-vfam` degenerates to a one-taxon scope -- the same references `-vtax` would
  search, but with `-mcb` still enforced.
- **`-vtax` accepts a family name** and routes it here, with a note on stderr.  Exact taxon match is
  tried first, so `-vtax Hantaviridae` keeps its existing meaning.
- Cannot be combined with `-vtax` or `-skip-classification`.

## Falling back to the GTO's own lineage with --fallback-viral-taxon

`-vtax` and `-vfam` both make you choose up front, for every genome, whether to trust BLASTn or to
override it.  `--fallback-viral-taxon` (GTO wrapper only) does not choose: it runs the normal
classification, and only if that finds no reference above `-mcb` does it retry once with `-vtax`
taken from the lineage already in the input GTO.

```
annotate_by_viral_pssm-GTO.pl -i in.gto -o out.gto --fallback-viral-taxon
```

```
No matching reference contigs with bit score greater than 150
--fallback-viral-taxon: classification found no reference above the cutoff; retrying with
-vtax Alphainfluenzavirus (from GTO lineage element 'Alphainfluenzavirus').
```

BLASTn gets first refusal, so a genome it can place is annotated exactly as it is without the flag
-- this is not an override.  The lineage is consulted only for genomes that would otherwise produce
no annotation at all, which in the failure-set rerun was the single largest recoverable bucket.

Details:

- **The lineage comes from the GTO's own `taxonomy` field**, not from any GenBank record it was
  built from.  Note that `rast-create-genome --from-genbank` does *not* populate `taxonomy` -- it
  writes only `ncbi_taxonomy_id` and `scientific_name` -- so a GTO built that way has nothing to
  fall back to and the flag will say so and do nothing.  (`ncbi_lineage`, if present, is read as an
  alternative.)
- The lineage is walked **most specific first**, so an `Orthopneumovirus muris` record lands on the
  species taxon rather than stopping at the genus one rank above it.  Matching is exact-string with
  spaces mapped to underscores; a Bunyavirales lineage finds nothing at species or genus and
  correctly lands on the family.
- A lineage that resolves no deeper than a rank LowVan has no taxon for -- the genus-less 7.2% --
  prints a note naming the lineage and does not retry.  `--viral-family` is the tool for those.
- **The trigger is the rejection message, not the exit code.**  The base script exits 1 both when it
  rejects the classification and when BLAST itself dies, and retrying a crash with a declared taxon
  would just crash again while claiming the taxon was chosen deliberately.
- Cannot be combined with `--viral-taxon`: the taxon is already declared, and under `--viral-taxon`
  a low `-mcb` is a warning rather than a failure, so nothing would ever trigger.  It *can* be
  combined with `--viral-family`, which can still fail `-mcb`; the retry drops `-vfam`, since a
  genus from the lineage is strictly more specific than a family scope.
- On the retry the GTO records where the taxon came from: `classification_scope`
  (`lineage_fallback`), `classification_taxon` and `classification_lineage` on the `close_genomes`
  record, `-vtax` in the analysis event's parameters, and `viral_family` set from the declared taxon
  so the downstream steps run even if no features were called.  The first attempt's stderr is kept
  as `<prefix>.stderr-classification.txt` rather than being overwritten by the retry.

## Checking how confident a classification is, with -margin

`-mcb` asks whether the best BLASTn hit is *strong enough*.  It does not ask whether it beat the
next taxon, and those are different questions: a fragment can score 183 bits -- comfortably over
the default cutoff of 150 -- while the runner-up taxon scores 145.  That is a coin-flip being
recorded as an annotation.

`-margin` scores the second question.  It takes the ratio of the winning taxon's bit score to the
runner-up's and warns when the ratio is below the value you give:

```
annotate_by_viral_pssm.pl -i contigs.fasta -p out -margin 1.3
annotate_by_viral_pssm-GTO.pl -i in.gto -o out.gto --margin 1.3
```

```
WARNING: low classification margin 1.19 (< -margin 1.3): Jeilongvirus scored 217.285 but
Morbillivirus scored 183.165. The winning taxon is only 19% ahead of the runner-up; treat this
call as unconfirmed.
```

Notes on the behaviour:

- **It is off by default.**  Without the flag nothing is computed, warned about, or written, and a
  run is bit-for-bit what it was before the option existed.
- **It never rejects.**  The genome is annotated either way; this marks the call, it does not gate
  it.  There is no default threshold, because the only value we have evidence for was fitted to a
  single case that has since been fixed by adding a reference.
- **It costs nothing.**  The per-taxon best score is already known from the classification sweep.
- **It is refused with `-skip-classification`**, which searches one taxon's references and so has
  no runner-up to compare against, and with any value `<= 1`, which would accept anything.

Whether or not it warns, `-margin` writes a `<prefix>.classification` sidecar next to the other
outputs -- tab-separated `key<TAB>value` lines, then the full per-taxon ranking:

```
taxon           Jeilongvirus
reference       Jeilongvirus.2.dna
bit             217.285
runner_up       Morbillivirus
runner_up_bit   183.165
margin          1.19
margin_threshold 1.3
below_threshold 1
annotated_as    Jeilongvirus
score   Jeilongvirus    217.285 Jeilongvirus.2.dna
score   Morbillivirus   183.165 Morbillivirus.2.dna
...
```

`margin` is the literal string `inf` when no other taxon scored above zero.  The GTO wrapper reads
the sidecar and adds `margin`, `runner_up`, `runner_up_value` and `margin_below_threshold` to the
`close_genomes` record, so the confidence of the call travels with the genome.

## Step 1.  Calling features based on PSSMs

The code is currently designed to work on the *Paramyxoviridae*, *Bunyavirales*, *Filoviridae*, and *Pneumoviridae*, although more taxa are planned.  As depicted in the image below, it first performs a BLASTn against a small set of representative genomes for each genus.  Then it sorts the results by bit score and chooses the best match.<br>

For each genus, there is a directory of PSSMs corresponding to each known protein for that genus. The PSSMs are derived from a set of hand curated alignments. In the next step, it cycles through each directory of PSSMs (there may be more than one PSSM per protein), choosing the best tBLASTn match per pssm. <br>

Note that it assumes your genome will have the same set of proteins as the nearest match. This is why it is not intended to be used as a discovery tool.  In the event that a new protein is found, a new PSSM must be added to the PSSM directory.  <br><br>

![Anno-Strategy](https://github.com/jimdavis1/Viral_Annotation/assets/7661533/0d6a3a44-47af-40bf-852d-5ddda250ad94)

<br><br>Finally it performs any special rules on the proteins/ORFs.  These rules are currently encoded in a JSON file called `Viral_PSSM.json`. The following is a description of the current JSON strucutre.<br><br>

```
  "Arenaviridae": {
    "family": "Arenaviridae",
    "segments": {
      "Small RNA Segment": {
        "max_len": 3741,
        "min_len": 3061,
        "replicon_geometry": "linear"
      },
      "Large RNA Segment": {
        "max_len": 8014,
        "min_len": 6556,
        "replicon_geometry": "linear"
      }
    "features": {
      "GPC": {
        "anno": "Pre-glycoprotein polyprotein GP complex (GPC protein)",
        "bit_cutoff": 100,
        "copy_num": 1,
        "coverage_cutoff": 0.65,
        "downstream_ext": 1,
        "feature_type": "CDS",
        "max_len": 617,
        "min_len": 455,
        "non_pssm_partner": ["Small Segment Stemloop"],
        "segment": "Small RNA Segment",
        "upstream_ext": 1
      },
      ...
```

The Viral_PSSM.json file is in a regular state of development, so this may change slightly, but the above shows an example for, *Arenaviridae*, and a single protein, GPC.  The three highest level keys are `family`, the ICTV family the annotation taxon belongs to (which is what `-vfam` groups on -- every taxon must have one, and for the Bunyavirales, which are themselves families, it repeats the taxon name); `segments`, which contains information on segments that are used for genome quality evaluation; and `features`, which currently contains information on CDS, mat_peptide, and RNA features.  A fourth, `close_genomes`, maps each `Viral-Rep-Contigs/` filename to the BV-BRC genome it came from. <br>

The following is a non-exhaustive description of fields that are used in the JSON<br><br>

`max_len and min_len` maximum or minimum length of a contig or feature that is expected and evaluated by the genome quality checker (not all features or contigs will have this).  The boundaries are currently very crude, and not bound by any sort of statistics, but effective<br><br>

`replicon_geometry` this is not currently used, but carries info on the genometry of the replicon and inserted into the GTO by the quality tool<br><br>

`copy_num` expected copy number of a feature<br><br>

`coverage_cutoff` blast subject coverage for calling a feature<br><br>

`upstream_ext and downstream_ext` tells the program if it can look upstream for a Met start or downstream for a stop codon<br><br>

`feature_type` currently CDS, mat_peptide, or RNA <br><br>

`segment` which segment a feature belongs to (used by the quality tool)<br><br>

`non_pssm_partner` used for placing a location based feature<br><br>

There are other fields that are not depicted in the example, including:<br>
`PMID` which contains the PubMed ID for one or more DLITS.  A DLIT is an example of an important paper that either defines the function or sequence of a feature. <br><br>

`"special": "transcript_edit"`  This field tells the program that an external program is being used to make a call.  In this case, `transcript_edit` is used to denote a feature that undergoes transcript editing and is found by using `get_transcript_edited_features.pl`.  `splice` is also a valid field and triggers the search for splice variant features.<br><br>


## Get Transcript Edited Features
Transcript editing is a phenomenon that occurs in the phosphoproteins of the *Paramyxoviridae* and the glycoproteins of the *Filoviridae*.  It occurs when the RNA-Dependent RNA polymerase encounters a region of low complexity and pauses.  The pause allows for the insertion of one or more new nucleotides into the transcript, which causes a frame shift. Thus, the amino acid sequence is not a direct translation of what is encoded in the genome.  We solve this problem by hand-curating a set of transcripts in their post-editing state. These are found in the `Transcript-Editing` directory.  We then BLAST these against the the contig, and for BLASTn matches with high enough scores, the alignment gap is filled in using the nucleotide sequence of closest curated transcript.  Currently in order to do this, the following strict BLASTn criteria must be met: <br>

1.  The match must have >= 95% nucleotide identity
2.  The match must have >= 95% query coverage
3.  The match must have <= 2 gap characters in the subject
4.  The gap characters must occur consecutively in a run <br><br>

These parameters are controlled using `--id`, `--cov`, and `--gaps` options, respectively.  The requirement for consecutive gap characters is hard-coded.
<br><br>

Because we may encounter a decent BLASTn match, but not have sufficient %identity, %query coverage, or there may be additional naturally-occurring gaps in the subject, this program will call a feature covering the BLASTn match when the above inclusion criteria are not met.  However, it will not attempt to correct the subject sequence.  Instead it will call a `partial_cds` feature and will not attempt a translation.  Parameters setting the minimum BLAST requirements for this type of feature call are `--eval`, `--lower_pid`, and `--lower_pcov`, which set the maximum BLAST e-value, and the minimum percent identity, and query coverage for consideration. <br><br>

Full usage for this program is as follows:


```	--input STR (or -i)    Input GTO
	--output STR (or -o)   Output GTO
	--cov INT (or -c)      Minimum BLASTn percent query coverage (D = 95)
	--id INT (or -p)       Minimum BLASTn percent identity  (D = 95)
	--gaps INT (or -g)     Maximum number of allowable gaps (D = 2)
	--e_val NUM (or -e)    Maximum BLASTn evalue for considering any HSP
	                       (D = 0.5)
	--lower_pid            Lower percent identity threshold for a feature
	                       call without transcript editing correction (D
	                       = 80)
	                       aka --lpi
	--lower_pcov           Lower percent query coverage for for a feature
	                       call without transcritp editing correction (D
	                       = 80)
	                       aka --lpi
	--threads INT (or -a)  Threads for the BLASTN (D = 24))
	--json STR (or -j)     Full path to the JSON opts file
	--dir STR (or -d)      Full path to the directory hand curated
	                       transcripts
	--tmp STR (or -t)      Declare name for temp dir (D = randomly named
	                       in cwd)
	--help (or -h)         Show this help message
	--debug (or -b)        Enable debugging
	
```
<br><br>


## Get Splice Variant Features
The use of splicing is fairly common in viruses, and is currently necessary for calling many proteins in *Influenza*.

In order to support splice variant calling we maintain a directory of hand-curated DNA sequences with the coordinates of the splice site carefully delineated.  These are found in the `Splice-Variants` directory.  We find these these by performing a BLASTn search against the contig using our curated sequences as the query.  Then for the BLASTn matches with high enough scores, the splice is made using the curated coordinates in the fasta header.<br> 

It is worth noting that:
1.  Query sequences don't need to be be aligned, but it is easier to deal with them if they are.
2.  The coordinates are sequence specific.  
3.  If query sequences are derived from an alignment they should not contain gap characters upstream of the splice junction. Gaps should be removed and coordinates updated accordingly.
4.  All query sequences must be in the forward direction  
<br>

Fasta headers are formatted in the following way:<br>

`>valid_sequence_ID SD:SD_Region_Start-SD_Region_End;Last_nt_of_SD  SA:SA_Region_Start-SA_Region_End;First_nt_of_SA`
<br>
where SD is sequence donor, and SA is sequence acceptor.  Here is what one looks like:
`>1413195.5 SD:371-381;373 SA:491-504;503`
 <br>
The fasta header is parsed based on a regex, so colons and semicolons should not exist in the idenifier.<br>

Finally, the code will not call a splice variant unless a valid sequence donor and sequence acceptor site exist somewhere in the set of query sequences.  That is, they do not have to match exactly between the best query-subject match, but they do need to have been seen before. <br>


Full usage for this program is as follows:

```
	--input STR (or -i)    Input GTO
	--output STR (or -o)   Output GTO
	--cov INT (or -c)      Overall Minimum BLASTn percent query coverage
	                       (D = 95)
	--id INT (or -p)       Overall Minimum BLASTn percent identity  (D =
	                       95)
	--threads INT (or -a)  Threads for the BLASTN (D = 24))
	--json STR (or -j)     Full path to the JSON opts file
	--dir STR (or -d)      Full path to the directory hand curated
	                       transcripts
	--tmp STR (or -t)      Declare name for temp dir (D = randomly named
	                       in cwd)
	--help (or -h)         Show this help message
	--debug (or -b)        Enable debugging
```

## Genome Quality Tool
As described above, the JSON file that contains information about the features also contains information about copy number of features and contigs.  The quality tool assess the the following things:<br>

1.  The number of ambigous bases per contig
2.  The number of expected segments 
3.  The length of each segment relative to what is expected
3.  The number of expected occurrences of each non-variable feature
4.  The legnth of each non-variable feature<br><br>

Currently, the tool mostly looks for CDS features, but it looks for some mat_peptides in the Filoviridae. 

The output is two tables: one is a contig report and the other is a feature report.  If any given contig or feature causes the genome quality to be "poor" the reason for the call is provided.  

Usage statement for the tool:
```viral_genome_quality.pl [-ahijop] [long options...]
	--ambiguous NUM (or -a)  Fraction of ambiguous bases, (Default = 0.01)
	--input STR (or -i)      Input GTO
	--output STR (or -o)     Output GTO
	--prefix STR (or -p)     Genome Quality File Prefix
	--json STR (or -j)       Full path to the JSON opts file
	--help (or -h)           Show this help message
```

## General remarks on the curation and development of PSSMs and the current state of the annotations
### Paramyxoviridae

I have recently updated the way transcript-edited features are called by adding `get_transcript_edited_features.pl`.  This is up-to-date and evaluated for the glycoproteins of Ebola, and the phosphoproteins in the Paramyxos.  They were originally called by splicing two BLAST HSPs, which turned out to be problematic in a few cases. DLITs that either describe the editing site, or the subsequent amino acid sequence for the transcript-edited proteins have been added to the json.  There are a handful, like Narmovirus, where I do not think protein work has been done to prove V and W, but the predicted editing site is supported by literature. At this point, all editing sites are backstopped by literature references. <br>

### Respirovirus
Note that in the respiroviruses, there is a nomenclature discrepancy regarding whether the third phosphoprotein (+2 G) is called W or D.  Currently these are all called W by the system and I have not enountered a compelling reason (other than the historical naming) to maintain the distinction between W and D.

### Tospoviridae:
I was unable to find any acceptable publications that unambiguously define the coordinates of Gn and Gc.<br>

### Fimoviridae:
I also could not find any publications clearly showing Gn and Gc.<br>

The Fimoviridae are the most poorly characterized family that I have encountered so far.  They are  multi-segmented and variable in their smaller segments. Proteins from these segments including P5, 6, 6a, 6b, 7, and 27 are all essentially uncharacterized.  They are numbered based on appearance in the genome in which they are described, but their ordering may or may not hold up as more are sequenced.  Furthermore, the proteins that have been called P5 and P6 have little to no similarity amongst themselves (usually < 35% identity) and could all have different functions in their own right.  I chose to split these into individual sets of pssms with the annotation "Fimoviridae uncharacterized protein."  We can hang an annotation on each when we learn what it does.  It is worth noting that due to the infrequency of these proteins, there are many low-occurrence uncharacterized proteins that did not get PSSMs and are not getting called.   The "P5" protein of Raspberry leaf blotch emaravirus is a good example here (fig|1980431.35.CDS.1).<br> 

In this family the quality checker will look for Segments 1-4 only, which correspond to the individual proteins L, GPC, N, and MOV, respectively.  Their segment lengths are highly variable, so the lenght cutoffs for segments 1-4 are based on the the lower length limit of the corresponding protein, and (the longest allowable gene + 0.5 X longest allowable gene) (this is arbitrary and could  be tuned).<br>

## Phasmaviridae
These are mostly insect virueses.  The set of genomes is highly diverse with few representatives in each genus, so the pssms only represent a fraction of the true diversity.  There were a considerable number of proteins that I could not get to cluster at 50% identity. I am currently dissatisfied with this family, so as more exemplars come in, this set should eventually get recomputed. 

## Pneumoviridae
The cleaved forms of the fusion glycoprotein differ between ortho- and metapneumoviride.  All of the orthos, except murine and close relatives, have a p27 peptide that is a real protein. This necessitated the insertion of three taxon-level directories (ortho, meta, and murine orthos).  I kept the orginal all-pneumo alignment directory which has everything and has seprate subdirectories for the mature F proteins.   The pssm directories for the three taxa contain the pssms that I had originally built for all pneumos. This means that there are a few extra pssms that won't match and can be cleaned up on a rainy day.    

## Orthomyxoviridae
For spliced proteins, whenever possible, the exact regions of the splice donor and splice acceptor sites are used as they have been described in the literature.  In many cases these correspond to specific hairpins in the RNA.  However, in C and D, the splicing is proven, but an exact splice donor site is not described in the literature, so 10nts upstream and downstream is used.  This may result in overly conservative calling, but the splices are correct and backed by literature.

In FluC the CM2 protein is a mat peptide that is cleaved from a pre-peptide.  This has not been studied in FluD, so DM2 is currently called from a met start [PMID: 31261944].

I have built a module for Isavirus, but i stopped at Thogotovirus and Quaranjavirus, becuase they lack sequenced examples. 

## Alphacoronavirus
There is a lack of bench data supporting the clevage sites of ORF1ab in the alphacoronaviruses.  The Mpro clevage sites are well conserved, and less troublesome, but the PLpro cleavage sites proved difficult.  The two residues upstream of the Mpro clevage sites for NSP4-NSP16 are well conserved and are usually L,I,M,or V followed by Q, but it's usually L.  However, there is considerable variation in the PLpro clevage sites and protein lengths for NSPs 1-4.  Becuase this was so difficult, I am providing the following table of references from which I based my PSSMs.  There were a few other papers that made predictions for FCoV and SADS-Cov, which I also used as guides [PMIDs: 38624231, PMID: 33913206].  However, I should note that I struggled with the predicted coordinates published for swine strain, PCV777 [PMID: 33913206], and the HKU10 and related bat strains [PMID: 22933277], as well as the Lucheng rat alphas, for which I currently lack a reference with clevage sites.  In summary, the C-terminus of NSP4 is clean, and the N- and C- termini of NSP5-16 are also clean (excluding the ribosomal slippage site).



| Cleavage site | Enzyme | P2 residue | P2 residue | P2 residue | Reference |
|---------------|--------|----|----|----|-----------|
| | | HCoV-NL63 | HCoV-229E | TGEV | |
| nsp1-nsp2 | Plpro | A | G | R | PMID: 11431476,16911043, 16306572 |
| nsp2-nsp3 | Plpro | A |  | G | PMID: 16476987,16911043 |
| nsp3-nsp4 | Plpro |  | A* | S* | PMID: 16476987 |
| nsp4-nsp5 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp5-nsp6 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp6-nsp7 | Mpro | V | V |  | PMID: 26948040,16911043 |
| nsp7-nsp8 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp8-nsp9 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp9-nsp10 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp10-nsp11 | Mpro | I | I |  | PMID: 26948040,16911043 |
| nsp12-nsp13 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp13-nsp14 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp14-nsp15 | Mpro | L | L |  | PMID: 26948040,16911043 |
| nsp15-nsp16 | Mpro | L | L |  | PMID: 26948040 |

### ORF4

There is an open reading frame (usually called ORF4) that occurs between the spike and envelope proteins.  There is a handful of papers suggesting that this ORF accumulates mutations in culture.  In the cultured 229E strains, this has resulted in the apperance of two small separate ORFs called ORF4a and ORF4b.  I was unable to find papers that suggesting that 4a or 4b are functional.  My current system only calls ORF4 starting from the original start position.     

## Betacoronavirus

There is considerable diversity in the genome structure of the Betacornoaviruses, but very little that conflicted across the embecos, hibecos, merbecos, nobecos, and sarbecos.  Overall, there are a few issues that desrve mentioning.  First, the PSSMs that are specific to these groups co-exist within the same module, which means that, for example, the handful of nobeco-specific PSSMs get blasted aginst your SARS-CoV-2 genome, but will have no match.  Second, in order to differentiate between the MERS and SARS genomes with and without furin cleavage sites in their spike proteins, I dialed up the bit scores of the blast matches for the S1 and S2 subunits, which is currently doing well. Third, there were only a small number of hibeco genomes, so PSSMs will likely need to be developed if we acquire more genomes.  Finally, it is worth noting that the embecos have a repeat in the NSP3 (PLPro) region of their ORF1a/b polyprotein.  In a test set of approximately 1056 reasonably good quality embeco genomes, there were thirty cases where the modules failed to call a contiguous PLPro because the BLAST result broke over multiple HSPs.  At this time, I am not going to develop a special module to handle these edge cases, but this is worth noting as a placeholder for future work.  If I begin to encounter many instances of long tandem repeats in the viruses, it may become necessary to develop a bespoke repeat detection module.

## Deltacoronavirus
There is a notable lack of experimental evidence supporting the mature peptide clevage sites in the Deltacoronaviruses.  In order to make a reasonable prediction of the sites, I aligned them against representative alphas and betas and looked for the boundaries.  In the cases of the main protease cleavage sites, these were fairly straighforward and correspond to the predicted boundaries published in (PMID: 2227823).  There is a crystal structure for a section of the PLpro (NSP3) (PMID: 33327866), but it notably excludes a well defined sequence for the C- terminus. I used the N-terminal clevage site for PLpro from that paper (and also predicted in PMID: 2227823), and the predicted clevage site of the NSP3/4 boundary from (PMID: 33327866).  Both regions are highly conserved in the alignment, but dont correspond as well with non-delta strains.  There are also some diverse deltas that lack homology in these regions.  I had to make an educated guess for the cleavage sites in those.
<br><br>
Secondly, there appears to be only one protien upstream of the PLpro (NSP3) in PDCoV.  A lot of papers call this NSP2.  I am calling this NSP1-2 to avoid the numeric confusion.  
<br>
Finally, my representative NSP11s are 7 amino acids in length and are too small to be reliably picked up by BLAST.  Even the longer NSP11s in the alphas and betas are dubious, so I am at peace with not detecting these.  If a paper shows that NSP11 actually has a function, I will go back and add a location-based assignment. 

## Gammacoronavirus
I encountered a lot of small uncharacterized ORFs in the gammas.  I didn't come across any papers showing that they actually encode proteins, but rather than lose them I created PSSMs and  annotated them as, "Gammacoronavirus uncharacterized ORF."   



