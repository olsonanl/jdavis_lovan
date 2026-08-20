# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains tools for viral genome annotation using Position-Specific Scoring Matrices (PSSMs). It covers Bunyavirales, Coronaviridae, Filoviridae, Orthomyxoviridae, Paramyxoviridae, and Pneumoviridae. The system identifies known viral proteins by BLAST-matching against curated PSSMs—it is not a de novo ORF discovery tool.

## Running the Annotation Pipeline

### Basic annotation from FASTA:
```bash
annotate_by_viral_pssm.pl -i contigs.fasta -p output_prefix
```

### Full pipeline with GTO (Genome Type Object) input/output:
```bash
# Step 1: Initial annotation
annotate_by_viral_pssm-GTO.pl -x prefix -i input.gto -o step1.gto

# Step 2: Transcript editing (for Paramyxoviridae phosphoproteins, Filoviridae glycoproteins)
get_transcript_edited_features.pl -i step1.gto -o step2.gto

# Step 3: Splice variants (for Influenza)
get_splice_variant_features.pl -i step2.gto -o step3.gto

# Step 4: Quality assessment
viral_genome_quality.pl -i step3.gto -o final.gto -p quality_report
```

Note: production (`p3x-annotate-lowvan.pl` in the `bvbrc_lowvan` repo) runs steps 2 and 3 in the
opposite order — annotate → **splice** → **transcript-editing** → quality. See `UPSTREAM-ISSUES.md`.

All three post-processing steps refuse to run on a GTO with no features (`No features in GTO`) or
no `viral_family` field, and because the stages are piped, one death cascades to rc=255 across all
of them. This is the single most common production failure mode — see
`LOWVAN-FAILURE-ANALYSIS.md`.

### Key options for annotate_by_viral_pssm.pl:
- `-j` JSON options file (default: Viral_PSSM.json)
- `-c` Representative contigs directory (default: Viral-Rep-Contigs/)
- `-pssm` PSSMs directory (default: Viral-PSSMs/)
- `-min`/`-max` Contig length bounds (default: 300/35000)
- `-mcb` Minimum contig bitscore (default: 150)
- `-vtax` Declare the annotation taxon instead of detecting it by BLASTn. Must be one of
  the 34 names from `-list-vtax`. The classification BLASTn still runs and warns on stderr if
  it disagrees; the `-mcb` rejection is downgraded to a warning. Closest-genome columns then
  report the best rep *within* the declared taxon.
- `-list-vtax` Print the valid `-vtax` values to stdout and exit (works without `-i`)
- `-vfam` Declare only the *family* and let BLASTn choose the genus within it, for records whose
  lineage stops above genus. Restricts the reference sweep to that family's taxa. **Unlike `-vtax`
  it still makes a call, so `-mcb` applies normally and a below-cutoff genome is still rejected.**
  See "Family-scoped classification" below. Incompatible with `-vtax` and `-skip-classification`.
- `-list-vfam` Print the valid `-vfam` values (13 families) and exit
- `-skip-classification` With `-vtax`, BLASTn only the declared taxon's own reference contigs
  (1–14 files) instead of all 87. Saves a flat ~12.8 s/genome (measured 2026-08-19: 0.149 s per
  blastn invocation, and the sweep grows with every reference added); the cost is that no
  cross-check against the other taxa is possible, so a wrong `-vtax` goes uncaught. For bulk reruns.
- `-margin` Minimum ratio of the winning taxon's bit score to the runner-up's, e.g. `-margin 1.3`.
  Off by default — without it nothing is computed, warned about, or written. See "Margin scoring"
  below.
- `-threads` BLAST threads
- Debug: `-tmp` (keep temp), `-no` with `-aa`/`-dna`/`-tbl` for stdout output

## Architecture

### Annotation Strategy
1. **Taxon identification**: BLASTn input contigs against `Viral-Rep-Contigs/` to find closest genus (or declare it with `-vtax` / `--viral-taxon` to skip detection)
2. **Feature calling**: tBLASTn each genus-specific PSSM from `Viral-PSSMs/{genus}.pssms/` against contigs
3. **Special processing**: Apply rules from `Viral_PSSM.json` (extensions, coverage cutoffs, copy numbers)
4. **Post-processing**: Transcript editing and splice variant detection for specific taxa

### Key Data Files
- `Viral_PSSM.json`: Feature parameters (bit cutoffs, coverage, copy numbers, segment info, PMIDs)
- `Viral-Rep-Contigs/`: Representative genomes per genus for initial classification
- `Viral-PSSMs/`: Hand-curated PSSMs organized as `{Genus}.pssms/` directories
- `Transcript-Editing/`: Post-edited transcript sequences for phosphoproteins/glycoproteins
- `Splice-Variants/`: Curated splice site sequences with header format `>ID SD:start-end;nt SA:start-end;nt`

**The 34 annotation taxa.** The same 34 names index all three data sets — `Viral-PSSMs/<name>.pssms/`,
the top-level keys of `Viral_PSSM.json`, and the `Viral-Rep-Contigs/<name>[.N].dna` prefixes — and
they are what `-vtax` accepts. The rank is *not* uniform, so "family" is the wrong word for most:

- **Family (8)**, all of Bunyavirales: `Arenaviridae`, `Fimoviridae`, `Hantaviridae`, `Nairoviridae`,
  `Peribunyaviridae`, `Phasmaviridae`, `Phenuiviridae`, `Tospoviridae`
- **Genus (25)**: the four `*coronavirus`, the four `*influenzavirus` plus `Isavirus`,
  `Orthoebolavirus` / `Orthomarburgvirus`, the twelve Paramyxoviridae genera, and
  `Metapneumovirus` / `Orthopneumovirus`
- **Species (1)**: `Orthopneumovirus_muris`, the binomial with the space written as `_`

Caveats: there is no `Coronaviridae` / `Paramyxoviridae` / `Filoviridae` entry (name the genus);
there is no `Orthohantavirus` (Bunyavirales resolve only to family); names are current ICTV, so
Influenza A is `Alphainfluenzavirus` and Ebola is `Orthoebolavirus`. `Orthopneumovirus` is a strict
prefix of `Orthopneumovirus_muris`, so any lookup on these names must be exact-string — never a
prefix, glob, or substring match. Segmented genomes take one taxon for the whole genome; segment
identity comes from the per-feature `segments` block in the JSON.

Each taxon also carries a **`family` scalar** in `Viral_PSSM.json` — the ICTV family it belongs to,
which is what `-vfam` groups on. Thirteen distinct values: Coronaviridae (4 taxa), Orthomyxoviridae
(5), Filoviridae (2), Paramyxoviridae (12), Pneumoviridae (3), and the eight Bunyavirales families,
which are annotation taxa themselves and so map to their own name. The field is **required** —
`taxon_family_map` dies naming every taxon that lacks one, so a new taxon cannot be added without
placing it in a family.

### PSSM Generation (Other_Scripts/)
The `fasta-cluster-pssm-2.pl` script builds PSSMs from protein FASTA:
```bash
cat proteins.faa | fasta-cluster-pssm-2.pl -a "Annotation" -p prefix
```
Uses MMSeqs2 clustering (~80% identity), MAFFT alignment, then PSI-BLAST PSSM generation. Output: `clusters/`, `alis/`, `corrected_alis/`, `pssms/` directories.

## Adding a reference genome

The classification BLASTn can only place a genome as close as its nearest reference, so extending
`Viral-Rep-Contigs/` is the main lever on classification failures. Use the script rather than doing
it by hand — it derives the taxon from the GenBank lineage, picks the next `<Taxon>.<N>.dna`,
matches the existing header convention, looks up the BV-BRC genome id, and adds the required
`Viral_PSSM.json` entry:

```bash
add_reference_genome.pl -i reference.gb --dry-run     # always dry-run first
add_reference_genome.pl -i reference.gb
add_reference_genome.pl -i seg1.gb -i seg2.gb         # segmented genome -> one .dna file
```

Two invariants it exists to protect:

- **`Viral_PSSM.json` → `<taxon>` → `close_genomes` is keyed by literal filename.** A `.dna` file
  with no matching entry produces empty closest-genome columns and an all-undef GTO `close_genomes`
  record (the defect `Viral-Rep-Contigs/Arenaviridae.6.dna` still has).
- **The file must be `<Taxon>.dna` or `<Taxon>.<N>.dna`.** The taxon is derived with `s/\..+//`, so
  the unnumbered file never needs renaming when numbered ones are added.

**`Viral_PSSM.json` is canonical JSON::XS output as of commit `e0a64bf`,** and re-encoding it
reproduces the file byte-for-byte. Edit it by decode → modify → re-encode with
`JSON::XS->new->pretty->canonical` and the diff stays proportional to the change; anything else
reformats all 7,400+ lines. `add_reference_genome.pl` depends on this — it re-encodes, compares
byte-for-byte, and refuses to write if the formatting has drifted, printing the JSON block to paste
by hand instead. Don't relax that check; restore the formatting instead.

Before `e0a64bf` the file was 2-space-indented and in insertion order, which no canonical encoder
reproduces. **Upstream `jimdavis1/Viral_Annotation` still has that format**, so a merge across that
boundary will conflict on the whole file. Reconcile by taking one side wholesale and re-applying the
other side's semantic changes, not by resolving hunks.

**Adding a reference changes classification for every genome, not just failing ones.** Re-run the
validation set and diff the per-genome taxon calls before relying on it — see "Validating a
reference change" below.

Reference files are not uniformly well-formed. `Viral-Rep-Contigs/Narmovirus.1.dna` has **no
trailing newline**, so naive concatenation glues the next `>` header onto its last sequence line,
silently dropping one sequence and corrupting another. The pipeline is unaffected (it BLASTs each
file separately), but anything that builds a combined fasta must re-join sequence lines rather than
copy them through.

## Validating a reference change

Adding or changing a reference shifts every classification call, so the change is only as good as
the before/after diff. The harness lives in `failure-analysis/`:

- `validation-set-500.tsv` — **use this one.** 500 genomes stratified by family × length band with
  per-row weights, drawn by `draw_stratified_sample.py`; 463 have a lineage-derived truth taxon and
  37 are genus-less. The weights recover population rates from a sample that deliberately
  over-represents rare strata.
- `validation-set.tsv` — the older 187-genome set. Kept for continuity with published results, but
  **too small for reference testing**: it reported "no regressions" for a change that the 500-genome
  set showed pushing a wrong call across `-mcb 150` into production. Prefer the 500.

Truth comes from each record's GenBank lineage, so both sets are drawn **only from failure-set
records whose lineage resolves to a LowVan taxon** — they structurally exclude the ~7% that stop at
family, and their accuracy is not population accuracy. Use them to detect regressions; use the
weights, not the raw counts, to estimate recovery.
- `build_ref_sweep.py` — concatenates `Viral-Rep-Contigs/*.dna` into one taxon-tagged BLASTn query
  (`>Taxon|file.dna|index`), so a single blastn scores all taxa at once. `--exclude` builds the
  "before" side, so both sides come out of the same script and differ only by the additions.
- `sweep_one.sh` — classify one genome, emitting winner and runner-up with bit scores. Its BLASTn
  parameters mirror the classification step in `annotate_by_viral_pssm.pl`; keep them in sync or the
  experiment stops predicting production.
- `gb2fa.py` — GenBank → fasta, for populating `$FA_DIR` (default `/tmp/gen/fa`) from
  `/home/mshukla/reannotation/genomes_gb/<family>/<accession>.gb`

```bash
python3 failure-analysis/build_ref_sweep.py /tmp/after.fa
python3 failure-analysis/build_ref_sweep.py /tmp/before.fa --exclude NewTaxon.2.dna
for side in before after; do
  awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' failure-analysis/validation-set.tsv \
  | xargs -P 8 -L1 bash -c "failure-analysis/sweep_one.sh /tmp/$side.fa \"\$@\"" _ > /tmp/$side.tsv
done
```

Then compare accuracy among calls clearing `-mcb 150` and, separately, list every genome whose call
changed. Judge both: a change that flips a call *below* the cutoff never reaches production, and
lumping the two together hides real regressions.

**Watch specifically for a call that crosses `-mcb 150` into a wrong taxon** — a below-cutoff wrong
call is harmless, the same call four bits higher is a mis-annotation. That is the failure mode a
new reference introduces, and it is rare enough that a 187-genome set missed it. Check the
winner/runner-up margin on anything newly accepted: at 86 references the one wrong accepted call of
the 500-genome run had the smallest margin of all 363 (1.26 against a median of 15.4). Also re-run
the full-genome regression
(`A-California-07-2009-H1N1-6.contigs.fasta`) — the feature table should be identical **apart from
the per-run allocated genome id in column 1**, which varies between runs regardless of references.

Current state: Tier 1 (`b4a7bc6`) and Tier 2 (`9e189f4`) of
`failure-analysis/proposed-reference-genomes.md` are both in, plus `Jeilongvirus.2.dna` (`19c94c7`,
not from the tier list) — **87 reference files**, 500-genome validation at 363/363 among calls
clearing `-mcb 150`, 78% of the set accepted (85% weighted). Tier 3 is unapplied and is not
failure-driven. The larger
remaining gap is not references at all: Dianlovirus, Thamnovirus, Oblavirus, Striavirus and
Parajeilongvirus are **missing taxa**, with no PSSM set or JSON entry, so adding a reference contig
for them would route records to a taxon that cannot annotate them.

## Family-scoped classification (`-vfam`)

For the 7.2% of the failure set whose GenBank lineage stops above genus, `-vtax` is unusable —
there is no genus to declare. `-vfam` restricts the classification sweep to the reference contigs
of one family's taxa and lets BLASTn pick among them.

**The arithmetic that governs when it helps.** Restricting the reference set can only *lower* the
maximum bit score, so `best_within_family ≤ best_overall`: a genome rejected by `-mcb` unrestricted
is still rejected under `-vfam`. Restriction rescues nothing on its own, and it was a mistake to
expect it to. What it changes is the *reliability of a call in the band below the cutoff* — the
validated rescue numbers in `LOWVAN-FAILURE-ANALYSIS.md` ("Rescuing the genus-less 7.2% by BLAST")
give ~93% accuracy for a family-restricted call in the 50–150 bit band against roughly 60% for an
unrestricted one, and the genus-less B1 subset has median best hit 114 (max 149, zero at ≥150), so
that band is exactly where these records live. The intended use is a second pass at `-vfam` with a
lowered `-mcb`, not `-vfam` alone.

**Measured over all 322 genus-less genomes of the 5,000 sample, 2026-08-20** (raw:
`failure-analysis/results-322-vfam150.tsv` / `-vfam50.tsv`, analysis in `rerun-5000-analysis.md`
§6A):

| arm | annotated | job exits 0 |
|---|---:|---:|
| pass 1, unrestricted | 28.0% | 28.0% |
| `-vfam -mcb 150` | 28.0% | **75.5%** |
| `-vfam -mcb 50` | **36.6%** | **92.9%** |

- **The control arm is where the value is.** At `-mcb 150` the pass-1 → `-vfam` status transition is
  the identity on all 322 rows and no genome changed genus, exactly as the arithmetic predicts — but
  job success still rose 28.0% → 75.5%, because the sidecar sets `viral_family` and guard 2 stops
  firing. 153 jobs recovered with no new classification calls, hence no accuracy risk. Adopt this.
- **`-mcb 50` is on hold.** It adds 28 annotated genomes, but of the 58 newly-classified only 9 have
  an organism name that decides the genus, and those go 6 right / 2 wrong / 1 out of scope. Do not
  reuse the ~93% band-validation figure here: that validation was drawn from records whose lineage
  resolves to a LowVan taxon, and this population is enriched for divergent and unrepresented
  viruses (two "Langya virus" records called Jeilongvirus; a letovirus called Alphacoronavirus).
- **Guard 1 counts all features, not LowVan features.** A genome with zero LowVan features still
  exits 0 when the input GTO carried GenBank features. Only `gb_feats = 0` records die there.
- 90/90 agreement with pass 1's unrestricted genus on the rows pass 1 annotated, at both cutoffs;
  0 of 322 chose a different genus at 50 than at 150.

Implementation notes:

- `-vtax` semantics are untouched. Exact taxon match is tried first, so the eight Bunyavirales names
  keep meaning what they meant; only a name that is *not* a taxon but *is* a family falls through to
  `-vfam`, with a note on stderr. That ordering is load-bearing — don't reverse it.
- `-mcb` is deliberately reused as the floor rather than adding a second knob, and the below-floor
  path exits 1 exactly as the undeclared path does. `-vtax`'s "downgrade the rejection to a warning"
  exception does **not** extend to `-vfam`.
- The `<prefix>.classification` sidecar is written under `-vfam` as well as `-margin`, with three
  extra lines (`scope family`, `family <name>`, `min_contig_bit <n>`). This is what lets the GTO
  wrapper recover the chosen genus when the annotation comes back empty: `viral_family` falls back
  to the sidecar's `annotated_as`, so guard 2 passes and the downstream stages run.
- The within-family margin is printed whether or not `-margin` is given, and is more meaningful than
  the unrestricted one — the runner-up is a sibling genus, not an unrelated family.
- The GTO wrapper gained `--min-contig-bit` (passed through as `-mcb`) at the same time; it had
  never surfaced the contig-bitscore gate, which made `--viral-family` useless from that path.

## Lineage fallback (`--fallback-viral-taxon`, GTO wrapper only)

`-vtax` and `-vfam` both decide up front, for the whole run, whether to trust BLASTn.
`--fallback-viral-taxon` decides per genome and only after the fact: run the normal classification,
and if it finds no reference above `-mcb`, retry once with `-vtax` derived from the lineage in the
**input GTO's own `taxonomy` field**. BLASTn keeps first refusal, so a genome it can place is
annotated byte-identically to a run without the flag. Off by default; the regression on
`A-California-07-2009-H1N1.gto` is identical apart from per-run timestamps and event UUIDs.

Implementation notes, each of which is a decision that could reasonably have gone the other way:

- **The trigger is the stderr text, not the exit code** — `classification_rejected()` matches the two
  strings the base script prints (`No matching reference contigs`, `No reference contig of family`).
  rc cannot tell a graceful rejection from a BLAST crash, and retrying a crash under a declared taxon
  would crash again with a provenance record claiming the taxon was chosen on purpose. Keep those two
  strings in step with `annotate_by_viral_pssm.pl` if the wording ever changes.
- **The valid taxa are asked of the base script** (`-list-vtax -pssm … -j …`), not re-derived, so a
  `LOWVAN_DATA_DIR` with a different taxon set answers for itself. Any failure returns the empty list
  and the wrapper simply does not fall back — the safe direction.
- **The lineage is walked most-specific-first, exact-string, spaces → underscores.** All three
  matter: `Orthopneumovirus` is a strict prefix of `Orthopneumovirus_muris` (so no prefix matching),
  the binomial is the only taxon with a space, and a Bunyavirales lineage must fall *past* genus and
  species to land on the family. Tested against all four cases plus a lineage stopping at family.
- **Refused with `--viral-taxon`** (nothing to fall back from; `-mcb` is only a warning there),
  **allowed with `--viral-family`** — that one can still fail `-mcb`, and the retry drops `-vfam`
  because a genus is strictly more specific than a family scope.
- `-vfam` is deliberately *not* used as the fallback: restriction can only lower the best bit score,
  so it cannot rescue a genome that already failed `-mcb` unrestricted.
- The retry reassigns `@params` before the analysis event is built, so the event records the
  invocation that produced the output. The first attempt's stderr is moved to
  `<prefix>.stderr-classification.txt` — the retry would otherwise destroy the only record of why the
  classification was rejected. `close_genomes` gets `classification_scope = lineage_fallback`,
  `classification_taxon`, `classification_lineage`; `viral_family` falls back to the declared taxon
  (plain `-vtax` writes no sidecar, so the wrapper has to record this itself).

## Margin scoring (`-margin`)

`-mcb` asks whether the best hit is strong enough; `-margin` asks whether it beat the runner-up.
The two come apart, and a near-tie in the wrong direction is a mis-annotation rather than a failed
job — see the evidence above and `failure-analysis/margin-scoring.md`.

- **Off by default and deliberately un-defaulted.** The 1.3 that separates the historical bad call
  is fitted to one counter-example, and `19c94c7` removed that counter-example. Don't promote it to
  a default without a second case.
- **At 5,000 genomes there is no second case, and no signal.** Over the 3,747 calls with a lineage
  truth, nine of the ten wrong calls sit at margins 1.52–2.14 and one at 19.75, while correct calls
  start at 1.04 (p1 = 1.70). A threshold of 1.5 catches **zero** of the ten; 2.25 catches nine but
  flags 117 genomes to do it, 108 of them correct. Margin does not separate right calls from wrong
  ones on this population — keep it as a per-genome diagnostic, don't treat it as a classifier, and
  don't spend a rerun testing a new threshold: the column is recorded for every genome, so any
  threshold can be evaluated after the fact from `failure-analysis/results-5000-pass1.tsv`.
- It **never rejects.** It warns on stderr and records. For the failure set a flagged annotation is
  worth more than another dead job.
- The runner-up is free: `%taxon_best` already holds each taxon's best score from the sweep. That
  is also why `-margin` is refused with `-skip-classification` — one taxon searched, no runner-up.
- `<prefix>.classification` is the sidecar: tab-separated key/value lines
  (`taxon`, `reference`, `bit`, `runner_up`, `runner_up_bit`, `margin`, `margin_threshold`,
  `below_threshold`, `annotated_as`) followed by repeated `score\t<taxon>\t<bit>\t<file>` ranking
  lines. `margin` is the literal string `inf` when nothing else scored above zero.
- The inner script has `chdir`ed into the temp dir by the time it writes the sidecar, so the path is
  rooted at `$base` (the original cwd) unless `-p` is absolute. The GTO wrapper reads it from `$here`
  and merges `margin` / `runner_up` / `runner_up_value` / `margin_below_threshold` onto the
  `close_genomes` record; with the flag off there is no file and the record keeps its original six
  keys. `read_classification` must stay silent when the file is absent — that is the normal case.

## Failure analysis

`failure-analysis/` holds the investigation of the ~160k-genome reannotation failure set, with
scripts and raw results:

- `../LOWVAN-FAILURE-ANALYSIS.md` — root cause, the failure buckets, and the recovery plan
- `rerun-5000-analysis.md` — **the full analysis of the 5,000-genome rerun**: sample design and
  weights, per-taxon breakdowns, both failure buckets in full, the ten misclassifications, threats
  to validity, operational findings, and the column format of every `results-5000-*.tsv`. Start here
  before cutting a new slice from the raw results
- `classification-failure-mechanism.md` — why the BLASTn step picks the wrong genus
- `margin-scoring.md` — using the winner/runner-up bit-score ratio to tell a real call from a tie
- `proposed-reference-genomes.md` — specific references to add, with expected impact (Tiers 1 and 2
  are applied; Tier 3 is not)
- `validation-set-500.tsv`, `validation-set.tsv`, `draw_stratified_sample.py`,
  `build_ref_sweep.py`, `sweep_one.sh`, `gb2fa.py` — the classification regression harness
  described above
- `results-500-refs78.tsv`, `results-500-refs86.tsv` — the raw before/after sweeps behind the
  500-genome numbers
- `job-log-scan-5000.tsv` — the production stdout/stderr scan that confirmed the rc cascade at
  scale, from `/home/olson/P3/dev-ubuntu/modules/bvbrc_lowvan/Viral_Annotation/jobs_failed_output/`

### Rerunning the pipeline in bulk

The `pipeline_*` scripts run the **whole four-stage pipeline** over a sample, not just
classification. `pipeline_env.sh` is the entry point for all of them; `pipeline_driver.sh` (plain
`--min 1`) and `pipeline_driver_vtax.sh` (`--viral-taxon` on whatever the first pass failed) are
resumable and default to `P=48`. `pipeline_summarize.py` / `pipeline_summarize_vtax.py` regenerate
every number in `LOWVAN-FAILURE-ANALYSIS.md`. Raw output for the 5,000-genome run of 2026-08-19/20
is in `results-5000-pass1.tsv`, `results-5000-pass2-vtax.tsv`, `results-5000-b1-subthreshold.tsv`
and `results-5000-b2-maxbit.tsv`.

`pipeline_driver_vfam.sh` / `pipeline_run_vfam.sh` / `pipeline_summarize_vfam.py` are the third
pass — family-scoped classification over the 322 genus-less rows, one `-mcb` per run
(`MCB=50 bash failure-analysis/pipeline_driver_vfam.sh`), defaulting to `P=32`. Raw output:
`results-322-vfam150.tsv`, `results-322-vfam50.tsv`, plus
`vfam-newly-classified-organisms.tsv` (the 58 newly-classified genomes with their GenBank organism
names, which is the only accuracy evidence this subset admits).

**Copy the tree and the GenBank inputs to `/tmp` first.** `/home/olson` and `/home/mshukla` are
NFS, one genome costs ~410 `tblastn` invocations each opening a PSSM file, and throughput at P=48
was *lower* than at P=20 with the 96 CPUs 96% idle — the filer, not the CPU, is the resource.
`pipeline_env.sh` encodes the tmpfs paths and the rsync commands; read its header before changing
`REPO`. The repo must also stay first on `PATH`: the P3 env ships a deployed
`annotate_by_viral_pssm.pl` predating `-vtax`.

The same trap applies to anything else that fans out over the PSSM set — process startup (0.142 s)
dwarfs the alignment (0.005 s), so one BLAST call per (genome, PSSM) pair is the wrong shape. Batch
all genomes of a taxon into one database and issue one call per PSSM instead, as
`pipeline_b2_maxbit.sh` does. If you do, raise `-evalue`: it scales with database size and the
default silently clips marginal HSPs a single-genome search would report.

### What the 5,000-genome rerun established (2026-08-19/20)

Baseline is 0% — every genome in the sample is a production failure.

| Lever | Effect | |
|---|---|---|
| `--min 1` | **67.5%** annotated | The gate excludes 94.1% of the failure set; nothing else comes close |
| `-vtax` from the GenBank lineage | +14.2 pts annotated, **+24.7 pts job success** | Recovers 76.7% of classification failures and **0%** of the other bucket, by construction |
| `-mcb 50` | 41.3% of the residual B1 at 92.7% accuracy | Accuracy is flat 89–94% from 125 down to 40, then collapses. Not for the normal production path |
| Lowering `bit_cutoff` | ~nothing | 606 of 620 "classified but no features" genomes have no HSP above 100 bits at all; median best is 24 |
| `-margin` | nothing | See above |

Combined: **81.7% annotated, 92.2% exit 0.** The 522 jobs that succeed with an empty annotation are
the point of `-vtax`, not a defect — for a 200 bp fragment "annotated nothing" is the right answer,
and it is what turns a rc=255 three stages downstream into a legible result. Everything still
failing after both passes dies at guard 1 (`No features in GTO`), which is now the only remaining
structural cause.

Two data defects surfaced, both recorded in `UPSTREAM-ISSUES.md` with the measured impact: the
post-stop-codon coverage recheck (fires 33 times in 20 of 5,000 genomes, all genuinely
low-coverage — the evidence issue #4 asked for before keeping the correction), and the orphaned
`Orthorubulavirus` `P-ORF2`/`I-ORF2` PSSM directories.

## Dependencies

- Perl 5.38+ with: JSON::XS, File::Slurp, IPC::Run, Getopt::Long, Getopt::Long::Descriptive
- `gjoseqlib.pm` from https://github.com/TheSEED/seed_gjo/
- BV-BRC modules: GenomeTypeObject.pm, P3DataAPI
- BLAST+ 2.13.0 (blastn, tblastn)—JSON output format is version-specific
- MMSeqs2, MAFFT (for PSSM generation)

For BV-BRC internal users: `source /vol/patric3/cli/ubuntu-cli/user-env.sh`

## Environment Variable

`LOWVAN_DATA_DIR`: Override default data directory paths (default: `/home/jjdavis/bin/Viral_Annotation`)
